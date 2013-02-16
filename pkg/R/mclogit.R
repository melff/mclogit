# require(methods)
# if(!require(Matrix)) warning("mclogit with random effects won't work without Matrix package")

quickInteraction <- function(by){
  if(is.list(by)){
    n.arg <- length(by)
    f <- 0L
    uf <- 0L
    for(i in rev(1:n.arg)){
      y <- by[[i]]
      y <- as.numeric(y)
      uy <- unique(na.omit(y))
      y <- match(y,uy,NA)
      l <- length(uy)
      f <- f*l + y - 1
      uf <- unique(na.omit(f))
      f <- match(f,uf,NA)
      uf <- seq(length(uf))
    }
  }
  else {
    by <- as.numeric(by)
    uf <- unique(na.omit(by))
    f <- match(by,uf,NA)
    uf <- seq(length(uf))
  }
  return(structure(f,unique=uf))
}

constInSets <- function(X,sets){
    ans <- integer(0)
    for(i in 1:ncol(X)){
        v <- tapply(X[,i],sets,var)
        if(all(v[is.finite(v)]==0)) ans <- c(ans,i)
    }
    names(ans) <- colnames(X)[ans]
    ans
}

mclogit <- function(
                formula,
                data=parent.frame(),
                random=NULL,
                subset,
                weights,
                offset=NULL,
                na.action = getOption("na.action"),
                model = TRUE, x = FALSE, y = TRUE,
                contrasts=NULL,
                start.theta=NULL,
                control=mclogit.control(...),
                ...
            ){
# Assumptions:
#   left hand side of formula: cbind(counts,  choice set index)
#   right hand side of the formula: attributes
#   intercepts are removed!

    call <- match.call(expand.dots = TRUE)

    if(length(random))
        random <- setupRandomFormula(random)

    if(missing(data)) data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "offset", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    na.action <- attr(mf,"na.action")
    weights <- as.vector(model.weights(mf))
    offset <- as.vector(model.offset(mf))
    if(!is.null(weights) && !is.numeric(weights))
        stop("'weights' must be a numeric vector")
    Y <- as.matrix(model.response(mf, "any"))
    if(ncol(Y)<2) stop("need response counts and choice set indicators")
    sets <- Y[,2]
    sets <- match(sets,unique(sets))
    Y <- Y[,1]
    if (is.null(weights)){
        prior.weights <- rep(1,length(Y))
        N <- rowsum(Y,sets,na.rm=TRUE)
        weights <- N[sets]
        }
    else{
        prior.weigths <- weights
        N <- rowsum(weights*Y,sets,na.rm=TRUE)
        weights <- N[sets]
        }
    N <- sum(N)
    good <- weights > 0
    w <- weights[good]
    Y <- Y[good]/weights[good]
    s <- sets[good]
    s <- match(s,unique(s))
    X <- model.matrix(mt,mf,contrasts)
    contrasts <- attr(X, "contrasts")
    xlevels <- .getXlevels(mt,mf)
    X <- X[good,,drop=FALSE]
    icpt <- match("(Intercept)",colnames(X),nomatch=0)
    if(icpt) X <- X[,-icpt,drop=FALSE]
    const <- constInSets(X,s)
    if(length(const)){
        warning("removing ",paste(names(const),collapse=",")," from model")
        X <- X[,-const,drop=FALSE]
    }
    fit <- mclogit.fit(Y,s,w,X,
                        control=control,
                        offset = offset)
    null.dev <- fit$null.deviance
    if(length(random)){ ## random effects
        if(length(all.vars(random$covariates))){
          warning("random slopes not yet implemented, ignoring covariates")
          random$covarates <- ~1
        }
        mfr <- match.call(expand.dots = FALSE)
        mr <- match(c("formula", "data", "subset", "weights", "na.action"), names(mfr), 0)
        mfr <- mfr[c(1, mr)]
        environment(random$all.vars) <-environment(formula)
        mfr$formula <- random$all.vars
        mfr$drop.unused.levels <- TRUE
        mfr[[1]] <- as.name("model.frame")
        mfr <- eval(mfr, parent.frame())
        if(!require(Matrix))
          stop("need Matrix package installed")
        Z <- reDesignMatrix(random,mfr,use=good)
        fit <- mclogit.fit.random(Y,s,w,X,Z,
                              start=fit$coef,
                              start.theta=start.theta,
                              control=control)
     }
    if(x) fit$x <- X
    if(x && length(random)) fit$z <- Z
    if(!y) {
        fit$y <- NULL
        ftt$s <- NULL
    }
    fit$null.deviance <- null.dev
    fit <- c(fit,list(call = call, formula = formula,
                      terms = mt,
                      random = random,
                      data = data,
                      contrasts = contrasts,
                      xlevels = xlevels,
                      na.action = na.action,
                      prior.weights=prior.weights,
                      weights=weights,
                      model=mf,
                      N=N))
    if(length(random))
        class(fit) <- c("mclogitRandeff","mclogit","lm")
    else
        class(fit) <- c("mclogit","lm")
    fit
}

mclogit.fit <- function(
      Y,
      s,
      w,
      X,
      start=NULL,
      offset=NULL,
      control=mclogit.control()
      ){
    nvar <- ncol(X)
    deviance <- Inf
    sw <- c(tapply(w,s,"[",1))
    sqrt.sw <- sqrt(sw)
    X.qr <- qr(X)
    Q <- qr.Q(X.qr)
    R <- qr.R(X.qr)
    p <- X.qr$rank
    pivot <- X.qr$pivot
    qr.select <- abs(diag(R)) > 1e-7
    Q <- Q[,qr.select,drop=FALSE]
    R <- R[qr.select,qr.select,drop=FALSE]
    pivot.sel <- pivot[qr.select]
    eta <- mclogitLinkInv(Y,s,w)
    nobs <- length(eta)
    if (is.null(offset))
        offset <- rep.int(0, nobs)

    P <- mclogitP(eta,s)
    if(length(start))
      last.coef1 <- R%*%start
    else last.coef1 <- NULL
    converged <- FALSE
    for(iter in 1:control$maxit){
        UPQ <- rowsum(Q*P,s)*sqrt.sw
        s.Pw <- sqrt(P*w)
        QWQ <- crossprod(Q*s.Pw) - crossprod(UPQ)
        y <- P*(eta - offset) + Y-P
        Py <-  w*y
        UPy <- rowsum(y,s)*sqrt.sw
        QWy <- crossprod(Q,Py) - crossprod(UPQ,UPy)
        coef1 <- try(solve(QWQ,QWy),silent=TRUE)
        if(inherits(coef1,"try-error"))
            coef1 <- ginv(QWQ)%*%QWy
        eta <- c(Q%*%coef1) + offset
        P <- mclogitP(eta,s)
        last.deviance <- deviance
        dev.resids <- mclogit.dev.resids(Y,P,w)
        deviance <- sum(dev.resids)
          ## check for divergence
          boundary <- FALSE
          if(!is.finite(deviance)){
            if(is.null(last.coef1))
                stop("no valid set of coefficients has been found: please supply starting values", call. = FALSE)
             warning("step size truncated due to divergence", call. = FALSE)
             ii <- 1
             while (!is.finite(deviance)){
                if(ii > control$maxit)
                  stop("inner loop; cannot correct step size")
                ii <- ii + 1
                coef1 <- (coef1 + last.coef1)/2
                eta <- c(Q %*% coef1) + offset
                P <- mclogitP(eta,s)
                dev.resids <- mclogit.dev.resids(Y,P,w)
                deviance <- sum(dev.resids)
             }
             boundary <- TRUE
             if (control$trace)
                  cat("Step halved: new deviance =", deviance, "\n")
          } ## inner loop
        if(control$trace)
            cat("\nIteration",iter,"- Deviance =",deviance)
        crit <- abs(deviance-last.deviance)/abs(0.1+deviance)
        if(crit < control$eps){
          converged <- TRUE
          cat("\nconverged\n")
          break
        }
    }
    if (!converged) warning("algorithm did not converge")
    if (boundary) warning("algorithm stopped at boundary value")
    eps <- 10*.Machine$double.eps
    if (any(P < eps) || any(1-P < eps))
        warning("fitted rates numerically 0 occurred")
    coef <- rep(NA,nvar)
    coef[pivot.sel] <- solve(R,coef1)
    covmat <- matrix(NA,nrow=nvar,ncol=nvar)
    UPQ <- rowsum(Q*P,s)*sqrt.sw
    s.Pw <- sqrt(P*w)
    QWQ <- crossprod(Q*s.Pw) - crossprod(UPQ)
    iR <- solve(R)
    covmat[pivot.sel,pivot.sel] <- iR %*% solve(QWQ,t(iR))
    names(coef) <- colnames(X)
    colnames(covmat) <- rownames(covmat) <- names(coef)
   ntot <- length(Y)
   P0 <- mclogitP(offset,s)
   null.deviance <- sum(mclogit.dev.resids(Y,P0,w))
   resid.df <- length(Y)#-length(unique(s))
   model.df <- ncol(X)
   resid.df <- resid.df - model.df
   phi <- sum(w/P*(Y-P)^2)/resid.df
   ll <- mclogit.logLik(Y,P,w)
   return(list(
      coefficients = drop(coef),
      linear.predictors = eta,
      residuals = (Y-P)/P,
      residual.df = resid.df,
      model.df = model.df,
      fitted.values = P,
      deviance=deviance,
      ll=ll,
      deviance.residuals=dev.resids,
      null.deviance=null.deviance,
      iter = iter,
      y = Y,
      s = s,
      offset = offset,
      converged = converged,
      control=control,
      covmat=covmat,
      phi=phi
      ))
}

tr <- function(x) sum(diag(x))




mclogit.fit.random <- function(
      Y,
      s,
      w,
      X,
      Z,
      start=NULL,
      start.theta=NULL,
      control=mclogit.control()
      ){
    if(!require(Matrix))
      stop("need Matrix package installed")
    #crossprod <- Matrix:::crossprod
    nvar <- ncol(X)
    deviance <- Inf
    lev.ics <- attr(Z,"col.indices")
    nlev <- length(lev.ics)
    sw <- c(tapply(w,s,"[",1))
    sqrt.sw <- sqrt(sw)
#    g1 <- lapply(groupings,function(g)lapply(g,function(g)g[good]))
    nvalid <- length(s)
#    nlev <- length(g1)
    X.qr <- qr(X)
    Q <- qr.Q(X.qr)
    R <- qr.R(X.qr)
    rank <- X.qr$rank
    pivot <- X.qr$pivot
    qr.select <- abs(diag(R)) > 1e-7
    Q <- Q[,qr.select,drop=FALSE]
    colnames(Q) <- paste("Q",seq(ncol(Q)),sep="")
    #R <- R[qr.select,qr.select,drop=FALSE]
    pivot.sel <- pivot[qr.select]
    coef1 <- R%*%start[pivot.sel]
    eta <- c(Q%*%coef1)
    P <- mclogitP(eta,s)
    s.Pw <- sqrt(P*w)
    dev.resids <- ifelse(Y>0,
            2*w*Y*(log(Y)-log(P)),
            0)
    deviance <- sum(dev.resids)
    PU <- Matrix(0,nrow=length(s),ncol=length(unique(s)))
    PU[cbind(seq(nvalid),s)] <- P
    converged <- FALSE
    lev.seq <- seq(length(lev.ics))

    ## Score Test for Variance Parameters
      UPr <- rowsum(Y-P,s)*sqrt.sw
      UPZ <- crossprod(PU,Z)*sqrt.sw
      u <- crossprod(Z,w*(Y - P)) #- crossprod(UPZ,UPr)
      A <- crossprod(Z*s.Pw) - crossprod(UPZ)
      l <- numeric(nlev)
      K <- matrix(nrow=nlev,ncol=nlev)
      for(i in lev.seq){
        ii <- lev.ics[[i]]
        A.i <- A[ii,ii]
        l[i] <- crossprod(u[ii]) - tr(A.i)
        K[i,i] <- sum(A.i*A.i)
        for(j in lev.seq[lev.seq > i]){
          jj <- lev.ics[[j]]
          A.ij <- A[ii,jj]
          K[i,j] <- K[j,i] <- sum(A.ij*A.ij)
        }
      }
    chisq.theta <- l^2/diag(K)
   ## Starting values for variance parameters
    if(length(start.theta)){
      if(length(start.theta)<length(lev.ics))
        stop("insufficient starting values for variance parameters")
      if(length(start.theta)>length(lev.ics))
        stop("to many starting values for variance parameters")
      last.theta <- theta <- start.theta
    }
    else{
      if(all(l < 0)) stop("insufficient residual variance; set some variance parameters to zero")
      theta <- solve(K,l)
      if(any(theta<0)){
        warning("Negative initial estimates found; correcting",call.=FALSE)
        cat("\ntheta=",theta)
        l[l < 0] <- 0
        cat("\nl=",l)
        ii <- 0
        bb <- 1
        I <- diag(x=diag(K))
        while(any(theta<0)){
          ii <- ii + 1
          cat("\nTrial ",ii)
          bb <- bb/2
          aa <- 1 - bb
#           if(ii > control$maxit)
#               stop("insufficient residual variance; set some variance parameters to zero")
          theta <- solve(aa*I+bb*K,l)
          cat("\ntheta=",theta)
        }
      }
      last.theta <- theta
    }
    if(control$trace)
      cat("\nInitial estimate of theta: ",theta)
    b <- rep(0,ncol(Z))
    rept <- sapply(lev.ics,length)
    #Sigma <- as(Diagonal(x=rep(theta,rept)),"sparseMatrix")
    ## Extended IWLS and Fisher-scoring for variance parameters
    converged <- FALSE
    for(iter in 1:control$maxit){
        ## Updating coef1
        UPQ <- rowsum(Q*P,s)*sqrt.sw
        s.Pw <- sqrt(P*w)
        QWQ <- crossprod(Q*s.Pw) - crossprod(UPQ)
        Py <-  w*(P*eta + Y-P)
        UPy <- rowsum(P*eta+Y-P,s)*sqrt.sw
        QWy <- crossprod(Q,Py) - crossprod(UPQ,UPy)
        PU[cbind(seq(nvalid),s)] <- P
        UPZ <- crossprod(PU,Z)*sqrt.sw
#         browser()
        A <- crossprod(Z*s.Pw) - crossprod(UPZ)
        diag(A) <- diag(A) + rep(1/theta,rept)
        A <- solve(A)
        ZWQ <- crossprod(Z,w*P*Q)
        # coercion necessary for some Matrix package weirdness ...
        ZWQ <- as(ZWQ,"sparseMatrix") - as(crossprod(UPZ,UPQ),"sparseMatrix")
        B <- crossprod(ZWQ,A)
        H <- QWQ - B%*%ZWQ
        ZWy <- crossprod(Z,Py)@x - crossprod(UPZ,UPy)@x # These should be vectors anyway
        h <- QWy - B%*%ZWy
        last.coef1 <- coef1
        coef1 <- solve(H,h)
        coef1 <- c(as(coef1,"matrix"))
        delta.coef1 <- coef1 - last.coef1
        last.b <- b
        b <- A%*%(ZWy-ZWQ%*%coef1)
        eta <- c(as(Q%*%coef1 + Z%*%b,"matrix"))
        P <- mclogitP(eta,s)
        dev.resids <- mclogit.dev.resids(Y,P,w)
        last.deviance <- deviance
        deviance <- sum(dev.resids)
          ## check for divergence
          boundary <- FALSE
          #if(!is.finite(deviance)){
          if(deviance > last.deviance){
            if(is.null(last.coef1) || is.null(last.theta))
                stop("no valid set of coefficients has been found: please supply starting values", call. = FALSE)
             if(!is.finite(deviance)) warning("step size truncated due to divergence", call. = FALSE)
             ii <- 1
             while (!is.finite(deviance) || deviance > last.deviance){
                if(ii > control$maxit){
                  if(is.finite(deviance)){
                    warning("cannot decrease deviance",call.=FALSE)
                    break
                  } else
                      stop("inner loop; cannot correct step size")
                  }
                ii <- ii + 1
                coef1 <- (coef1 + last.coef1)/2
                b <- (b + last.b)/2
                theta <- (theta + last.theta)/2
                eta <- c(as(Q%*%coef1 + Z%*%b,"matrix"))
                P <- mclogitP(eta,s)
                dev.resids <- mclogit.dev.resids(Y,P,w)
                deviance <- sum(dev.resids)
                if (control$trace){
                      cat("Step halved: new deviance =", deviance, "\n")
                      }
                if(is.finite(deviance) && deviance > last.deviance){
                  crit <- abs(deviance-last.deviance)/abs(0.1+deviance)
                  if(crit < control$eps) {
                    #warning("cannot decrease deviance",call.=FALSE)
                    break
                  }
                }
             }
             #boundary <- TRUE
          } ## inner loop

        ## Fisher-scoring step for variance parameters
        last.theta <- theta
        grad.theta <- numeric(nlev)
        theta. <- numeric(nlev)
        Info.theta <- matrix(0,ncol=nlev,nrow=nlev)
        for(i in lev.seq){
          ii <- lev.ics[[i]]
          m.i <- length(ii)
          A.i <- A[ii,ii]
          grad.theta[i] <- (crossprod(b[ii]) -theta[i]*m.i + tr(A.i))/theta[i]^2
          Info.theta[i,i] <- m.i/theta[i]^2 - 2*tr(A.i)/theta[i]^3 + sum((A.i@x)^2)/theta[i]^4
          for(j in lev.seq[lev.seq > i]){
            jj <- lev.ics[[j]]
            A.ij <- A[ii,jj]
            Info.theta[i,j] <- Info.theta[j,i] <- sum((A.ij)^2)/theta[i]^2/theta[j]^2
          }
        }
        delta.theta <- solve(Info.theta,grad.theta)
        theta <- theta + delta.theta
        if(any(theta<0)){
          ## Handle negative variances
          warning("negative values of variance parameters occured",call.=FALSE)
          ii <- 1
          while(any(theta<0)){
            if(ii > control$maxit)
                  stop("inner loop; cannot correct step size")
            if(control$trace)
              cat("\nSecond order inner iteration ",ii," - theta = ",theta)
            ii <- ii + 1
            theta <- (theta + last.theta)/2
          }
        }
        if(control$trace)
            cat("\nIteration",iter,"- Deviance =",deviance,#"theta =",theta,
            "criterion = ",abs(deviance-last.deviance)/abs(0.1+deviance)#,
            #"criterion[2] = ",max(abs(theta - last.theta))
            )

        ## Checking convergence
        #crit <- (control$eps)^(-2) *max(abs(delta.theta/(theta+.1)),abs(delta.coef1/(coef1+.1)))
        #crit.theta <- max(abs(delta.theta/(theta+.1)))
        crit <- abs(deviance-last.deviance)/abs(0.1+deviance)
        if(crit < control$eps){
          converged <- TRUE
          if(control$trace) cat("\nconverged\n")
          break
        }
    }
    if (!converged) warning("algorithm did not converge")
    if (boundary) warning("algorithm stopped at boundary value")
    eps <- 10*.Machine$double.eps
    if (any(P < eps) || any(1-P < eps))
        warning("fitted rates numerically 0 occurred")
    coef <- rep(NA,nvar)
    coef[pivot.sel] <- solve(R,coef1)
    covmat <- matrix(NA,nrow=nvar,ncol=nvar)
    UPQ <- rowsum(Q*P,s)*sqrt.sw
    s.Pw <- sqrt(P*w)
    QWQ <- crossprod(Q*s.Pw) - crossprod(UPQ)
    PU[cbind(seq(nvalid),s)] <- P
    UPZ <- crossprod(PU,Z)*sqrt.sw
    A <- crossprod(Z*s.Pw) - crossprod(UPZ)
    diag(A) <- diag(A) + rep(1/theta,rept)
    A <- solve(A)
    ZWQ <- crossprod(Z,w*P*Q)
    # coercion necessary for some Matrix package weirdness ...
    ZWQ <- as(ZWQ,"sparseMatrix") - as(crossprod(UPZ,UPQ),"sparseMatrix")
    B <- crossprod(ZWQ,A)
    H <- QWQ - B%*%ZWQ
    iR <- solve(R)
    covmat[pivot.sel,pivot.sel] <- as(iR %*% solve(H,t(iR)),"matrix")
    names(coef) <- colnames(X)
    colnames(covmat) <- rownames(covmat) <- names(coef)
    H.theta <- matrix(ncol=nlev,nrow=nlev)
    for(i in lev.seq){
      ii <- lev.ics[[i]]
      m.i <- length(ii)
      A.i <- A[ii,ii]
      Info.theta[i,i] <- m.i/theta[i]^2 - 2*tr(A.i)/theta[i]^3 + sum((A.i@x)^2)/theta[i]^4
      for(j in lev.seq[lev.seq > i]){
        jj <- lev.ics[[j]]
        A.ij <- A[ii,jj]
        Info.theta[i,j] <- Info.theta[j,i] <- sum((A.ij)^2)/theta[i]^2/theta[j]^2
      }
    }
   itheta2 <- ifelse(theta==0,0,1/theta^2)
   H.theta <- sweep(K,2,itheta2,"*")*itheta2
   H.theta.eigen <- eigen(H.theta)
   if(any(H.theta.eigen$values <= 0)){
    U.Ht <- H.theta.eigen$vectors
    id.Ht <- ifelse(H.theta.eigen$values <= 0,0,1/H.theta.eigen$values)
    covmat.theta <- U.Ht%*%(id.Ht*t(U.Ht))
    }
   else
    covmat.theta <- solve(H.theta)
   names(theta) <- names(lev.ics)
   colnames(covmat.theta) <- rownames(covmat.theta) <- names(theta)
   reff <- b
   reff <- lapply(lev.ics,function(ii)reff[ii])
   B <- as.matrix(B)
   resid.df <- length(Y)#-length(unique(s))
   model.df <- ncol(X) + length(theta)
   resid.df <- resid.df-model.df
   phi <- sum(w/P*(Y-P)^2)/resid.df
   return(list(
      coefficients = drop(coef),
      varPar = theta,
      chisq.theta = chisq.theta,
      random.effects = reff,
      linear.predictors = eta,
      fitted.values = P,
      ll=NA,
      deviance=deviance,
      deviance.residuals=dev.resids,
      residuals=(Y-P)/P,
      residual.df = resid.df,
      model.df = model.df,
      iter = iter,
      y = Y,
      s = s,
      converged = converged,
      control=control,
      covmat=covmat,
      covmat.varPar = covmat.theta,
      rank=rank,
      phi=phi
      ))
}


setupRandomFormula <- function(formula){
  fo <- delete.response(terms(formula))
  attributes(fo) <- NULL
  if(length(fo[[2]]) < 2 || as.character(fo[[2]][1])!="|")
    stop("missing '|' operator")
  covariates <- fo
  groups <- fo
  covariates[2] <- fo[[2]][2]
  groups[2] <- fo[[2]][3]
  list(
    covariates=covariates,
    groups=groups,
    all.vars=reformulate(all.vars(fo))
    )
}

reDesignMatrix <- function(random,data,use=NULL){
  covariates <- all.vars(random$covariates)
  if(length(covariates)) warning("covariates not yet implemented")
  if(!length(use)) use <- TRUE
  groups <- data[use,all.vars(random$groups),drop=FALSE]
  gnames <- names(groups)
  n <- length(groups[[1]])
  nlev <- length(groups)
  groups[[1]] <- quickInteraction(groups[[1]])
  if(nlev>1)
    for(i in 2:nlev)
      groups[[i]] <- quickInteraction(groups[c(i-1,i)])
  un <- length(attr(groups[[1]],"unique"))
  Z <- Matrix(0,nrow=n,ncol=un,dimnames=list(NULL,paste(gnames[1],seq(un),sep="")))
  ij <- cbind(1:n,groups[[1]])
  Z[ij] <- 1
  lev.ics <- list()
  lev.ics[[1]] <- seq.int(un)
  if(nlev>1)
    for(i in 2:nlev){
      un.i <- length(attr(groups[[i]],"unique"))
      Z.i <- Matrix(0,nrow=n,ncol=un.i,dimnames=list(NULL,paste(gnames[i],seq(un.i),sep="")))
      ij <- cbind(1:n,groups[[i]])
      Z.i[ij] <- 1
      Z <- cbind2(Z,Z.i)
      lev.ics.i <- seq.int(un.i) + max(lev.ics[[i-1]])
      lev.ics <- c(lev.ics,list(lev.ics.i))
    }
  rand.names <- all.vars(random$groups)
  if(length(rand.names) > 1)
    for(i in 2:length(rand.names))
      rand.names[i] <- paste(rand.names[i-1],rand.names[i],sep=":")
  names(lev.ics) <- rand.names
  structure(Z,col.indices=lev.ics)
}

mclogit.control <- function(
                            epsilon = 1e-08,
                            maxit = 25,
                            trace=TRUE
                            ) {
    if (!is.numeric(epsilon) || epsilon <= 0)
        stop("value of epsilon must be > 0")
    if (!is.numeric(maxit) || maxit <= 0)
        stop("maximum number of iterations must be > 0")
    list(epsilon = epsilon, maxit = maxit, trace = trace)
}

log.Det <- function(x) determinant(x,logarithm=TRUE)$modulus

mclogitP <- function(eta,s){
  expeta <- exp(eta)
  sum.expeta <- rowsum(expeta,s)
  expeta/sum.expeta[s]
}

mclogit.dev.resids <- function(y,p,w)
      ifelse(y>0,
                2*w*y*(log(y)-log(p)),
                0)

mclogit.logLik <- function(y,p,w) sum(w*y*log(p))
                
                
mclogitLinkInv <- function(y,s,w){
  n.alt <- tapply(y,s,length)
  c(log(sqrt(w)*y+1/n.alt[s])-log(w)/2)
}

print.mclogit <- function(x,digits= max(3, getOption("digits") - 3), ...){
    cat("\nCall: ", deparse(x$call), "\n\n")
    if(length(coef(x))) {
        cat("Coefficients")
        if(is.character(co <- x$contrasts))
            cat("  [contrasts: ",
                apply(cbind(names(co),co), 1, paste, collapse="="), "]")
        cat(":\n")
        print.default(format(x$coefficients, digits=digits),
                      print.gap = 2, quote = FALSE)
    } else cat("No coefficients\n\n")
    if(length(x$varPar)) {
        cat("Variance paremeters")
        cat(":\n")
        print.default(format(x$varPar, digits=digits),
                      print.gap = 2, quote = FALSE)
    }
    if(nchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep="")
    cat("\nDispersion:       ",   format(signif(x$phi,digits)),
        "\nNull Deviance:    ",   format(signif(x$null.deviance, digits)),
        "\nResidual Deviance:", format(signif(x$deviance, digits)),
        "\n")
    invisible(x)
}

vcov.mclogit <- function(object,...){
    return(object$covmat)
}

weights.mclogit <- function(object,...){
  return(object$weights)
}

deviance.mclogit <- function(object,...){
  return(object$deviance)
}

summary.mclogit <- function(object,dispersion=NULL,correlation = FALSE, symbolic.cor = FALSE,...){

    if(!length(dispersion))
        dispersion <- max(object$phi,1)


    ## calculate coef table

    coef <- object$coefficients
    covmat.scaled <- object$covmat * dispersion
    var.cf <- diag(covmat.scaled)
    s.err <- sqrt(var.cf)
    zvalue <- coef/s.err
    pvalue <- 2*pnorm(-abs(zvalue))

    coef.table <- array(NA,dim=c(length(coef),4))
    dimnames(coef.table) <- list(names(coef),
            c("Estimate", "Std. Error","z value","Pr(>|z|)"))
    coef.table[,1] <- coef
    coef.table[,2] <- s.err
    coef.table[,3] <- zvalue
    coef.table[,4] <- pvalue

    if(length(object$varPar)){
      covmat.vp.scaled <- object$covmat.varPar * dispersion
      varPar <- object$varPar
      var.vp <- diag(covmat.vp.scaled)
      s.err.vp <- sqrt(var.vp)
      #zvalue.vp <- varPar/s.err.vp
      #pvalue.vp <- 2*pnorm(-abs(zvalue.vp))
      zvalue.vp <- sqrt(object$chisq.theta)
      pvalue.vp <- pchisq(object$chisq.theta,1,lower.tail=FALSE)


      varPar.table <- array(NA,dim=c(length(varPar),4))
      dimnames(varPar.table) <- list(names(varPar),
            c("Estimate", "Std. Error","z value","Pr(>|z|)"))
      varPar.table[,1] <- varPar
      varPar.table[,2] <- s.err.vp
      varPar.table[,3] <- zvalue.vp
      varPar.table[,4] <- pvalue.vp
    } else varPar.table <- NULL

    ans <- c(object[c("call","terms","baseline","deviance","contrasts",
                       "null.deviance","iter","na.action","phi","model.df","residual.df")],
              list(coefficients = coef.table,
                    varPar = varPar.table,
                    cov.scaled=covmat.scaled,
                    cov.unscaled=object$covmat,
                    cov.varPar.scaled=object$covmat.varPar * dispersion,
                    cov.varPar.unscaled=object$covmat.varPar))
    p <- length(coef)
    if(correlation && p > 0) {
        dd <- sqrt(diag(ans$covmat.unscaled))
        ans$correlation <-
            ans$covmat.unscaled/outer(dd,dd)
        ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.mclogit"
    return(ans)
}

print.summary.mclogit <-
    function (x, digits = max(3, getOption("digits") - 3),
              symbolic.cor = x$symbolic.cor,
              signif.stars = getOption("show.signif.stars"), ...){
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep="\n", collapse="\n"), "\n\n", sep="")

    coefs <- x$coefficients
    if(length(x$varPar)){
      varPar <- x$varPar
      rownames(varPar) <- paste("Var(",rownames(varPar),")",sep="")
      coefs <- rbind(coefs,varPar)
    }
    printCoefmat(coefs, digits=digits, signif.stars=signif.stars,
                     na.print="NA", ...)
#     cat("\n")
#     cat("AIC: ", format(x$aic, digits= max(4, digits+1)),"\n\n",
#         "Number of Fisher Scoring iterations: ", x$iter,
#         "\n", sep="")
    cat("\nDispersion:       ",   format(signif(x$phi,digits)), "on", x$residual.df, "degrees of freedom",
        "\nNull Deviance:    ",   format(signif(x$null.deviance, digits)),
        "\nResidual Deviance:", format(signif(x$deviance, digits)),
        "\nNumber of Fisher Scoring iterations: ", x$iter,
        "\n")
    correl <- x$correlation
    if(!is.null(correl)) {
        p <- NCOL(correl)
        if(p > 1) {
            cat("\nCorrelation of Coefficients:\n")
            if(is.logical(symbolic.cor) && symbolic.cor) {
                print(symnum(correl, abbr.colnames = NULL))
            } else {
                correl <- format(round(correl, 2), nsmall = 2, digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -p, drop=FALSE], quote = FALSE)
            }
        }
    }

    cat("\n")
    invisible(x)
}



getSummary.mclogit <- function(obj,
            alpha=.05,
            rearrange=NULL,
            #as.columns=NULL,
            ...){

  smry <- summary(obj)
  N <- obj$N
  coef <- smry$coefficients
  varPar <- smry$varPar

  lower.cf <- qnorm(p=alpha/2,mean=coef[,1],sd=coef[,2])
  upper.cf <- qnorm(p=1-alpha/2,mean=coef[,1],sd=coef[,2])
  coef <- cbind(coef,lower.cf,upper.cf)
  colnames(coef) <- c("est","se","stat","p","lwr","upr")
  if(length(varPar)){
    se.log.varPar <- varPar[,1]*varPar[,2]
    lower.log.varPar <- qnorm(p=alpha/2,mean=log(varPar[,1]),sd=se.log.varPar[2])
    upper.log.varPar <- qnorm(p=1-alpha/2,mean=log(varPar[,1]),sd=se.log.varPar[2])
    varPar <- cbind(varPar,exp(lower.log.varPar),exp(upper.log.varPar))
    colnames(varPar) <- c("est","se","stat","p","lwr","upr")
    rownames(varPar) <- paste("Var(",rownames(varPar),")",sep="")
  }
  if(length(rearrange)){
      coef.grps <- lapply(rearrange,function(ii){
          if(is.character(ii) && any(ii %nin% rownames(coef)))
               stop("coefficient(s) ",dQuote(unname(ii[ii %nin% rownames(coef)]))," do not exist")
          structure(coef[ii,],
              dimnames=list(names(ii),dimnames(coef)[[2]])
              )
          })
      grp.titles <- names(rearrange)
      coef.grps <- memisc:::clct.arrays(coef.grps)
      coef <- array(NA,dim=c(
          dim(coef.grps)[1] + NROW(varPar),
          dim(coef.grps)[2],
          dim(coef.grps)[3]
          ))
      coef[seq(dim(coef.grps)[1]),,] <- coef.grps
      if(length(varPar))
       coef[dim(coef.grps)[1]+seq(nrow(varPar)),,1] <- varPar
      dimnames(coef) <- list(
         c(dimnames(coef.grps)[[1]],rownames(varPar)),
         dimnames(coef.grps)[[2]],
         grp.titles
         )
    }
#     else if(length(as.columns)){
#       #groupix <- sapply(as.columns,grep,rownames(coef),simplify=FALSE)
#       groupix <- sapply(as.columns,function(patrn) if(is.atomic(patrn))
#                                        which(memisc:::str.has(rownames(coef),patrn))
#                                       else
#                                        which(do.call(memisc:::str.has,
#                                            c(list(rownames(coef)),patrn)
#                                            )),
#                           simplify=FALSE
#                           )
#       grp.titles <- if(length(names(as.columns))) names(as.columns)
#           else if(is.atomic(as.columns)) as.columns
#           else sapply(as.columns,function(patrn)
#                                   if(is.atomic(patrn))
#                                       paste(patrn,collapse=":")
#                                   else paste(patrn[[1]],collapse=":")
#                               )
#       ii <- sapply(groupix,length)>0
#       groupix <- groupix[ii]
#       as.columns <- as.columns[ii]
#       grouped <- unlist(groupix)
#       ungrouped <- setdiff(seq(nrow(coef)),grouped)
#       max.glen <- max(sapply(groupix,length))
#       n.grp <- length(groupix)
#       coef.groups <- lapply(groupix,function(i)coef[i,,drop=FALSE])
#       for(g in 1:n.grp){
#           newnames <- rownames(coef.groups[[g]])
#           newnames <- strsplit(newnames,":")
#           patrn <- as.columns[[g]]
#           newnames <- lapply(newnames,function(x){
#                 if(is.atomic(patrn))
#                   hasit <- memisc:::str.has(x,patrn)
#                 else
#                   hasit <- memisc:::str.has(x,patrn[[1]],not=patrn$not)
#                 x <- x[!hasit]
#                 if(length(x)) paste(x,collapse=":") else "Main effect"
#               })
#           rownames(coef.groups[[g]]) <- newnames
#         }
#       coef.groups <- memisc:::clct.arrays(coef.groups)
#       coef.ungrouped <- coef[ungrouped,]
#       coef <- array(NA,dim=c(
#         dim(coef.groups)[1] + nrow(coef.ungrouped)+NROW(varPar),
#         dim(coef.groups)[2],
#         dim(coef.groups)[3]
#         ))
#       coef[seq(dim(coef.groups)[1]),,] <- coef.groups
#       if(nrow(coef.ungrouped))
#         coef[dim(coef.groups)[1]+seq(nrow(coef.ungrouped)),,1] <- coef.ungrouped
#       if(length(varPar))
#         coef[dim(coef.groups)[1]+nrow(coef.ungrouped)+seq(nrow(varPar)),,1] <- varPar
#       dimnames(coef) <- list(
#           c(dimnames(coef.groups)[[1]],rownames(coef.ungrouped),rownames(varPar)),
#           dimnames(coef.groups)[[2]],
#           grp.titles
#           )
#       }
    else {
      .coef <- coef
      coef <- matrix(NA,nrow=nrow(.coef)+NROW(varPar),ncol=ncol(.coef))
      coef[seq(nrow(.coef)),] <- .coef
      if(length(varPar))
        coef[nrow(.coef)+seq(nrow(varPar)),] <- varPar
      rownames(coef) <- c(rownames(.coef),rownames(varPar))
      colnames(coef) <- colnames(.coef)
    }


   phi <- smry$phi
   LR <- smry$null.deviance - smry$deviance
#   df <- smry$df.null - smry$df.residual
  df <- obj$model.df
  deviance <- deviance(obj)


  if(df > 0){
    p <- pchisq(LR,df,lower.tail=FALSE)
    L0.pwr <- exp(-smry$null.deviance/N)
    LM.pwr <- exp(-smry$deviance/N)

    McFadden <- 1- smry$deviance/smry$null.deviance
    Cox.Snell <- 1 - exp(-LR/N)
    Nagelkerke <- Cox.Snell/(1-L0.pwr)
    }
  else {
    LR <- NA
    df <- NA
    p <- NA
    McFadden <- NA
    Cox.Snell <- NA
    Nagelkerke <- NA
    }

  ll <- obj$ll
  AIC <- AIC(obj)
  BIC <- AIC(obj,k=log(N))
  sumstat <- c(
          phi         = phi,
          LR             = LR,
          df         = df,
          #p             = p,
          logLik        = ll,
          deviance      = deviance,
          McFadden      = McFadden,
          Cox.Snell       = Cox.Snell,
          Nagelkerke    = Nagelkerke,
          AIC           = AIC,
          BIC           = BIC,
          N             = N
          )

  #coef <- apply(coef,1,applyTemplate,template=coef.template)

  #sumstat <- drop(applyTemplate(sumstat,template=sumstat.template))
  list(coef=coef,sumstat=sumstat)
}



fitted.mclogit <- function(object,type=c("probabilities","counts"),...){
  weights <- object$weights
  nobs <- length(weights)
  res <- rep(NA,nobs)
  res[weights>0] <- object$fitted.values
  type <- match.arg(type)
  switch(type,
    probabilities=res,
    counts=weights*res)
}

predict.mclogit <- function(object, newdata,type=c("link","response"),se=FALSE,...){

  type <- match.arg(type)
  rhs <- object$formula[-2]
  m <- model.frame(rhs,data=newdata)
  X <- model.matrix(rhs,m,
          contasts.arg=object$contrasts,
          xlev=object$xlevels
          )
  drop <- match("(Intercept)",colnames(X))
  X <- X[,-drop]
  eta <- c(X %*% object$coef)
  if(se){
    stopifnot(dim(X) == dim(X %*% vcov(object)))
    se.eta <- sqrt(rowSums(X * (X %*% vcov(object))))
  }

  if(type=="response") {
    lhs <- object$formula[[2]]
    set <- lhs[[3]]
    set <- eval(set,newdata,parent.frame())
    set <- match(set,unique(set))
    exp.eta <- exp(eta)
    sum.exp.eta <- rowsum(exp.eta,set)
    p <- exp.eta/sum.exp.eta[set]
    if(se)
      list(pred=p,se.pred=p*(1-p)*se.eta)
    else p
  }
  else if(se) list(pred=eta,se.pred=se.eta) else eta
}

logLik.mclogit <- function(object,...){
    if (length(list(...)))
        warning("extra arguments discarded")
    val <- if(length(object$ll))
            object$ll
           else NA
    attr(val, "nobs") <- object$N
    attr(val, "df") <- object$model.df
    class(val) <- "logLik"
    return(val)
}