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
                weights=NULL,
                offset=NULL,
                na.action = getOption("na.action"),
                model = TRUE, x = FALSE, y = TRUE,
                contrasts=NULL,
                start=NULL,
                control=mclogit.control(...),
                ...
            ){
# Assumptions:
#   left hand side of formula: cbind(counts,  choice set index)
#   right hand side of the formula: attributes
#   intercepts are removed!

    call <- match.call(expand.dots = TRUE)

    if(missing(data)) data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "offset", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")

    if(length(random)){
        rf <- paste(c(".~.",all.vars(random)),collapse="+")
        rf <- as.formula(rf)
        mff <- structure(mf$formula,class="formula")
        mf$formula <- update(mff,rf)
    }
    
    mf <- eval(mf, parent.frame())
    mt <- terms(formula)
    
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
        prior.weights <- weights
        N <- rowsum(weights*Y,sets,na.rm=TRUE)
        weights <- N[sets]
        }
    N <- sum(N)
    Y <- Y/weights
    Y[weights==0] <- 0
    
    X <- model.matrix(mt,mf,contrasts)
    contrasts <- attr(X, "contrasts")
    xlevels <- .getXlevels(mt,mf)
    icpt <- match("(Intercept)",colnames(X),nomatch=0)
    if(icpt) X <- X[,-icpt,drop=FALSE]
    const <- constInSets(X,sets)
    if(length(const)){
        warning("removing ",paste(names(const),collapse=",")," from model")
        X <- X[,-const,drop=FALSE]
    }
    if(!length(start)){
      drop.coefs <- check.mclogit.drop.coefs(Y,sets,weights,X,
                                             offset = offset)
      if(any(drop.coefs)){
        warning("removing ",paste(colnames(X)[drop.coefs],collapse=",")," from model")
        X <- X[,!drop.coefs,drop=FALSE]
      }
    }
    
    
    fit <- mclogit.fit(Y,sets,weights,X,
                        control=control,
                        start = start,
                        offset = offset)
    null.dev <- fit$null.deviance
    if(length(random)){ ## random effects
        
        random <- setupRandomFormula(random)
        rt <- terms(random$formula)
        
        groups <- random$groups
        rX <- model.matrix(rt,mf,contrasts)
        
        groups <- mf[groups]
        
        nlev <- length(groups)
        groups[[1]] <- quickInteraction(groups[1])

        if(nlev > 1){
            for(i in 2:nlev)
                groups[[i]] <- quickInteraction(groups[c(i-1,i)])
        }
        
        Z <- lapply(groups,mkZ,rX=rX)
        G <- mkG(rX)
        G <- rep(list(G),length(groups))
        names(Z) <- names(groups)
        names(G) <- names(groups)
        
        fit <- mmclogit.fitPQL(Y,sets,weights,
                               X,Z,G,groups,
                               start=fit$coef,
                               control=control,
                               offset = offset)
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
        class(fit) <- c("mmclogit","mclogit","lm")
    else
        class(fit) <- c("mclogit","lm")
    fit
}


check.mclogit.drop.coefs <- function(y,
                                     s,
                                     w,
                                     X,
                                     offset){
  nvar <- ncol(X)
  nobs <- length(y)
  if(!length(offset))
    offset <- rep.int(0, nobs)
  eta <- mclogitLinkInv(y,s,w)
  pi <- mclogitP(eta,s)
  y.star <- eta - offset + (y-pi)/pi
  yP.star <- y.star - rowsum(pi*y.star,s)[s]
  XP <- X - rowsum(pi*X,s)[s,,drop=FALSE]
  ww <- w*pi
  good <- ww > 0
  wlsFit <- lm.wfit(x=XP[good,,drop=FALSE],y=yP.star[good],w=ww[good])  
  is.na(wlsFit$coef)
}

mclogit.fit <- function(
      y,
      s,
      w,
      X,
      start=NULL,
      offset=NULL,
      control=mclogit.control()
      ){
    nvar <- ncol(X)
    nobs <- length(y)
    if(!length(offset))
      offset <- rep.int(0, nobs)
    if(length(start)){
      stopifnot(length(start)==ncol(X))
      eta <- c(X%*%start) + offset
    }
    else
      eta <- mclogitLinkInv(y,s,w)
    pi <- mclogitP(eta,s)
    dev.resids <- ifelse(y>0,
                         2*w*y*(log(y)-log(pi)),
                         0)
    deviance <- sum(dev.resids)
    if(length(start))
      last.coef <- start
    else last.coef <- NULL
    converged <- FALSE
    for(iter in 1:control$maxit){
        y.star <- eta - offset + (y-pi)/pi
        yP.star <- y.star - rowsum(pi*y.star,s)[s]
        XP <- X - as.matrix(rowsum(pi*X,s))[s,,drop=FALSE]
        ww <- w*pi
        good <- ww > 0
        wlsFit <- lm.wfit(x=XP[good,,drop=FALSE],y=yP.star[good],w=ww[good])
        coef <- wlsFit$coefficients
        
        eta <- c(X%*%coef) + offset
        pi <- mclogitP(eta,s)
        last.deviance <- deviance
        dev.resids <- ifelse(y>0,
                2*w*y*(log(y)-log(pi)),
                0)
        deviance <- sum(dev.resids)
          ## check for divergence
          boundary <- FALSE
          if(!is.finite(deviance)){
            if(is.null(last.coef))
                stop("no valid set of coefficients has been found: please supply starting values", call. = FALSE)
             warning("step size truncated due to divergence", call. = FALSE)
             ii <- 1
             while (!is.finite(deviance)){
                if(ii > control$maxit)
                  stop("inner loop; cannot correct step size")
                ii <- ii + 1
                coef <- (coef + last.coef)/2
                eta <- c(X %*% coef) + offset
                pi <- mclogitP(eta,s)
                dev.resids <- ifelse(y>0,2*w*y*(log(y)-log(pi)),0)
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
          if(control$trace)
            cat("\nconverged\n")
          break
        }
    }
    if (!converged) warning("algorithm did not converge")
    if (boundary) warning("algorithm stopped at boundary value")
    #eps <- 10*.Machine$double.eps
    #if (any(pi < eps) || any(1-pi < eps))
    #    warning("fitted probabilities numerically 0 occurred")

    XP <- X - as.matrix(rowsum(pi*X,s))[s,,drop=FALSE]
    ww <- w*pi
    Information <- crossprod(XP,ww*XP)
        
    ntot <- length(y)
    pi0 <- mclogitP(offset,s)
    null.deviance <- sum(ifelse(y>0,
                    2*w*y*(log(y)-log(pi0)),
                    0))
    resid.df <- length(y)#-length(unique(s))
    model.df <- ncol(X)
    resid.df <- resid.df - model.df
    ll <- mclogit.logLik(y,pi,w)
    return(list(
        coefficients = drop(coef),
        linear.predictors = eta,
        working.residuals = (y-pi)/pi,
        working.weights = w,
        response.residuals = y-pi,
        residual.df = resid.df,
        model.df = model.df,
        fitted.values = pi,
        deviance=deviance,
        ll=ll,
        deviance.residuals=dev.resids,
        null.deviance=null.deviance,
        iter = iter,
        y = y,
        s = s,
        offset = offset,
        converged = converged,
        control=control,
        covmat=solve(Information)
        ))
}


mmclogit.fitPQL <- function(
                            y,
                            s,
                            w,
                            X,
                            Z,
                            G,
                            groups,
                            start,
                            offset=NULL,
                            control=mclogit.control()
      ){

    nvar <- ncol(X)
    nobs <- length(y)
    nsets <- length(unique(s))
    
    nlevs <- length(G)
    nvpar <- sapply(G,length)
    vpar.selector <- rep(1:nlevs,nvpar)

    if(!length(offset))
      offset <- rep.int(0, nobs)
    
    coef <- start
    eta <- c(X%*%coef) + offset
    pi <- mclogitP(eta,s)
    
    dev.resids <- ifelse(y>0,2*w*y*(log(y)-log(pi)),0)
    deviance <- sum(dev.resids)
    converged <- FALSE

    #
    mk <- lapply(groups,lunq)
    I.mk <- lapply(mk,Diagonal)
    G.star <- Map(G.star1,I.mk,G)

    ## Starting values for variance parameters

    sqrt.w <- sqrt(w)

    i <- 1:nobs
    W <- Matrix(0,nrow=nobs,ncol=nsets)
    W[cbind(i,s)] <- sqrt.w*pi
    W <- Diagonal(x=w*pi)-tcrossprod(W)
    
    y.star <- eta + (y-pi)/pi

    Wy <- W%*%y.star
    WX <- W%*%X
    XWX <- crossprod(X,WX)
    iXWX <- solve(XWX)

    XWy <- crossprod(X,Wy)

    WZ <- lapply(Z,`%*%`,x=W)
    ZWX <- lapply(Z,crossprod,y=WX)
    ZWy <- lapply(Z,crossprod,y=Wy)
    dim(ZWX) <- c(nlevs,1)
    dim(ZWy) <- c(nlevs,1)
    
    u <- list()
    S <- matrix(list(),nlevs,nlevs)

    for(k in 1:nlevs){
        ZWX.iXWX <- ZWX[[k]]%*%iXWX
        v.k <- ZWy[[k]] - ZWX.iXWX%*%XWy
        A.kk <- tcrossprod(ZWX.iXWX,ZWX[[k]])

        G.star.k <- G.star[[k]]
        u1 <- sapply(G.star.k,quadform,x=v.k)
        u2 <- sapply(G.star.k,tr.crossprod,B=A.kk)
        u[[k]] <- u1 - u2

        G.star.A.kk <- lapply(G.star.k,`%*%`,y=A.kk)
        
        hh.k <- length(G[[k]])
        S.kk <- matrix(0,hh.k,hh.k)

        lwrTri <- lower.tri(S.kk,diag=TRUE)
        h1 <- row(S.kk)[lwrTri]
        h2 <- col(S.kk)[lwrTri]
        tmp <- mapply(tr.crossprod,G.star.A.kk[h1],G.star.A.kk[h2])
        S.kk[lwrTri] <- tmp
        S.kk <- lwr2sym(S.kk)
        S[[k,k]] <- S.kk
        
        if(k<nlevs){
            for(r in (k+1):nlevs){

                A.kr <- tcrossprod(ZWX.iXWX,ZWX[[r]])
                G.star.A.kr <- lapply(G.star.k,`%*%`,y=A.kr)

                p.r <- ncol(Z[[r]])
                S.kr <- matrix(0,p.k,p.r)
                h1 <- rows(S.kr)
                h2 <- cols(S.kr)
                S.kr[] <- mapply(tr.crossprod,G.star.A.kr[h1],G.star.A.kr[h2])
                S[[k,r]] <- S.kr
                S[[k,r]] <- t(S.kr)                
            }
        }
    }
    u <- unlist(u)
    S <- fuseMat(S)
    
    theta <- solve(S,u)
    theta <- split(theta,vpar.selector)
    Phi <- list()
    for(k in 1:nlevs)
        Phi[[k]] <- fillG(G[[k]],theta[[k]])

    iPhi <- lapply(Phi,solve)

    Sigma <- Map(`%x%`,I.mk,Phi)
    Sigma <- do.call(bdiag,unname(Sigma))

    iSigma <- Map(`%x%`,I.mk,iPhi)
    iSigma <- do.call(bdiag,unname(iSigma))

    b.split <- rep(1:nlevs,sapply(Z,ncol))
    
    ## Extended IWLS and Fisher-scoring for variance parameters
    converged <- FALSE
    ZWZ <- matrix(list(),nlevs,nlevs)
    for(k in 1:nlevs){
        Z.k <- Z[[k]]
        ZWZ[[k,k]] <- crossprod(Z.k,W%*%Z.k)
        if(k<nlevs){
            for(r in (k+1):nlevs){
                Z.r <- Z[[r]]                    
                ZWZ.kr <- crossprod(Z.k,W%*%Z.r)
                ZWZ[[k,r]] <- ZWZ.kr
                ZWZ[[r,k]] <- t(ZWZ.kr)
            }
        }
    }
    
    ZWX <- fuseMat(ZWX)
    ZWy <- fuseMat(ZWy)
    ZWZ.iSigma <- fuseMat(ZWZ)+iSigma
    
    for(iter in 1:control$maxit){
        
        K <- solve(ZWZ.iSigma)
        XiVX <- XWX - crossprod(ZWX,K%*%ZWX)
        XiVy <- XWy - crossprod(ZWX,K%*%ZWy)
        coef <- solve(XiVX,XiVy)

        b <- K%*%(ZWy-ZWX%*%coef)

        b. <- split(b,b.split)
        Zb <- Map(`%*%`,Z,b.)
        Zb <- do.call(rbind,Zb)

        eta <- as.vector(X%*%coef) + as.vector(Zb) + offset
        pi <- mclogitP(eta,s)

        dev.resids <- ifelse(y>0,2*w*y*(log(y)-log(pi)),0)

        W <- Matrix(0,nrow=nobs,ncol=nsets)
        W[cbind(i,s)] <- sqrt.w*pi
        W <- Diagonal(x=w*pi)-tcrossprod(W)

        y.star <- eta + (y-pi)/pi

        Wy <- W%*%y.star
        WX <- W%*%X
        XWX <- crossprod(X,WX)
        iXWX <- solve(XWX)

        XWy <- crossprod(X,Wy)

        WZ <- lapply(Z,`%*%`,x=W)
        ZWX <- lapply(Z,crossprod,y=WX)
        ZWy <- lapply(Z,crossprod,y=Wy)
        dim(ZWX) <- c(nlevs,1)
        dim(ZWy) <- c(nlevs,1)

        for(k in 1:nlevs){
            Z.k <- Z[[k]]
            ZWZ[[k,k]] <- crossprod(Z.k,W%*%Z.k)
            if(k<nlevs){
                for(r in (k+1):nlevs){
                    Z.r <- Z[[r]]                    
                    ZWZ.kr <- crossprod(Z.k,W%*%Z.r)
                    ZWZ[[k,r]] <- ZWZ.kr
                    ZWZ[[r,k]] <- t(ZWZ.kr)
                }
            }
        }

        ## Fisher-scoring step for variance parameters
        ZWZ.iSigma <- fuseMat(ZWZ)+iSigma
        ZWX <- fuseMat(ZWX)
        ZWy <- fuseMat(ZWy)
        K <- solve(ZWZ.iSigma)

        ZWZ. <- list()
        for(k in 1:nlevs){
            ZWZ.[[k]] <- do.call(rbind,ZWZ[,k])
        }
        
        u <- list()
        S <- matrix(list(),nlevs,nlevs)
        iSigma.b <- iSigma%*%b
        iSigma.b <- split(iSigma.b,b.split)
        
        for(k in 1:nlevs){

            G.star.k <- G.star[[k]]
            iSigma.b.k <- iSigma.b[[k]]
            u[[k]] <- sapply(G.star.k,quadform,x=iSigma.b.k)

            ZWZ.kk <- ZWZ[[k,k]]
            ZWZ.k <- ZWZ.[[k]]
            
            A.kk <- ZWZ.kk - crossprod(ZWZ.k,K%*%ZWZ.k)
            G.star.A.kk <- lapply(G.star.k,`%*%`,y=A.kk)
            
            hh.k <- length(G[[k]])
            S.kk <- matrix(0,hh.k,hh.k)

            lwrTri <- lower.tri(S.kk,diag=TRUE)
            h1 <- row(S.kk)[lwrTri]
            h2 <- col(S.kk)[lwrTri]
            tmp <- mapply(tr.crossprod,G.star.A.kk[h1],G.star.A.kk[h2])
            S.kk[lwrTri] <- tmp
            S.kk <- lwr2sym(S.kk)
            S[[k,k]] <- S.kk

            if(k<nlevs){

                for(r in (k+1):nlevs){
                    
                    ZWZ.kr <- ZWZ[[k,r]]
                    ZWZ.r <- ZWZ.[[r]]

                    A.kr <- ZWZ.kr - crossprod(ZWZ.k,K%*%ZWZ.r)
                    G.star.A.kr <- lapply(G.star.k,`%*%`,y=A.kr)

                    p.r <- ncol(Z[[r]])
                    S.kr <- matrix(0,p.k,p.r)
                    h1 <- rows(S.kr)
                    h2 <- cols(S.kr)
                    S.kr[] <- mapply(tr.crossprod,G.star.A.kr[h1],G.star.A.kr[h2])
                    S[[k,r]] <- S.kr
                    S[[k,r]] <- t(S.kr)
                }
            }
        }

        u <- unlist(u)
        Info.theta <- fuseMat(S)

        last.theta <- theta
        theta <- solve(Info.theta,u)
        theta <- split(theta,vpar.selector)
        Phi <- list()
        for(k in 1:nlevs)
            Phi[[k]] <- fillG(G[[k]],theta[[k]])

        logDetPhi <- lapply(Phi,log.Det)
        logDetSigma <- Map(`*`,logDetPhi,mk)
        logDetSigma <- sum(unlist(logDetSigma))
        
        iPhi <- lapply(Phi,solve)

        Sigma <- Map(`%x%`,I.mk,Phi)
        Sigma <- do.call(bdiag,unname(Sigma))

        iSigma <- Map(`%x%`,I.mk,iPhi)
        iSigma <- do.call(bdiag,unname(iSigma))

        ZWZ.iSigma <- fuseMat(ZWZ)+iSigma
        b.iSigma.b <- as.numeric(crossprod(b,iSigma%*%b))

        logDet.ZWZ.iSigma <- log.Det(ZWZ.iSigma)

        last.deviance <- deviance
        deviance <- sum(dev.resids) + logDetSigma + logDet.ZWZ.iSigma + b.iSigma.b

        crit <- abs(deviance-last.deviance)/abs(0.1+deviance)
        
        if(control$trace)
            cat("\nIteration",iter,"- Deviance =",deviance#,"theta =",theta,
                #"criterion = ",abs(deviance-last.deviance)/abs(0.1+deviance),
                #                        "criterion[2] = ",max(abs(unlist(theta) - unlist(last.theta)))
                )

        if(crit < control$eps){
            converged <- TRUE
            if(control$trace) cat("\nconverged\n")
            break
        }
    }
    if (!converged) warning("algorithm did not converge")

    K <- solve(ZWZ.iSigma)
    XiVW <- XWX - crossprod(ZWX,K%*%ZWX)
    covmat.coef <- solve(XiVX)
    covmat.coef <- as.matrix(covmat.coef)
    
    covmat.theta <- solve(Info.theta)
    
    coef <- drop(coef)
    colnames(covmat.coef) <- rownames(covmat.coef) <- names(covmat.coef)

    Phi <- structure(lapply(Phi,as.matrix),
                     names=names(G))
                     
    resid.df <- length(y)#-length(unique(s))
    model.df <- ncol(X) + length(theta)
    resid.df <- resid.df-model.df

    return(list(
        coefficients = coef,
        VarCov = Phi,
        varPar = theta,
        linear.predictors = eta,
        fitted.values = pi,
        ll=NA,
        deviance=deviance,
        deviance.residuals=dev.resids,
        working.residuals=(y-pi)/pi,
        response.residuals=y-pi,
        residual.df = resid.df,
        model.df = model.df,
        iter = iter,
        y = y,
        s = s,
        groups = groups,
        G = G,
        converged = converged,
        control=control,
        covmat=covmat.coef,
        covmat.varPar = covmat.theta,
        rank=rank
    ))
}


setupRandomFormula <- function(formula){

    trms <- terms(formula)
    fo <- delete.response(trms)
    
    attributes(fo) <- NULL
    if(length(fo[[2]]) < 2 || as.character(fo[[2]][1])!="|")
        stop("missing '|' operator")

    groups <- fo
    fo[2] <- fo[[2]][2]
    groups[2] <- groups[[2]][3]
    list(
        formula=structure(fo,class="formula"),
        groups=all.vars(groups)
    )
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

# mclogit.dev.resids <- function(y,p,w)
#       ifelse(y>0,
#                 2*w*y*(log(y)-log(p)),
#                 0)

mclogit.logLik <- function(y,p,w) sum(w*y*log(p))
                
                
mclogitLinkInv <- function(y,s,w){
  #n.alt <- tapply(y,s,length)
  #c(log(sqrt(w)*y+1/n.alt[s])-log(w)/2)
  n <- w*y+0.5
  f <- n/(rowsum(n,s)[s])
  log(f) - ave(log(f),s)
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
    cat("\nNull Deviance:    ",   format(signif(x$null.deviance, digits)),
        "\nResidual Deviance:", format(signif(x$deviance, digits)))
    if(nchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep="")
    else cat("\n")
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

    ## calculate coef table

    coef <- object$coefficients
    covmat.scaled <- object$covmat 
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

    ans <- c(object[c("call","terms","deviance","contrasts",
                       "null.deviance","iter","na.action","model.df","residual.df","N")],
              list(coefficients = coef.table,
                    cov.coef=object$covmat))
    p <- length(coef)
    if(correlation && p > 0) {
        dd <- sqrt(diag(ans$cov.coef))
        ans$correlation <-
            ans$cov.coef/outer(dd,dd)
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
    printCoefmat(coefs, digits=digits, signif.stars=signif.stars,
                     na.print="NA", ...)

    cat("\nNull Deviance:    ",   format(signif(x$null.deviance, digits)),
        "\nResidual Deviance:", format(signif(x$deviance, digits)),
        "\nNumber of Fisher Scoring iterations: ", x$iter,
        "\nNumber of observations: ",x$N,
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

    if(nchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n\n", sep="")
    else cat("\n\n")
    invisible(x)
}




fitted.mclogit <- function(object,type=c("probabilities","counts"),...){
  weights <- object$weights
  
  res <- object$fitted.values
  type <- match.arg(type)
  
  na.act <- object$na.action
  
  res <- switch(type,
          probabilities=res,
          counts=weights*res)
          
  if(is.null(na.act))
    res
  else 
    napredict(na.act,res)
}

predict.mclogit <- function(object, newdata=NULL,type=c("link","response"),se.fit=FALSE,...){

  type <- match.arg(type)
  rhs <- object$formula[-2]
  if(missing(newdata)){
    m <- model.frame(rhs,data=object$model)
    na.act <- object$na.action
  }
  else{
    m <- model.frame(rhs,data=newdata)
    na.act <- attr(m,"na.action")
  }
  X <- model.matrix(rhs,m,
          contasts.arg=object$contrasts,
          xlev=object$xlevels
          )
  drop <- match("(Intercept)",colnames(X))
  X <- X[,-drop,drop=FALSE]
  eta <- c(X %*% coef(object))
  if(se.fit){
    V <- vcov(object)
    stopifnot(ncol(X)==ncol(V))
  }

  
  if(type=="response") {
    lhs <- object$formula[[2]]
    set <- lhs[[3]]
    set <- eval(set,newdata,parent.frame())
    if(length(na.act))
        set <- set[-na.act]
    set <- match(set,unique(set))
    exp.eta <- exp(eta)
    sum.exp.eta <- rowsum(exp.eta,set)
    p <- exp.eta/sum.exp.eta[set]
    if(se.fit){
      wX <- p*(X - rowsum(p*X,set)[set,,drop=FALSE])
      se.p <- sqrt(rowSums(wX * (wX %*% V)))
      if(is.null(na.act))
        list(fit=p,se.fit=se.p)
      else
        list(fit=napredict(na.act,p),
             se.fit=napredict(na.act,se.p))
    }
    else {
      if(is.null(na.act)) p
      else napredict(na.act,p)
    }
  }
  else if(se.fit) {
    se.eta <- sqrt(rowSums(X * (X %*% V)))
    if(is.null(na.act))
      list(fit=eta,se.fit=se.eta) 
    else
      list(fit=napredict(na.act,eta),
           se.fit=napredict(na.act,se.eta))
    }
  else {
    if(is.null(na.act)) eta
    else napredict(na.act,eta)
  }
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

residuals.mclogit <- 
 function(object, type = c("deviance", "pearson", "working",
                           "response", "partial"), ...){
   type <- match.arg(type)
   
   resid <- switch(type,
                 deviance=mclogit.dev.resids(object),
                 pearson=stop("not yet implemented"),
                 working=object$working.residuals,
                 response=object$response.residuals,
                 partial=stop("not yet implemented")
                 )
   naresid(object$na.action,resid)
 }
 
mclogit.dev.resids <- function(obj){
 y <- obj$y
 s <- obj$s
 w <- obj$weights
 pi <- obj$fitted.values
 
 n <- w*y+0.5
 f <- n/(rowsum(n,s)[s])
 #sign(y-p)*sqrt(2*abs(log(f)-log(y)))
 r <- 2*(f*log(f/pi))
 r - ave(r,s)
}


nobs.mclogit <- function(object,...) object$N

extractAIC.mclogit <- function(fit, scale = 0, k = 2, ...) 
{
  N <- fit$N
  edf <- N - fit$residual.df
  aic <- AIC(fit)
  c(edf, aic + (k - 2) * edf)
}

weights.mclogit <- function(object, type = c("prior", "working"),...) {
  type <- match.arg(type)
  res <- if (type == "prior") 
    object$prior.weights
  else object$weights
  if (is.null(object$na.action)) 
    res
  else naresid(object$na.action, res)
}


tr <- function(x) sum(diag(x))

mkZ <- function(groups,rX){

    n <- length(groups)
    ug <- unique(groups)
    m <- length(ug)
    p <- ncol(rX)
    n.j <- tabulate(groups,m)
    
    Z <- Matrix(0,nrow=n,ncol=m*p)

    i <- 1:n
    k <- 1:p
    j <- groups
    
    i <- rep(i,p)
    jk <- rep((j-1)*p,p)+rep(k,each=n)
    i.jk <- cbind(i,jk)

    Z[i.jk] <- rX
    Z
}

mkG <- function(rX){

    p <- ncol(rX)
    nms <- colnames(rX)
    
    G <- matrix(0,p,p)
    ltT <- lower.tri(G,diag=TRUE)
    ltF <- lower.tri(G,diag=FALSE)

    n <- p*(p+1)/2
    m <- p*(p-1)/2
    
    diag(G) <- 1:p
    G[ltF] <- p + 1:m
    G <- lwr2sym(G)
    rownames(G) <- colnames(G) <- nms
    
    lapply(1:n,mkG1,G)
}

mkG1 <- function(i,G) Matrix(array(as.integer(i==G),
                                   dim=dim(G),
                                   dimnames=dimnames(G)
                                   ))

fillG <- function(G,theta){

    Phi <- Map(`*`,theta,G)

    if(length(Phi)>1){
        for(i in 2:length(Phi))
            Phi[[1]] <- Phi[[1]] + Phi[[i]]
    }
    Phi[[1]]
}


lunq <- function(x)length(attr(x,"unique"))

G.star1 <- function(I,G)Map(`%x%`,list(I),G)
quadform <- function(A,x) as.numeric(crossprod(x,A%*%x))
tr.crossprod <- function(A,B) sum(A*B)

lwr2sym <- function(X){

    lwrX <- lower.tri(X)
    x.lwr <- X[lwrX]
    Y <- t(X)
    Y[lwrX] <- x.lwr
    Y
}


fuseMat <- function(X){
    if(ncol(X)>1)
        X <- apply(X,1,cbindList)
    do.call(rbind,X)
}


print.mmclogit <- function(x,digits= max(3, getOption("digits") - 3), ...){
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

    cat("\n(Co-)Variances:\n")
    VarCov <- x$VarCov
    for(k in 1:length(VarCov)){
        cat("Grouping level:",names(VarCov)[k],"\n")
        VarCov.k <- VarCov[[k]]
        VarCov.k[] <- format(VarCov.k, digits=digits)
        VarCov.k[upper.tri(VarCov.k)] <- ""
        print.default(VarCov.k, print.gap = 2, quote = FALSE)
    }
    
    cat("\nNull Deviance:    ",   format(signif(x$null.deviance, digits)),
        "\nResidual Deviance:", format(signif(x$deviance, digits)))
    if(nchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep="")
    else cat("\n")
    invisible(x)
}

summary.mmclogit <- function(object,dispersion=NULL,correlation = FALSE, symbolic.cor = FALSE,...){

    ## calculate coef table

    coef <- object$coefficients
    covmat.scaled <- object$covmat 
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

    G <- object$G
    nlevs <- length(G)
    nvpar <- sapply(G,length)
    vpar.selector <- rep(1:nlevs,nvpar)
    
    se.varPar <- sqrt(diag(object$covmat.varPar))
    se.varPar <- split(se.varPar,vpar.selector)
    VarCov.table <- list()
    for(k in 1:nlevs){
        VarCov.k <- object$VarCov[[k]]
        se.VarCov.k <- as.matrix(fillG(G[[k]],se.varPar[[k]]))
        VarCov.table.k <- array(NA,dim=c(dim(VarCov.k),2))
        VarCov.table.k[,,1] <- VarCov.k
        VarCov.table.k[,,2] <- se.VarCov.k
        dimnames(VarCov.table.k)[1:2] <- dimnames(VarCov.k)
        dimnames(VarCov.table.k)[[3]] <- c("Estimate", "Std. Error")
        VarCov.table[[k]] <- VarCov.table.k
    }
    names(VarCov.table) <- names(object$VarCov)
    
    ans <- c(object[c("call","terms","deviance","contrasts",
                       "null.deviance","iter","na.action","model.df","residual.df","N")],
              list(coefficients = coef.table,
                   cov.coef=object$covmat,
                   VarCov = VarCov.table))
    p <- length(coef)
    if(correlation && p > 0) {
        dd <- sqrt(diag(ans$cov.coef))
        ans$correlation <-
            ans$cov.coef/outer(dd,dd)
        ans$symbolic.cor <- symbolic.cor
    }

 
    class(ans) <- "summary.mmclogit"
    return(ans)
}


print.summary.mmclogit <-
    function (x, digits = max(3, getOption("digits") - 3),
              symbolic.cor = x$symbolic.cor,
              signif.stars = getOption("show.signif.stars"), ...){
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep="\n", collapse="\n"), "\n\n", sep="")

    coefs <- x$coefficients
    cat("Coefficents:\n")
    printCoefmat(coefs, digits=digits, signif.stars=signif.stars,
                     na.print="NA", ...)

    cat("\n(Co-)Variances:\n")
    VarCov <- x$VarCov
    for(k in 1:length(VarCov)){
        cat("Grouping level:",names(VarCov)[k],"\n")
        VarCov.k <- VarCov[[k]]

        utri <- rep(upper.tri(VarCov.k[,,1]),2)
        VarCov.k[utri] <- NA
        VarCov.k <- ftable(VarCov.k,col.vars=3:2)
        VarCov.k <- format(VarCov.k, digits=digits, quote=FALSE)[-3,-2]
        VarCov.k[-(1:2),-1] <- gsub("NA","  ",VarCov.k[-(1:2),-1],fixed=TRUE)
        
        VarCov.k[1,] <- format(trimws(VarCov.k[1,]),justify="left")
        VarCov.k <- format(VarCov.k)
        cat(paste(apply(VarCov.k,1,paste,collapse=" "),collapse="\n"),"\n")
    }

    cat("\nNull Deviance:    ",   format(signif(x$null.deviance, digits)),
        "\nResidual Deviance:", format(signif(x$deviance, digits)),
        "\nNumber of Fisher Scoring iterations: ", x$iter,
        "\nNumber of observations: ",x$N,
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
    
    if(nchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n\n", sep="")
    else cat("\n\n")
    invisible(x)
}
