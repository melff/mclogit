mmclogit.fitPQLMQL <- function(
                               y,
                               s,
                               w,
                               X,
                               Z,
                               d,
                               start = NULL,
                               start.Phi = NULL,
                               start.b = NULL,
                               offset = NULL,
                               method = c("PQL","MQL"),
                               estimator = c("ML","REML"),
                               control=mmclogit.control()
                               ){
    method <- match.arg(method)
    estimator <- match.arg(estimator)

    nvar <- ncol(X)
    nobs <- length(y)
    nsets <- length(unique(s))
    nlevs <- length(Z)
    m <- sapply(Z,ncol)/d

    sqrt.w <- sqrt(w)

    i <- 1:nobs
    
    if(!length(offset))
      offset <- rep.int(0, nobs)
    if(length(start)){
      stopifnot(length(start)==ncol(X))
      eta <- c(X%*%start) + offset
      if(method=="PQL"){
          if(length(start.b) == nlevs){
              for(k in 1:nlevs)
                  eta <- eta +  as.vector(Z[[k]]%*%start.b[[k]])
          }
          else stop("PQL requires starting values for random effects")
      }
    }
    else
      eta <- mclogitLinkInv(y,s,w)
    pi <- mclogitP(eta,s)
    dev.resids <- ifelse(y>0,
                         2*w*y*(log(y)-log(pi)),
                         0)
    deviance <- sum(dev.resids)

    # Outer iterations: update non-linear part of the model
    converged <- FALSE
    fit <- NULL
    do.backup <- FALSE
    step.truncated <- FALSE

    msg <- "Random effects design matrix at index %d has fewer rows than columns (%d < %d).
This will almost certainly lead to noncovergence or other numerical problems.
Please reconsider your model specification."
    
    for(k in 1:nlevs){
        Z.k <- Z[[k]]
        if(nrow(Z.k) < ncol(Z.k))
           warning(sprintf(msg,k,nrow(Z.k),ncol(Z.k)))
    }
    parms <- NULL
    last.parms <- NULL
    last.deviance <- deviance
    prev.last.deviance <- NULL
    last.eta <- eta

    model.struct <- list(y=y,
                         s=s,
                         nsets=nsets,
                         nobs=nobs,
                         i=i,
                         w=w,
                         sqrt.w=sqrt.w,
                         offset=offset,
                         X=X,
                         Z=Z,
                         d=d,
                         m=m,
                         nlevs=nlevs)
    
    parms$coefficients <- list(fixed=start,
                               random=start.b)
    parms$Phi <- start.Phi
    for(iter in 1:control$maxit){

        W <- Matrix(0,nrow=nobs,ncol=nsets)
        W[cbind(i,s)] <- sqrt.w*pi
        W <- Diagonal(x=w*pi)-tcrossprod(W)

        y.star <- eta - offset + (y-pi)/pi

        # cat("\n")
        # print(head(y.star))
        
        prev.last.parms <- last.parms
        last.parms <- parms
        
        aux <- list(y=y.star,W=W)

        parms <- PQLMQL_innerFit(parms,aux,model.struct,method,estimator,control)

        step.back <- FALSE
        if(inherits(parms,"try-error")){
            if(length(prev.last.deviance) && 
               last.deviance > prev.last.deviance && 
               length(prev.last.parms)){
                # Previous step increased the deviance, so we better step back twice
                warning("Numeric problems in inner iteration and previous step increased deviance,
  stepping back twice")
                parms <- prev.last.parms
            }
            else { # Previous step decreased the deviance
                warning("Numeric problems in inner iteration, stepping back")
                parms <- last.parms
            }
            step.back <- TRUE
        }
        
        last.fit <- fit
        fit <- PQLMQL_eval_parms(parms,model.struct,method,estimator)
        
        deviance <- fit$deviance
        if(control$trace){
            cat("\nIteration",iter,"- deviance =",deviance)
        }

        if(is.finite(deviance)){
            if(deviance > last.deviance && control$break.on.increase){
                warning("Cannot decrease the deviance, stepping back",call.=FALSE)
                step.back <- TRUE
                parms <- last.parms
                fit <- last.fit
                deviance <- fit$deviance
            }
            if(deviance < 0 && control$break.on.negative){
                warning("Negative deviance, backing up",call.=FALSE)
                step.back <- TRUE
                parms <- last.parms
                fit <- last.fit
                deviance <- fit$deviance
            }
        }
        else if(!is.finite(deviance)){
            warning("Non-finite deviance, backing up",call.=FALSE)
            step.back <- TRUE
            parms <- last.parms
            fit <- last.fit
            deviance <- fit$deviance
            
        }

        eta <- fit$eta
        pi <- fit$pi
        coef <- parms$coefficients$fixed
        Phi <- parms$Phi
        # print(start)
        # print(coef)
        # print(start.Phi)
        # print(Phi)
        if(step.back) {
            if(control$trace)
                cat(" - new deviance = ",deviance)
            break
        }
        else {
            if(length(last.fit))
                last.eta <- last.fit$eta
            crit <- sum((eta - last.eta)^2) /sum(eta^2)
            
            if(control$trace)
                cat(" - criterion =",crit)
            
            if(crit <= control$eps){
                converged <- TRUE
                if(control$trace)
                    cat("\nconverged\n")
                break
            }
        }
    }
    if(!converged && !step.back){
        # if(control$trace) cat("\n")
        warning("Algorithm did not converge",call.=FALSE)
    }
    if(step.back){
        # if(control$trace) cat("\n")
        warning("Algorithm stopped without convergence",call.=FALSE)
    }
    eps <- 10*.Machine$double.eps
    if (any(pi < eps) || any(1-pi < eps)){
        # if(control$trace) cat("\n")
        warning("Fitted probabilities numerically 0 or 1 occurred",call.=FALSE)
    }
    if(deviance < 0){
        # if(control$trace) cat("\n")
        warning("Approximate deviance is negative.\nYou might be overfitting your data or the group size is too small.",call.=FALSE)
    }
    
    ntot <- length(y)
    pi0 <- mclogitP(offset,s)
    null.deviance <- sum(ifelse(y>0,
                    2*w*y*(log(y)-log(pi0)),
                    0))
    resid.df <- length(y) - length(unique(s))
    model.df <- ncol(X) + length(parms$lambda)
    resid.df <- resid.df - model.df
    
    return(
          list(
              coefficients = parms$coefficients$fixed,
              random.effects = parms$coefficients$random,
              VarCov = parms$Phi,
              lambda = parms$lambda,
              linear.predictors = eta,
              working.residuals = (y-pi)/pi,
              response.residuals = y-pi,
              df.residual = resid.df,
              model.df = model.df,
              deviance=deviance,
              deviance.residuals=dev.resids,
              null.deviance=null.deviance,
              method = method,
              estimator = estimator,
              iter = iter,
              y = y,
              s = s,
              offset = offset,
              converged = converged,
              control=control,
              info.coef = parms$info.fixed,
              info.fixed.random = parms$info.fixed.random,
              info.lambda = parms$info.lambda,
              info.psi = parms$info.psi
          ))
}

matrank <- function(x) {
    qr(x)$rank
}


PQLMQL_innerFit <- function(parms,aux,model.struct,method,estimator,control){

    m <- model.struct$m
    d <- model.struct$d
    nlevs <- model.struct$nlevs
    X <- model.struct$X
    Z <- model.struct$Z

    y <- aux$y
    W <- aux$W

    # Naive starting values
    Wy <- W%*%y
    WX <- W%*%X
    XWX <- crossprod(X,WX)
    XWy <- crossprod(X,Wy)
    yWy <- crossprod(y,Wy)
    
    alpha.start <- parms$coefficients$fixed
    Phi.start <- parms$Phi

    if(!length(alpha.start))
        alpha.start <- solve(XWX,XWy)

    y_Xalpha <- as.vector(y - X%*%alpha.start)

    if(!length(Phi.start)){
        Phi.start <- list()
        for(k in 1:nlevs){
            Z.k <- Z[[k]]
            ZZ.k <- crossprod(Z.k)
            Zy_Xa.k <- crossprod(Z.k,y_Xalpha)
            ZZ.k <- ZZ.k + Diagonal(ncol(ZZ.k))
            b.k <- solve(ZZ.k,Zy_Xa.k)
            m.k <- m[k]
            d.k <- d[k]
            dim(b.k) <- c(d.k,m.k)
            dimnames(b.k) <- NULL
            S.k <- tcrossprod(b.k)
            if(matrank(S.k) < d.k){
            #warning(sprintf("Singular initial covariance matrix at level %d in inner fitting routine",k))
                S.k <- diag(S.k)
                S.k <- diag(x=S.k,nrow=d)
            }
            Phi.start[[k]] <- S.k/(m.k-1)
        }
    }
    Psi.start <- lapply(Phi.start,safeInverse)
    Lambda.start <- lapply(Psi.start,chol)
    lambda.start <- unlist(lapply(Lambda.start,uvech))
    
    WZ <- bMatProd(W,Z)
    ZWZ <- bMatCrsProd(WZ,Z)
    ZWX <- bMatCrsProd(WZ,X)
    ZWy <- bMatCrsProd(WZ,y)

    aux <- list(yWy=yWy,
                XWy=XWy,
                ZWy=ZWy,
                XWX=XWX,
                ZWX=ZWX,
                ZWZ=ZWZ)

    if(control$trace.inner) cat("\n")

    devfunc <- function(lambda) 
       -2*as.vector(PQLMQL_pseudoLogLik(lambda,
                           model.struct=model.struct,
                           estimator=estimator,
                           aux=aux)$logLik)
    gradfunc <- function(lambda) 
       -2*as.vector(PQLMQL_pseudoLogLik(lambda,
                           model.struct=model.struct,
                           estimator=estimator,
                           aux=aux,
                           gradient=TRUE)$gradient) 

    if(control$inner.optimizer=="nlminb"){
    
        res.port <- nlminb(lambda.start,
                           objective = devfunc,
                           gradient = if(control$use.gradient == "analytic") gradfunc,
                           control = list(trace = as.integer(control$trace.inner),
                                          iter.max=control$maxit.inner)
                           )
        if(res.port$convergence != 0){
            cat("\n")
            warning(sprintf("Inner iterations did not coverge - nlminb message: %s",res.port$message),
                    call.=FALSE,immediate.=TRUE)
        }
        
        lambda <- res.port$par
    }
    else if(control$inner.optimizer=="nlm") {

        # 'nlminb' seems to be more stable - but this allows to check the analyticals.
        
        dev_f <- function(lambda){
            res <- PQLMQL_pseudoLogLik(lambda,
                                       model.struct=model.struct,
                                       estimator=estimator,
                                       aux=aux,
                                       gradient=TRUE)
            structure(-2*res$logLik,
                      gradient=if(control$use.gradient == "analytic") -2*res$gradient)
        }
        
        res.nlm <- nlm(f = dev_f,
                       p = lambda.start,
                       check.analyticals = TRUE,
                       print.level = if(control$trace.inner) 2 else 0,
                       iterlim = control$maxit.inner)
        
        if(res.nlm$code > 2){
            nlm.messages <- c("","",
                              paste("Last global step failed to locate a point lower than",
                                    "'estimate'.  Either 'estimate' is an approximate local",
                                    "minimum of the function or 'steptol' is too small.",sep="\n"),
                              "Iteration limit exceeded.",
                              paste("Maximum step size 'stepmax' exceeded five consecutive",
                                    "times.  Either the function is unbounded below, becomes",
                                    "asymptotic to a finite value from above in some direction",
                                    "or 'stepmax' is too small.",sep="\n"))
            retcode <- res.nlm$code
            cat("\n")
            warning(sprintf("Possible non-convergence of inner iterations - nlm code indicates:\n %s",
                            nlm.messages[retcode]),
                    call.=FALSE,immediate.=TRUE)
        }
        lambda <- res.nlm$estimate
    } else if(control$inner.optimizer %in% 
              c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN")){
        optim.method <- control$inner.optimizer
        optim.control <- list(
            trace = as.integer(control$trace.inner),
            maxit = control$maxit.inner,
            REPORT = switch(control$inner.optimizer,
                            SANN          = 100,
                            `Nelder-Mead` = 100,
                                            1),
            type  = if(optim.method == "CG") control$CG.type,
            alpha = if(optim.method == "Nelder-Mead") control$NM.alpha,
            beta  = if(optim.method == "Nelder-Mead") control$NM.beta,
            gamma = if(optim.method == "Nelder-Mead") control$NM.gamma,
            temp  = if(optim.method == "SANN") control$SANN.temp,
            tmax  = if(optim.method == "SANN") control$SANN.tmax
        )
        res.optim <- optim(par     = lambda.start,
                           fn      = devfunc,
                           gr      = if(control$use.gradient == "analytic") gradfunc,
                           method  = optim.method,
                           control = optim.control
                           )
        if(res.optim$convergence > 0){
            cat("\n")
            if(res.optim$convergence == 1) 
                warning("Inner iterations did not converge",
                        call.=FALSE,immediate.=TRUE)
            if(res.optim$convergence == 10) 
                warning("Degeneracy of the Nelder-Mead simplex",
                        call.=FALSE,immediate.=TRUE)
            if(length(res.optim$message))
                warning(sprintf("Message from 'optim':\n%s",
                                res.optim$message),
                        call.=FALSE,immediate.=TRUE)
        }
        lambda <- res.optim$par
    }
    else if(control$inner.optimizer == "ucminf" && 
            requireNamespace("ucminf", quietly = TRUE)){
        ucminf.control <- list(
            trace = as.integer(control$trace.inner)
            )
        for(nn in c("grtol","xtol","stepmax","maxeval","grad"))
            if(length(control[nn])) ucminf.control[[nn]] <- control[[nn]]
        res.ucminf <- ucminf::ucminf(par     = lambda.start,
                             fn      = devfunc,
                             gr      = if(control$use.gradient == "analytic") gradfunc,
                             control = ucminf.control
                           )
        if(res.ucminf$convergence > 2){
            cat("\n")
            if(length(res.ucminf$message))
                warning(sprintf("Message from 'ucminf':\n%s",
                                res.ucminf$message),
                        call.=FALSE,immediate.=TRUE)
        }
        else if(ucminf.control$trace > 0){
            cat("\n")
            if(length(res.ucminf$message))
                message(sprintf("Message from 'ucminf':\n%s",
                                res.ucminf$message))
        }
        lambda <- res.ucminf$par
    }
    else
        stop(sprintf("Unknown optimizer '%s'",control$inner.optimizer))

    info.varPar <- PQLMQL_pseudoLogLik(lambda,
                                       model.struct=model.struct,
                                       estimator=estimator,
                                       aux=aux,
                                       info.lambda=TRUE,
                                       info.psi=TRUE)$info

    Lambda <- lambda2Mat(lambda,m,d)
    Psi <- lapply(Lambda,crossprod)
    iSigma <- Psi2iSigma(Psi,m)
    Phi <- lapply(Psi,safeInverse)
    
    ZWZiSigma <- ZWZ + iSigma
    K <- solve(ZWZiSigma)


    log.det.iSigma <- Lambda2log.det.iSigma(Lambda,m)
    
    log.det.ZWZiSigma <- 2*sum(log(diag(chol_blockMatrix(ZWZiSigma,resplit=FALSE))))

    XiVX <- XWX - fuseMat(bMatCrsProd(ZWX,bMatProd(K,ZWX)))
    XiVy <- XWy - fuseMat(bMatCrsProd(ZWX,bMatProd(K,ZWy)))

    alpha <- solve(XiVX,XiVy)
    alpha <- drop(as.matrix(alpha))
    b <- bMatProd(K,ZWy-bMatProd(ZWX,alpha))
    b[] <- lapply(b[],as.matrix) 

    XZWiSZX <- structure(rbind(cbind(blockMatrix(XWX),bMatTrns(ZWX)),
                               cbind(ZWX,ZWZiSigma)),class="blockMatrix")

    list(
        lambda = lambda,
        coefficients = list(fixed = alpha,
                            random = b),
        Psi = Psi,
        Phi = Phi,
        info.fixed = as.matrix(XiVX),
        info.fixed.random = XZWiSZX,
        info.lambda = info.varPar$lambda,
        info.psi = info.varPar$psi,
        log.det.iSigma = log.det.iSigma,
        log.det.ZiVZ = log.det.ZWZiSigma,
        ZiVZ = ZWZiSigma
    )
 }

PQLMQL_eval_parms <- function(parms,model.struct,method,estimator){

    nlevs <- model.struct$nlevs
    d <- model.struct$d
    s <- model.struct$s
    y <- model.struct$y
    w <- model.struct$w

    X <- model.struct$X
    Z <- model.struct$Z
    offset <- model.struct$offset

    alpha <- parms$coefficients$fixed
    b <- parms$coefficients$random
    Psi <- parms$Psi
    ZiVZ <- parms$ZiVZ

    eta <- as.vector(X%*%alpha) + offset

    if(method=="PQL"){
        rand.ssq <- 0
        for(k in 1:nlevs){
            eta <- eta +  as.vector(Z[[k]]%*%b[[k]])
            B.k <- matrix(b[[k]],nrow=d[k])
            Psi.k <- Psi[[k]]
            rand.ssq <- rand.ssq + sum(B.k * (Psi.k%*%B.k))
        }
    } else {
        b_ <- blockMatrix(b,nrow=nlevs,ncol=1)
        rand.ssq <- as.vector(fuseMat(bMatCrsProd(b_,bMatProd(ZiVZ,b_))))
    }
    
    pi <- mclogitP(eta,s)
    dev.resids <- ifelse(y>0,
                         2*w*y*(log(y)-log(pi)),
                         0)
    
    deviance <-  -parms$log.det.iSigma + parms$log.det.ZiVZ + sum(dev.resids) + rand.ssq
    
    list(
        eta = eta,
        pi = pi,
        deviance = deviance
    )
}

log_Det <- function(x) determinant(x,logarithm=TRUE)$modulus


PQLMQL_pseudoLogLik <- function(lambda,
                                model.struct,
                                estimator,
                                aux,
                                gradient=FALSE,
                                info.lambda=FALSE,
                                info.psi=FALSE
                   ){

    nlevs <- model.struct$nlevs
    d <- model.struct$d
    m <- model.struct$m

    yWy <- aux$yWy
    XWy <- aux$XWy
    ZWy <- aux$ZWy
    XWX <- aux$XWX
    ZWX <- aux$ZWX
    ZWZ <- aux$ZWZ

    Lambda <- lambda2Mat(lambda,m,d)
    Psi <- lapply(Lambda,crossprod)
    iSigma <- Psi2iSigma(Psi,m)

    H <- ZWZ + iSigma
    if(getOption("mclogit.use_blkinv", TRUE)) {
        K <- blk_inv.squareBlockMatrix(H)
    }
    else {
        K <- solve(H)
    }

    XiVX <- XWX - fuseMat(bMatCrsProd(ZWX,bMatProd(K,ZWX)))
    XiVy <- XWy - fuseMat(bMatCrsProd(ZWX,bMatProd(K,ZWy)))
    XiVX <- symmpart(XiVX)
    
    alpha <- solve(XiVX,XiVy)
    b <- bMatProd(K,ZWy-bMatProd(ZWX,alpha))

    y.aXiVXa.y <- yWy - crossprod(XWy,alpha) - fuseMat(bMatCrsProd(ZWy,b))

    log.det.iSigma <- Lambda2log.det.iSigma(Lambda,m)
    
    log.det.H <- 2*sum(log(diag(chol_blockMatrix(H,resplit=FALSE))))
    logLik <- (log.det.iSigma - log.det.H - y.aXiVXa.y)/2
    if(estimator == "REML"){
        log.det.XiVX <- log_Det(XiVX)
        logLik <- logLik - log.det.XiVX/2
    }
    res <- list(
        logLik=as.vector(logLik),
        coefficients=as.vector(alpha),
        random.effects=b,
        Psi=Psi
        )

    if(gradient || info.lambda || info.psi){
        if(estimator=="REML"){
            iA <- solve(XiVX)
            XWZK <- bMatCrsProd(ZWX,K)
            iAXWZK <- bMatProd(blockMatrix(iA),XWZK)
            M <- bMatCrsProd(XWZK,iAXWZK)
        }
    }
    
    if(gradient){
        if(estimator=="REML"){
            K <- K + M
        }
        Phi <- lapply(Psi,safeInverse)
        S <- mapply(v_bCrossprod,b,d,SIMPLIFY=FALSE)
        K.kk <- diag(K)
        SumK.k <- mapply(sum_blockDiag,K.kk,d,SIMPLIFY=FALSE)
        Gr <- list()
        for(k in 1:nlevs)
            Gr[[k]] <- Lambda[[k]]%*%(m[k]*Phi[[k]] - SumK.k[[k]] - S[[k]])
        res$gradient <- unlist(lapply(Gr,uvech))
    }
    if(info.lambda || info.psi){
        res$info <- list()
        T <- iSigma - K
        if(estimator=="REML"){
            T <- T - M
        }
        if(info.lambda){
            G.lambda <- d.psi.d.lambda(Lambda)
            I.lambda <- blockMatrix(list(matrix(0,0,0)),nlevs,nlevs)
        }
        if(info.psi)
            I.psi <- blockMatrix(list(matrix(0,0,0)),nlevs,nlevs)
        for(k in 1:nlevs){
            T.k <- T[[k,k]]
            B.kk <- block_kronSum(T.k,m[k],m[k])
            if(info.lambda){
                G.k <- G.lambda[[k]]
                I.lambda[[k,k]] <- crossprod(G.k,B.kk%*%G.k)
            }
            if(info.psi){
                I.psi[[k,k]] <- B.kk/2
            }
            if(k < nlevs){
                for(k_ in seq(from=k+1,to=nlevs)){
                    T.kk_ <- T[[k,k_]]
                    B.kk_ <- block_kronSum(T.kk_,m[k],m[k_])
                    if(info.lambda){
                        G.k_ <- G.lambda[[k_]]
                        I.lambda[[k,k_]] <- crossprod(G.k,B.kk_%*%G.k_)
                        I.lambda[[k_,k]] <- t(I.lambda[[k,k_]])
                    }
                    if(info.psi){
                        I.psi[[k,k_]] <- B.kk_/2
                        I.psi[[k_,k]] <- t(I.psi[[k,k_]])
                    }
                }
            }
        }
        if(info.lambda)
            res$info$lambda <- as.matrix(fuseMat(I.lambda))
        if(info.psi)
            res$info$psi <- as.matrix(fuseMat(I.psi))
    }
    return(res)
}



vech <- function(x) x[lower.tri(x,diag=TRUE)]
setVech <- function(x,v) {
    ij <- lower.tri(x,diag=TRUE)
    x[ij] <- v
    x <- t(x)
    x[ij] <- v
    x
}

uvech <- function(x) x[upper.tri(x,diag=TRUE)]
set_uvech <- function(x,v,symm=FALSE) {
    ij <- upper.tri(x,diag=TRUE)
    x[ij] <- v
    if(symm){
        x <- t(x)
        x[ij] <- v
    }
    x
}
lambda2Mat <- function(lambda,m,d){
    nlevs <- length(m)
    dd2 <- d*(d+1)/2
    lambda <- split_(lambda,dd2)
    D <- lapply(d,diag)
    Map(set_uvech,D,lambda)
}

Psi2iSigma <- function(Psi,m){
    iSigma <- mapply(mk.iSigma.k,Psi,m,SIMPLIFY=FALSE)
    blockDiag(iSigma)
}

mk.iSigma.k <- function(Psi,m){
    Diagonal(m) %x% Psi
}

split_ <- function(x,d){
    m <- length(x)
    n <- length(d)
    i <- rep(1:n,d)
    split(x,i)
}

mmclogit.control <- function(
                             epsilon = 1e-08,
                             maxit = 25,
                             trace = TRUE,
                             trace.inner = FALSE,
                             avoid.increase = FALSE,
                             break.on.increase = FALSE,
                             break.on.infinite = FALSE,
                             break.on.negative = FALSE,
                             inner.optimizer = "nlminb",
                             maxit.inner = switch(inner.optimizer,
                                                  SANN          = 10000,
                                                  `Nelder-Mead` = 500,
                                                                  100),
                             CG.type = 1,
                             NM.alpha = 1,
                             NM.beta = 0.5,
                             NM.gamma = 2.0,
                             SANN.temp = 10,
                             SANN.tmax = 10,
                             grtol = 1e-6,
                             xtol = 1e-8,
                             maxeval = 100,
                             gradstep = c(1e-6, 1e-8),
                             use.gradient = c("analytic","numeric")) {
    if (!is.numeric(epsilon) || epsilon <= 0)
        stop("value of epsilon must be > 0")
    if (!is.numeric(maxit) || maxit <= 0)
        stop("maximum number of iterations must be > 0")
    m <- match.call()
    
    use.gradient <- match.arg(use.gradient)

    list(epsilon = epsilon, maxit = maxit,
         trace = trace, trace.inner = trace.inner,
         avoid.increase = avoid.increase,
         break.on.increase = break.on.increase,
         break.on.infinite = break.on.infinite,
         break.on.negative = break.on.negative,
         inner.optimizer = inner.optimizer,
         maxit.inner = maxit.inner,
         CG.type = CG.type,
         NM.alpha = NM.alpha,
         NM.beta = NM.beta,
         NM.gamma = NM.gamma,
         SANN.temp = SANN.temp,
         SANN.tmax = SANN.tmax,
         grtol = grtol,
         xtol = xtol,
         maxeval = maxeval,
         gradstep = gradstep,
         use.gradient = use.gradient
         )
}

split_bdiag1 <- function(x,n){
    m0 <- ncol(x)
    stopifnot(nrow(x)==m0)
    m <- m0%/%n
    i <- rep(1:m,each=n)
    j <- rep(1:m0)
    j <- split(j,i)
    y <- list()
    for(k in 1:m){
        j.k <- j[[k]]
        y[[k]] <- x[j.k,j.k]
    }
    y
}

split_bdiag <- function(x,d){
    m <- length(d)
    n <- ncol(x)
    s <- 1:m
    s <- rep(s,d)
    j <- 1:n
    j <- split(j,s)
    y <- list()
    for(k in 1:m){
        j.k <- j[[k]]
        y[[k]] <- x[j.k,j.k]
    }
    y
}


se_Phi <- function(Phi,info.lambda){
    d <- sapply(Phi,ncol)
    dd2 <- d*(d+1)/2
    info.lambda <- split_bdiag(info.lambda,dd2)
    Map(se_Phi_,Phi,info.lambda)
}

block_kronSum <- function(A,m1,m2){
    nr <- nrow(A)
    nc <- ncol(A)
    d1 <- nr%/%m1
    d2 <- nc%/%m2
    A <- as.array(A)
    dim(A) <- c(d1,m1,d2,m2)
    A <- aperm(A,c(2,4,1,3)) # dim = m1,m2,d1,d2
    dim(A) <- c(m1*m2,d1*d2)
    B <- crossprod(A) # dim = d1*d2,d1*d2
    dim(B) <- c(d1,d2,d1,d2)
    B <- aperm(B,c(1,3,2,4)) # dim = d1,d1,d2,d2
    dim(B) <- c(d1*d1,d2*d2)
    return(B)
}


d.psi.d.lambda <- function(Lambda) {
    lapply(Lambda,d.psi.d.lambda.1)
}

d.psi.d.lambda.1 <- function(Lambda){
    d <- ncol(Lambda)
    d_2 <- d*(d+1)/2
    G <- array(0,c(d,d,d,d))

    g <- rep(1:d,d*d*d)
    h <- rep(1:d,each=d,d*d)
    i <- rep(1:d,each=d*d,d)
    j <- rep(1:d,each=d*d*d)

    delta <- diag(d)
    
    G[cbind(g,h,i,j)] <- delta[cbind(g,j)]*Lambda[cbind(i,h)] + Lambda[cbind(i,g)]*delta[cbind(h,j)]
    
    dim(G) <- c(d*d,d*d)
    keep.lambda <- as.vector(upper.tri(Lambda,diag=TRUE))
    G[,keep.lambda]
}

solve_ <- function(x){
    res <- try(solve(x),silent=TRUE)
    if(inherits(res,"try-error")){
        warning("Singlular matrix encountered, trying a Moore-Penrose inverse")
        return(ginv(x))
    } else return(res)
        
}


se_Phi_ <- function(Phi,info.lambda){
    d <- ncol(Phi)
    Psi <- solve(Phi)
    Lambda <- chol(Psi)
    G <- d.psi.d.lambda.1(Lambda)
    vcov.lambda <- solve_(info.lambda)
    vcov.psi <- G%*%tcrossprod(vcov.lambda,G)
    PhiPhi <- Phi%x%Phi
    vcov.phi <- PhiPhi%*%vcov.psi%*%PhiPhi
    se.phi <- sqrt(diag(vcov.phi))
    matrix(se.phi,d,d,dimnames=dimnames(Phi))
}

Lambda2log.det.iSigma <- function(Lambda,m){
    res <- Map(Lambda2log.det.iSigma_1,Lambda,m)
    sum(unlist(res))
}

Lambda2log.det.iSigma_1 <- function(Lambda,m){
    dLambda <- diag(Lambda)
    if(any(dLambda < 0)){
        Psi <- crossprod(Lambda)
        svd.Psi <- svd(Psi)
        dLambda <- svd.Psi$d/2
    }
    2*m*sum(log(dLambda))
}


reff <- function(object){
    b <- object$random.effects
    Phi <- object$VarCov
    nlev <- length(b)
    B <- list()
    for(k in 1:nlev){
        d <- ncol(Phi[[k]])
        B_k <- matrix(b[[k]],nrow=d)
        B_k <- t(B_k)
        colnames(B_k) <- colnames(Phi[[k]])
        B[[k]] <- B_k
    }
    B
}
