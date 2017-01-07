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
            p.k <- ncol(Z[[k]])
            for(r in (k+1):nlevs){

                A.kr <- tcrossprod(ZWX.iXWX,ZWX[[r]])
                G.star.A.kr <- lapply(G.star.k,`%*%`,y=A.kr)

                hh.r <- length(G[[r]])
                S.kr <- matrix(0,hh.k,hh.r)

                h1 <- row(S.kr)
                h2 <- col(S.kr)
                S.kr[] <- mapply(tr.crossprod,G.star.A.kr[h1],G.star.A.kr[h2])
                S[[k,r]] <- S.kr
                S[[r,k]] <- t(S.kr)                
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

    ## Check for negative variances
    diagPhi <- unlist(lapply(Phi,diag))
    if(any(diagPhi < 0)){
        warning("Moment equations give negative variances.
  Your model appears to be misspecified.
  I will use a dummy covariance matrix.",call.=FALSE)
        for(k in 1:nlevs){
            Phi[[k]] <- diag(ncol(Phi[[k]]))
        }
    }
    
    iPhi <- lapply(Phi,solve)

    Sigma <- Map(`%x%`,I.mk,Phi)
    Sigma <- do.call(bdiag,unname(Sigma))

    iSigma <- Map(`%x%`,I.mk,iPhi)
    iSigma <- do.call(bdiag,unname(iSigma))

    b.split <- rep(1:nlevs,sapply(Z,ncol))
    
    converged <- FALSE
   
    ## IWLS for coefficients and variance parameters    

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

    ZWZ.iSigma <- fuseMat(ZWZ) + iSigma
    K <- solve(ZWZ.iSigma)
    b <- K%*%(ZWy-ZWX%*%coef)    
    
    for(iter in 1:control$maxit){
        
        last.deviance <- deviance
        last.coef <- coef
        last.theta <- theta

        ## General preparations

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

        ZWZ.iSigma <- fuseMat(ZWZ)+iSigma
        ZWX <- fuseMat(ZWX)
        ZWy <- fuseMat(ZWy)
        K <- solve(ZWZ.iSigma)
        
        ## Update coef 

        XiVX <- XWX - crossprod(ZWX,K%*%ZWX)
        XiVy <- XWy - crossprod(ZWX,K%*%ZWy)
        coef <- solve(XiVX,XiVy)

        ## Update b

        b <- K%*%(ZWy-ZWX%*%coef)
        
        ## Update theta

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

                p.k <- ncol(Z[[k]])
                for(r in (k+1):nlevs){
                    
                    ZWZ.kr <- ZWZ[[k,r]]
                    ZWZ.r <- ZWZ.[[r]]

                    A.kr <- ZWZ.kr - crossprod(ZWZ.k,K%*%ZWZ.r)
                    G.star.A.kr <- lapply(G.star.k,`%*%`,y=A.kr)

                    hh.r <- length(G[[r]])
                    S.kr <- matrix(0,hh.k,hh.r)

                    h1 <- row(S.kr)
                    h2 <- col(S.kr)
                    S.kr[] <- mapply(tr.crossprod,G.star.A.kr[h1],G.star.A.kr[h2])
                    S[[k,r]] <- S.kr
                    S[[r,k]] <- t(S.kr)
                }
            }
        }
        
        u <- unlist(u)
        Info.theta <- fuseMat(S)

        theta <- solve(Info.theta,u)
        theta <- split(theta,vpar.selector)

        
        ## Compute deviance and determine step size

        stepsize.found <- FALSE

        ZWZ <- fuseMat(ZWZ)
        for(iiter in 1:control$maxit) {

            K <- solve(ZWZ.iSigma)
            b <- K%*%(ZWy-ZWX%*%coef)
        
            b. <- split(b,b.split)
            Zb <- Map(`%*%`,Z,b.)
            Zb <- Reduce(`+`,Zb)
            
            eta <- as.vector(X%*%coef) + as.vector(Zb) + offset
            pi <- mclogitP(eta,s)

            dev.resids <- ifelse(y>0,2*w*y*(log(y)-log(pi)),0)

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
            
            ZWZ.iSigma <- ZWZ + iSigma
            b.iSigma.b <- as.numeric(crossprod(b,iSigma%*%b))

            logDet.ZWZ.iSigma <- log.Det(ZWZ.iSigma)
            
            deviance <- sum(dev.resids) + logDetSigma + logDet.ZWZ.iSigma + b.iSigma.b
            crit <- abs(deviance-last.deviance)/abs(0.1+deviance)

            if(deviance < last.deviance || crit < control$eps){
                stepsize.found <- TRUE
                break
            }

            coef <- (coef + last.coef)/2
            theta <- Map(`/`,Map(`+`,theta,last.theta),2)
            #if(control$trace) cat("\n\tDeviance =",deviance)
        }
        if(!stepsize.found){
            cat("\n")
            warning("Cannot find an appropriate step size, giving up",
                    call.=FALSE)
            break
        }
        
        
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
    if (!converged) warning("Algorithm did not converge.",
                    call.=FALSE)

    W <- Matrix(0,nrow=nobs,ncol=nsets)
    W[cbind(i,s)] <- sqrt.w*pi
    W <- Diagonal(x=w*pi)-tcrossprod(W)

    WX <- W%*%X
    XWX <- crossprod(X,WX)
    iXWX <- solve(XWX)

    WZ <- lapply(Z,`%*%`,x=W)
    ZWX <- lapply(Z,crossprod,y=WX)
    ZWy <- lapply(Z,crossprod,y=Wy)
    dim(ZWX) <- c(nlevs,1)
    dim(ZWy) <- c(nlevs,1)
    
    ZWX <- fuseMat(ZWX)
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

