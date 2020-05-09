mmclogit.fitPQL <- function(
                            y,
                            s,
                            w,
                            X,
                            Z,
                            groups,
                            start,
                            offset=NULL,
                            control=mmclogit.control()
      ){

    nvar <- ncol(X)
    nobs <- length(y)
    nsets <- length(unique(s))
    nlevs <- length(groups)

    Znm <- colnames(Z)
    
    sqrt.w <- sqrt(w)

    i <- 1:nobs
    
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

    Z0 <- Z
    d <- ncol(Z0)
    # Expand random effects design matrix from intecepts and slopes for random effects
    # for each level
    Z <- blockMatrix(lapply(groups,mkZ,rX=Z0))
    # Outer iterations: update non-linear part of the model
    converged <- FALSE
    for(iter in 1:control$maxit){

        W <- Matrix(0,nrow=nobs,ncol=nsets)
        W[cbind(i,s)] <- sqrt.w*pi
        W <- Diagonal(x=w*pi)-tcrossprod(W)

        y.star <- eta - offset + (y-pi)/pi
        ww <- w*pi
        good <- ww > 0

        fit <- PQLinnerFit(y.star,X,Z,W,d,groups,offset,control)

        coef <- fit$coefficients
        eta <- as.vector(X%*%coef$fixed) + offset
        for(k in 1:nlevs){
            eta <- eta +  as.vector(Z[[k]]%*%coef$random[[k]])
        }
        
        pi <-   mclogitP(eta,s)
        last.deviance <- deviance
        dev.resids <- ifelse(y>0,
                2*w*y*(log(y)-log(pi)),
                0)
        log.det.iSigma <- fit$log.det.iSigma
        log.det.ZWZiSigma <- fit$log.det.ZWZiSigma
        deviance <- sum(dev.resids) - log.det.iSigma + log.det.ZWZiSigma
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
    #eps <- 10*.Machine$double.eps
    #if (any(pi < eps) || any(1-pi < eps))
    #    warning("fitted probabilities numerically 0 occurred")

    info.coef <- fit$info.fixed
    info.theta <- fit$info.theta
    G <- fit$G
    
    Phi <- fit$Phi
    for(k in seq_along(Phi))
        dimnames(Phi[[k]]) <- list(Znm,Znm)

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
        coefficients = coef$fixed,
        VarCov = Phi,
        linear.predictors = eta,
        working.residuals = (y-pi)/pi,
        working.weights = w,
        response.residuals = y-pi,
        df.residual = resid.df,
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
        info.coef = info.coef,
        info.theta = info.theta,
        G = G
        ))
}


PQLinnerFit <- function(y,X,Z,W,d,groups,offset,control){

    nlevs <- length(groups)
    m <- lapply(groups,lunq)

    # Design matrix for variance parameters
    G <- mapply(Gfunc,d,m,SIMPLIFY=TRUE)
    
    # Naive starting values
    
    Wy <- W%*%y
    WX <- W%*%X
    XWX <- crossprod(X,WX)
    XWy <- crossprod(X,Wy)
    yWy <- crossprod(y,Wy)
    
    alpha.start <- solve(XWX,XWy)

    y.Xalpha <- as.vector(y - X%*%alpha.start)

    Phi.start <- list()
    for(k in 1:nlevs){
        Z.k <- Z[[k]]
        b.k <- solve(crossprod(Z.k),crossprod(Z.k,y.Xalpha))
        m.k <- m[[k]]
        dim(b.k) <- c(d,m.k)
        S.k <- tcrossprod(b.k)
        Phi.start[[k]] <- S.k/(m.k-1)
    }

    phi.start <- unlist(lapply(Phi.start,vech))
    Psi.start <- lapply(Phi.start,solve)
    theta.start <- unlist(lapply(Psi.start,vech))

    WZ <- bMatProd(W,Z)
    ZWZ <- bMatCrsProd(WZ,Z)
    ZWX <- bMatCrsProd(WZ,X)
    ZWy <- bMatCrsProd(WZ,y)

    qval <- qfunc(theta.start,y,d,yWy,XWy,ZWy,XWX,ZWX,ZWZ)
    qfunc1 <- function(theta) qfunc(theta,y,d,yWy,XWy,ZWy,XWX,ZWX,ZWZ)
    nqfunc1 <- function(theta) -qfunc(theta,y,d,yWy,XWy,ZWy,XWX,ZWX,ZWZ)
    #qfuncv <- Vectorize(qfunc1)

    gval <- grfunc(theta.start,y,d,yWy,XWy,ZWy,XWX,ZWX,ZWZ)
    grfunc1 <- function(theta) grfunc(theta,y,d,yWy,XWy,ZWy,XWX,ZWX,ZWZ)
    ngrfunc1 <- function(theta) -grfunc(theta,y,d,yWy,XWy,ZWy,XWX,ZWX,ZWZ)

    #grfuncv <- Vectorize(grfunc1)
    # layout(1:3)
    # curve(qfuncv(x),from=2.75,to=3)
    # points(theta.start,qval)
    # curve(grfuncv(x),from=2.75,to=3)
    # points(theta.start,gval)
    # abline(h=0)

    info_func1 <- function(theta) info_func(theta,d,ZWZ,G)
    #hess_func1 <- function(theta) -info_func(theta,d,ZWZ,G)
    #info_funcv <- function(theta) {s1 <- lapply(theta,info_func1);sapply(s1,as.vector)}
    #curve(info_funcv(x),from=2.75,to=3)

    if(control$trace.inner) cat("\n")
    res.port <- nlminb(theta.start,
                       objective = nqfunc1,
                       gradient = ngrfunc1,
                       control = list(trace = as.integer(control$trace.inner))
                       )

    theta <- res.port$par
    info.theta <- info_func1(theta)

    Psi <- theta2Psi(theta,m,d)
    iSigma <- Psi2iSigma(Psi,m)

    Phi <- lapply(Psi,solve)
    
    H <- ZWZ + iSigma
    K <- solve(H)

    log.det.iSigma <- 2*sum(log(diag(chol_blockMatrix(iSigma,resplit=FALSE))))
    log.det.ZWZiSigma <- 2*sum(log(diag(chol_blockMatrix(H,resplit=FALSE))))

    XiVX <- XWX - fuseMat(bMatCrsProd(ZWX,bMatProd(K,ZWX)))
    XiVy <- XWy - fuseMat(bMatCrsProd(ZWX,bMatProd(K,ZWy)))

    alpha <- solve(XiVX,XiVy)
    alpha <- drop(as.matrix(alpha))
    b <- bMatProd(K,ZWy-bMatProd(ZWX,alpha))
    b[] <- lapply(b[],as.matrix) 
    
    list(
        theta = theta,
        coefficients = list(fixed = alpha,
                            random = b),
        Psi = Psi,
        Phi = Phi,
        info.fixed = as.matrix(XiVX),
        info.theta = info.theta,
        log.det.iSigma   = log.det.iSigma,
        log.det.ZWZiSigma = log.det.ZWZiSigma,
        G = G
    )
 }

qfunc <- function(theta,y,d,yWy,XWy,ZWy,XWX,ZWX,ZWZ){
    nlevs <- ncol(ZWZ)
    m <- bM_ncol(ZWZ)%/%d
    Psi <- theta2Psi(theta,m,d)
    iSigma <- Psi2iSigma(Psi,m)

    H <- ZWZ + iSigma
    K <- solve(H)

    XiVX <- XWX - fuseMat(bMatCrsProd(ZWX,bMatProd(K,ZWX)))
    XiVy <- XWy - fuseMat(bMatCrsProd(ZWX,bMatProd(K,ZWy)))

    alpha <- solve(XiVX,XiVy)
    b <- bMatProd(K,ZWy-bMatProd(ZWX,alpha))

    y.aXiVXa.y <- yWy - crossprod(XWy,alpha) - fuseMat(bMatCrsProd(ZWy,b))
    S <- mapply(v_bCrossprod,b,d)
    bPsib <- mapply(`%*%`,S,Psi)
    bPsib <- Reduce(`+`,bPsib)
    
    log.det.iSigma <- 2*sum(log(diag(chol_blockMatrix(iSigma,resplit=FALSE))))
    log.det.H <- 2*sum(log(diag(chol_blockMatrix(H,resplit=FALSE))))
    res <- (log.det.iSigma - log.det.H - y.aXiVXa.y)/2
    as.vector(res)
}

grfunc <- function(theta,y,d,yWy,XWy,ZWy,XWX,ZWX,ZWZ){
    nlevs <- ncol(ZWZ)
    m <- bM_ncol(ZWZ)%/%d
    Psi <- theta2Psi(theta,m,d)
    iSigma <- Psi2iSigma(Psi,m)

    H <- ZWZ + iSigma
    K <- solve(H)

    XiVX <- XWX - fuseMat(bMatCrsProd(ZWX,bMatProd(K,ZWX)))
    XiVy <- XWy - fuseMat(bMatCrsProd(ZWX,bMatProd(K,ZWy)))

    alpha <- solve(XiVX,XiVy)
    b <- bMatProd(K,ZWy-bMatProd(ZWX,alpha))

    Phi <- lapply(Psi,solve)
    S <- mapply(v_bCrossprod,b,d)
    K.kk <- diag(K)
    SumK.k <- mapply(sum_blockDiag,K.kk,d)
    Gr <- list()
    for(k in 1:nlevs)
        Gr[[k]] <- (m[k]*Phi[[k]] - SumK.k[[k]] - S[[k]])/2
    gr <- lapply(Gr,vech)
    unlist(gr)
}

info_func <- function(theta,d,ZWZ,G){
    nlevs <- ncol(ZWZ)
    m <- bM_ncol(ZWZ)%/%d
    Psi <- theta2Psi(theta,m,d)
    iSigma <- Psi2iSigma(Psi,m)

    H <- ZWZ + iSigma
    K <- solve(H)
    T <- iSigma - K
    Imat <- blockMatrix(list(matrix(0,0,0)),nlevs,nlevs)
    for(k in 1:nlevs) {
        T.k <- T[[k,k]]
        G.k <- G[[k]]
        TG.k <- T.k%*%G.k
        Imat[[k,k]] <- crossprod(TG.k)
        if(k < nlevs)
            for(k_ in seq(from=k+1,to=nlevs)){
                G.k_ <- G[[k_]]
                TG.kk_ <- T[[k,k_]]%*%G.k_
                TG.k_k <- T[[k_,k]]%*%G.k
                Imat[[k_,k]] <- crossprod(TG.kk_,TG.k_k)
            }
    }
    as.matrix(fuseMat(Imat))
} 

vech <- function(x) x[lower.tri(x,diag=TRUE)]
setVech <- function(x,v) {
    ij <- lower.tri(x,diag=TRUE)
    x[ij] <- v
    x <- t(x)
    x[ij] <- v
    x
}

theta2Psi <- function(theta,m,d){
    nlevs <- length(m)
    dd2 <- d*(d+1)%/%2
    theta <- split_(theta,dd2)
    lapply(theta,setVech,x=diag(d))
}

Psi2iSigma <- function(Psi,m){
    iSigma <- mapply(mk.iSigma.k,Psi,m,SIMPLIFY=FALSE)
    blockDiag(iSigma)
}

theta2iSigma <- function(theta,m,d){
    nlevs <- length(m)
    dd2 <- d*(d+1)%/%2
    theta <- split_(theta,dd2)
    Psi <- lapply(theta,setVech,x=diag(d))
    iSigma <- mapply(mk.iSigma.k,Psi,m,SIMPLIFY=FALSE)
    blockDiag(iSigma)
}

mk.iSigma.k <- function(Psi,m){
    d <- ncol(Psi)
    iSigma.k <- Matrix(0,d*m,d*m)
    set_blockDiag(iSigma.k,Psi)
}

split_ <- function(x,n){
    m <- length(x)
    m <- m%/%n
    i <- rep(1:m,each=n)
    split(x,i)
}

vec.vech <- function(d){
    dd <- d*d
    dd2 <- (d*(d+1))%/%2
    I <- matrix(1:dd,d,d)
    J <- I*0
    J[lower.tri(J,diag=TRUE)] <- 1:dd2
    J[upper.tri(J)] <- t(J)[upper.tri(J)]
    ij <- cbind(as.vector(I),as.vector(J))
    res <- Matrix(0,dd,dd2)
    res[ij] <- 1
    res
}

Gfunc <- function(d,m){
    G_ <- vec.vech(d)
    do.call(rbind,rep(list(G_),m))
}

mmclogit.control <- function(
                             epsilon = 1e-08,
                             maxit = 25,
                             trace = TRUE,
                             trace.inner = FALSE
                            ) {
    if (!is.numeric(epsilon) || epsilon <= 0)
        stop("value of epsilon must be > 0")
    if (!is.numeric(maxit) || maxit <= 0)
        stop("maximum number of iterations must be > 0")
    list(epsilon = epsilon, maxit = maxit,
         trace = trace, trace.inner = trace.inner)
}

se_Phi_ <- function(Phi,info.theta){
    d <- ncol(Phi)
    Phi.Phi <- Phi %x% Phi
    G.theta <- vec.vech(d)
    vcov.theta <- solve(info.theta)
    vcov.psi <- G.theta %*% tcrossprod(vcov.theta,G.theta)
    vcov.phi <- Phi.Phi %*% vcov.psi %*% Phi.Phi
    se.phi <- sqrt(diag(vcov.phi))
    matrix(se.phi,d,d,dimnames=dimnames(Phi))
}

split_bdiag <- function(x,n){
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

se_Phi <- function(Phi,info.theta){
    d <- ncol(Phi[[1]])
    info.theta <- split_bdiag(info.theta,d)
    Map(se_Phi_,Phi,info.theta)
}
