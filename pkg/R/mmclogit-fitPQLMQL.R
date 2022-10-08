mmclogit.fitPQLMQL <- function(
                               y,
                               s,
                               w,
                               X,
                               Z,
                               d,
                               start = NULL,
                               offset = NULL,
                               method = c("PQL","MQL"),
                               estimator = c("ML","REML"),
                               control=mmclogit.control()
                               ){

    method <- match.arg(method)
    
    nvar <- ncol(X)
    nobs <- length(y)
    nsets <- length(unique(s))
    nlevs <- length(Z)

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
    
    
    for(iter in 1:control$maxit){

        W <- Matrix(0,nrow=nobs,ncol=nsets)
        W[cbind(i,s)] <- sqrt.w*pi
        W <- Diagonal(x=w*pi)-tcrossprod(W)

        y.star <- eta - offset + (y-pi)/pi

        prev.last.parms <- last.parms
        last.parms <- parms
        
        parms <- try(PQLMQL_innerFit(y.star,X,Z,W,d,offset,method,estimator,control),
                   silent=TRUE)

        step.back <- FALSE
        if(inherits(parms,"try-error")){
            if(length(prev.last.deviance) && 
               prev.deviance > prev.last.deviance && 
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
        fit <- PQLMQL_eval_parms(y,s,w,parms,X,Z,d,offset,method,estimator)
        
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


PQLMQL_innerFit <- function(y,X,Z,W,d,offset,method,estimator,control){

    nlevs <- length(Z)

    m <- sapply(Z,ncol)/d

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
        ZZ.k <- crossprod(Z.k)
        ZyXa.k <- crossprod(Z.k,y.Xalpha)
        ZZ.k <- ZZ.k + Diagonal(ncol(ZZ.k))
        b.k <- solve(ZZ.k,ZyXa.k)
        m.k <- m[k]
        d.k <- d[k]
        dim(b.k) <- c(d.k,m.k)
        S.k <- tcrossprod(b.k)
        if(matrank(S.k) < d.k){
            #warning(sprintf("Singular initial covariance matrix at level %d in inner fitting routine",k))
            S.k <- diag(S.k)
            S.k <- diag(x=S.k,nrow=d)
        }
        Phi.start[[k]] <- S.k/(m.k-1)
    }

    Psi.start <- lapply(Phi.start,safeInverse)
    Lambda.start <- lapply(Psi.start,chol)
    lambda.start <- unlist(lapply(Lambda.start,uvech))
    
    WZ <- bMatProd(W,Z)
    ZWZ <- bMatCrsProd(WZ,Z)
    ZWX <- bMatCrsProd(WZ,X)
    ZWy <- bMatCrsProd(WZ,y)

    nqfunc1 <- function(lambda) -qfunc(lambda,y,d,yWy,XWy,ZWy,XWX,ZWX,ZWZ,estimator)
    ngrfunc1 <- function(lambda) -grfunc(lambda,y,d,yWy,XWy,ZWy,XWX,ZWX,ZWZ,estimator)

    info_func1 <- function(lambda) info_func(lambda,d,XWX,ZWX,ZWZ,estimator)

    if(control$trace.inner) cat("\n")
    res.port <- nlminb(lambda.start,
                       objective = nqfunc1,
                       gradient = ngrfunc1,
                       control = list(trace = as.integer(control$trace.inner))
                       )

    lambda <- res.port$par
    info.lambda <- info_func1(lambda)
    info.psi <- info_func_psi(lambda,d,XWX,ZWX,ZWZ,estimator)

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
        info.lambda = info.lambda,
        info.psi = info.psi,
        log.det.iSigma = log.det.iSigma,
        log.det.ZiVZ = log.det.ZWZiSigma,
        ZiVZ = ZWZiSigma
    )
 }

PQLMQL_eval_parms <- function(y,s,w,parms,X,Z,d,offset,method,estimator){

    nlevs <- length(Z)
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
        b_ <- blockMatrix(b,nrow=nlevs)
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


qfunc <- function(lambda,y,d,yWy,XWy,ZWy,XWX,ZWX,ZWZ,estimator){

    estimator <- match.arg(estimator,c("ML","REML"))
    nlevs <- ncol(ZWZ)
    m <- bM_ncol(ZWZ)%/%d
    Lambda <- lambda2Mat(lambda,m,d)
    Psi <- lapply(Lambda,crossprod)
    iSigma <- Psi2iSigma(Psi,m)

    H <- ZWZ + iSigma
    K <- solve(H)

    XiVX <- XWX - fuseMat(bMatCrsProd(ZWX,bMatProd(K,ZWX)))
    XiVy <- XWy - fuseMat(bMatCrsProd(ZWX,bMatProd(K,ZWy)))
    XiVX <- symmpart(XiVX)
    
    alpha <- solve(XiVX,XiVy)
    b <- bMatProd(K,ZWy-bMatProd(ZWX,alpha))

    y.aXiVXa.y <- yWy - crossprod(XWy,alpha) - fuseMat(bMatCrsProd(ZWy,b))
    #S <- mapply(v_bCrossprod,b,d)
    #bPsib <- mapply(`%*%`,S,Psi)
    #bPsib <- Reduce(`+`,bPsib)

    
    #log.det.iSigma <- 2*sum(log(diag(chol_blockMatrix(iSigma,resplit=FALSE))))
    log.det.iSigma <- Lambda2log.det.iSigma(Lambda,m)
    
    log.det.H <- 2*sum(log(diag(chol_blockMatrix(H,resplit=FALSE))))
    res <- (log.det.iSigma - log.det.H - y.aXiVXa.y)/2
    if(estimator == "REML"){
        log.det.XiVX <- log.Det(XiVX)
        res <- res - log.det.XiVX/2
    }
    as.vector(res)
}

grfunc <- function(lambda,y,d,yWy,XWy,ZWy,XWX,ZWX,ZWZ,estimator){

    estimator <- match.arg(estimator,c("ML","REML"))
    nlevs <- ncol(ZWZ)
    m <- bM_ncol(ZWZ)%/%d
    Lambda <- lambda2Mat(lambda,m,d)
    Psi <- lapply(Lambda,crossprod)
    iSigma <- Psi2iSigma(Psi,m)

    H <- ZWZ + iSigma
    K <- solve(H)

    XiVX <- XWX - fuseMat(bMatCrsProd(ZWX,bMatProd(K,ZWX)))
    XiVy <- XWy - fuseMat(bMatCrsProd(ZWX,bMatProd(K,ZWy)))
    XiVX <- symmpart(XiVX)

    if(estimator=="REML"){
        iA <- solve(XiVX)
        alpha <- iA%*%XiVy
    }
    else
        alpha <- solve(XiVX,XiVy)
    b <- bMatProd(K,ZWy-bMatProd(ZWX,alpha))

    if(estimator=="REML"){
        XWZK <- bMatCrsProd(ZWX,K)
        iAXWZK <- bMatProd(blockMatrix(iA),XWZK)
        M <- bMatCrsProd(XWZK,iAXWZK)
        K <- K + M
    }
    Phi <- lapply(Psi,safeInverse)
    S <- mapply(v_bCrossprod,b,d,SIMPLIFY=FALSE)
    K.kk <- diag(K)
    SumK.k <- mapply(sum_blockDiag,K.kk,d,SIMPLIFY=FALSE)
    Gr <- list()
    for(k in 1:nlevs)
        Gr[[k]] <- Lambda[[k]]%*%(m[k]*Phi[[k]] - SumK.k[[k]] - S[[k]])
    gr <- lapply(Gr,uvech)
    unlist(gr)
}

info_func_psi <- function(lambda,d,XWX,ZWX,ZWZ,estimator){

    estimator <- match.arg(estimator,c("ML","REML"))
    nlevs <- ncol(ZWZ)
    m <- bM_ncol(ZWZ)%/%d

    Lambda <- lambda2Mat(lambda,m,d)
    Psi <- lapply(Lambda,crossprod)
    iSigma <- Psi2iSigma(Psi,m)

    H <- ZWZ + iSigma
    K <- solve(H)
    T <- iSigma - K

    if(estimator=="REML"){
        XiVX <- XWX - fuseMat(bMatCrsProd(ZWX,bMatProd(K,ZWX)))
        XiVX <- symmpart(XiVX)
        iA <- solve(XiVX)
        XWZK <- bMatCrsProd(ZWX,K)
        iAXWZK <- bMatProd(blockMatrix(iA),XWZK)
        M <- bMatCrsProd(XWZK,iAXWZK)
        T <- T - M
    }
    
    Imat <- blockMatrix(list(matrix(0,0,0)),nlevs,nlevs)
    for(k in 1:nlevs) {
        T.k <- T[[k,k]]
        B.kk <- block_kronSum(T.k,m[k],m[k])
        Imat[[k,k]] <- B.kk/2
        if(k < nlevs)
            for(k_ in seq(from=k+1,to=nlevs)){
                T.kk_ <- T[[k,k_]]
                B.kk_ <- block_kronSum(T.kk_,m[k],m[k_])
                Imat[[k,k_]] <- B.kk_/2
                Imat[[k_,k]] <- t(Imat[[k,k_]])
            }
    }
    as.matrix(fuseMat(Imat))
} 


info_func <- function(lambda,d,XWX,ZWX,ZWZ,estimator){

    estimator <- match.arg(estimator,c("ML","REML"))
    nlevs <- ncol(ZWZ)
    m <- bM_ncol(ZWZ)%/%d

    Lambda <- lambda2Mat(lambda,m,d)
    Psi <- lapply(Lambda,crossprod)
    iSigma <- Psi2iSigma(Psi,m)

    G.lambda <- d.psi.d.lambda(Lambda)
    
    H <- ZWZ + iSigma
    K <- solve(H)
    T <- iSigma - K

    if(estimator=="REML"){
        XiVX <- XWX - fuseMat(bMatCrsProd(ZWX,bMatProd(K,ZWX)))
        XiVX <- symmpart(XiVX)
        iA <- solve(XiVX)
        XWZK <- bMatCrsProd(ZWX,K)
        iAXWZK <- bMatProd(blockMatrix(iA),XWZK)
        M <- bMatCrsProd(XWZK,iAXWZK)
        T <- T - M
    }
    
    Imat <- blockMatrix(list(matrix(0,0,0)),nlevs,nlevs)
    for(k in 1:nlevs) {
        T.k <- T[[k,k]]
        B.kk <- block_kronSum(T.k,m[k],m[k])
        G.k <- G.lambda[[k]]
        Imat[[k,k]] <- crossprod(G.k,B.kk%*%G.k)/2
        if(k < nlevs)
            for(k_ in seq(from=k+1,to=nlevs)){
                T.kk_ <- T[[k,k_]]
                B.kk_ <- block_kronSum(T.kk_,m[k],m[k_])
                G.k_ <- G.lambda[[k_]]
                Imat[[k,k_]] <- crossprod(G.k,B.kk_%*%G.k_)/2
                Imat[[k_,k]] <- t(Imat[[k,k_]])
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
                             break.on.negative = FALSE
                            ) {
    if (!is.numeric(epsilon) || epsilon <= 0)
        stop("value of epsilon must be > 0")
    if (!is.numeric(maxit) || maxit <= 0)
        stop("maximum number of iterations must be > 0")
    list(epsilon = epsilon, maxit = maxit,
         trace = trace, trace.inner = trace.inner,
         avoid.increase = avoid.increase,
         break.on.increase = break.on.increase,
         break.on.infinite = break.on.infinite,
         break.on.negative = break.on.negative
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

se_Phi_ <- function(Phi,info.lambda){
    d <- ncol(Phi)
    Psi <- solve(Phi)
    Lambda <- chol(Psi)
    G <- d.psi.d.lambda.1(Lambda)
    vcov.lambda <- solve(info.lambda)
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
    m*sum(log(dLambda))
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
