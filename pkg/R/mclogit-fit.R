mclogit.fit <- function(
      y,
      s,
      w,
      X,
      dispersion=FALSE,
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
      coef <- start
    else coef <- NULL
    converged <- FALSE
    for(iter in 1:control$maxit){
        y.star <- eta - offset + (y-pi)/pi
        yP.star <- y.star - rowsum(pi*y.star,s)[s]
        XP <- X - as.matrix(rowsum(pi*X,s))[s,,drop=FALSE]
        ww <- w*pi
        good <- ww > 0 & is.finite(yP.star)
        wlsFit <- lm.wfit(x=XP[good,,drop=FALSE],y=yP.star[good],w=ww[good])
        last.coef <- coef
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
          if(!is.finite(deviance) || deviance > last.deviance && iter > 1){
            if(is.null(last.coef))
                stop("no valid set of coefficients has been found: please supply starting values", call. = FALSE)
             warning("step size truncated due to divergence", call. = FALSE)
             ii <- 1
             while (!is.finite(deviance) || deviance > last.deviance){
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
        crit <- abs(deviance-last.deviance)/abs(0.1+deviance)
        if(control$trace)
            cat("\nIteration",iter,"- deviance =",deviance,"- criterion =",crit)
        if(crit < control$eps){
          converged <- TRUE
          if(control$trace)
            cat("\nconverged\n")
          break
        }
    }
    if (!converged) warning("algorithm did not converge",call.=FALSE)
    if (boundary) warning("algorithm stopped at boundary value",call.=FALSE)
    eps <- 10*.Machine$double.eps
    if (any(pi < eps) || any(1-pi < eps))
       warning("fitted probabilities numerically 0 occurred",call.=FALSE)

    XP <- X - as.matrix(rowsum(pi*X,s))[s,,drop=FALSE]
    ww <- w*pi
    XWX <- crossprod(XP,ww*XP)
        
    ntot <- length(y)
    pi0 <- mclogitP(offset,s)
    null.deviance <- sum(ifelse(y>0,
                    2*w*y*(log(y)-log(pi0)),
                    0))
    resid.df <- length(y)-length(unique(s))
    model.df <- ncol(X)
    resid.df <- resid.df - model.df
    ll <- mclogit.logLik(y,pi,w)

    if(!isFALSE(dispersion)){
        if(isTRUE(dispersion))
            odisp.method <- "Afroz"
        else
            odisp.method <- match.arg(dispersion,
                                      c("Afroz",
                                        "Fletcher",
                                        "Pearson",
                                        "Deviance"))
        phi <- mclogit.dispersion(y,w,s,pi,coef,
                                      method=odisp.method)
    }
    else phi <- 1


    return(list(
        coefficients = drop(coef),
        phi = phi,
        linear.predictors = eta,
        working.residuals = (y-pi)/pi,
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
        information.matrix=XWX
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

