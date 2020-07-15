mclogit.dispersion <- function(y,w,s,pi,coef,method){
    N <- length(w)
    n <- length(unique(s))
    p <- length(coef)
    res.df <- N - n -p
    if(method=="Deviance"){
        Dresid <- 2*w*y*(log(y)-log(pi))
        Dresid[w==0 | y== 0] <- 0
        D <- sum(Dresid)
        phi <- D/res.df
    }
    else {
        X2 <- sum(w*(y - pi)^2/pi)
        phi.pearson <- X2/(N - n - p)
        if(method %in% c("Afroz","Fletcher"))
            s.bar <- sum((y - pi)/pi)/(N - n)
        phi <- switch(method,
                      Pearson = phi.pearson,
                      Afroz = phi.pearson/(1 + s.bar),
                      Fletcher = phi.pearson - (N - n)*s.bar/(N - n - p))
    }
    return(phi)
}

update_mclogit_dispersion <- function(object,dispersion){

    if(!missing(dispersion)){
        if(is.numeric(dispersion))
            phi <- dispersion
        else {
        if(isTRUE(dispersion))
            method <- "Afroz"
        else 
            method <- match.arg(dispersion,
                                      c("Afroz",
                                        "Fletcher",
                                        "Pearson",
                                        "Deviance"))
        phi <- dispersion(object,method=method)
        }
    }
    else phi <- 1

    object$phi <- phi
    return(object)
}

dispersion <- function(object,method,...)
    UseMethod("dispersion")

dispersion.mclogit <- function(object,method=NULL,...){
    if(is.null(method))
        return(object$phi)
    else {
        y <- object$y
        s <- object$s
        w <- object$weights
        pi <- object$fitted.values
        coef <- object$coefficients
        method <- match.arg(method,c("Afroz",
                                     "Fletcher",
                                     "Pearson",
                                     "Deviance"))
        phi <- mclogit.dispersion(y,w,s,pi,coef,
                                      method=method)
        return(phi)
    }
}
