mclogit.overdispersion <- function(y,w,s,pi,coef,method){
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

update_mclogit_overdispersion <- function(object,overdispersion){

    if(!isFALSE(overdispersion)){
        if(is.numeric(overdispersion))
            phi <- overdispersion
        else {
        if(isTRUE(overdispersion))
            odisp.method <- "Afroz"
        else 
            odisp.method <- match.arg(overdispersion,
                                      c("Afroz",
                                        "Fletcher",
                                        "Pearson",
                                        "Deviance"))
        phi <- overdispersion(object,method=odisp.method)
        }
    }
    else phi <- 1

    object$phi <- phi
    return(object)
}

overdispersion <- function(object,method=NULL,...){
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
        phi <- mclogit.overdispersion(y,w,s,pi,coef,
                                      method=method)
        return(phi)
    }
}
