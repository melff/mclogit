mclogit.dispersion <- function(y,w,s,pi,coef,method,groups=NULL){
    if(!length(groups))
        mclogit.dispersion1(y,w,s,pi,coef,method)
    else {
        agg <- aggregate_responses(y,w,s,pi,groups)
        y <- agg$y
        w <- agg$w
        s <- agg$s
        pi <- agg$pi
        mclogit.dispersion1(y,w,s,pi,coef,method)
    } 
}

mclogit.dispersion1 <- function(y,w,s,pi,coef,method){
    w1 <- split(w,s)
    w1 <- unlist(lapply(w1,"[",1))
    N <- sum(w1)
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

dispersion.mclogit <- function(object,method=NULL,groups=NULL,...){
    if(is.null(method))
        return(object$phi)
    else {
        y <- object$y
        s <- object$s
        w <- object$weights
        pi <- object$fitted.values
        coef <- object$coefficients
        if(!length(groups))
            groups <- object$groups
        method <- match.arg(method,c("Afroz",
                                     "Fletcher",
                                     "Pearson",
                                     "Deviance"))
        phi <- mclogit.dispersion(y,w,s,pi,coef,
                                  method=method,
                                  groups=groups)
        return(phi)
    }
}

aggregate_responses <- function(y,w,s,pi,groups){
    y <- split(y,s)
    w <- split(w,s)
    pi <- split(pi,s)
    pi <- split(pi,groups)
    y <- split(y,groups)
    w <- split(w,groups)
    agg <- try(Map(aggregate_responses_uq,y,pi,w),silent=TRUE)
    if(inherits(agg,"try-error")){
        msg <- "Predictions vary within groups"
        if(grepl(msg,agg,fixed=TRUE)) 
            warning(
                paste("While estimating overdispersion -",
                      tolower(msg)),
                call.=FALSE,immediate.=TRUE)
        else stop("Cannot aggregate heterogenous responses")
        agg <- Map(aggregate_responses_av,y,pi,w)
    }
    y <- lapply(agg,"[[","y")
    pi <- lapply(agg,"[[","pi")
    w <- lapply(agg,"[[","w")
    k <- lapply(agg,"[[","k")
    y <- unlist(y)
    pi <- unlist(pi)
    w <- unlist(w)
    k <- unlist(k)
    s <- rep(seq.int(length(k)),k)
    list(y=y,w=w,pi=pi,s=s)
}

aggregate_responses_uq <- function(y,pi,w){
    Pi <- do.call(rbind,pi)
    diffPi <- diff(Pi)
    if(any(diffPi != 0)){
        stop("Predictions vary within groups")
    }
    else 
        pi <- Pi[1,]
    Y <- do.call(rbind,y)
    W <- do.call(rbind,w)
    w <- colSums(W)
    y <- colSums(Y*W)/w
    k <- length(pi)
    list(y=y,pi=pi,w=w,k=k)
}

aggregate_responses_av <- function(y,pi,w){
    Y <- do.call(rbind,y)
    W <- do.call(rbind,w)
    w <- colSums(W)
    y <- colSums(Y*W)/w
    Pi <- do.call(rbind,pi)
    pi <- Pi[1,]
    k <- length(pi)
    list(y=y,pi=pi,w=w,k=k)
}
