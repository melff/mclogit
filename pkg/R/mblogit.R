mblogit <- function(formula,
                    data=parent.frame(),
                    subset,
                    weights=NULL,
                    na.action = getOption("na.action"),
                    model = TRUE, x = FALSE, y = TRUE,
                    contrasts=NULL,
                    start.theta=NULL,
                    control=mclogit.control(...),
                    ...){
  
  call <- match.call(expand.dots = TRUE)
  
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
  
  y <- model.response(mf, "any")
  X <- model.matrix(mt,mf,contrasts)
  contrasts <- attr(X, "contrasts")
  xlevels <- .getXlevels(mt,mf)
  
  if(is.null(weights))
    weights <- rep(1,nrow(X))
  N <- sum(weights)
  prior.weights <- weights
  
  if(is.factor(y)){
    
    D <- contrasts(y)
    I <- diag(nlevels(y))
    dimnames(I) <- list(levels(y),levels(y))
    yy <- c(I[,y])
    
    weights <- rep(weights,each=nlevels(y))
  } else if(is.matrix(y)){
    
    D <- diag(ncol(y))[,-1]
    if(length(colnames(y))){
      rownames(D) <- colnames(y)
      colnames(D) <- colnames(y)[-1]
    }
    else {
      rownames(D) <- 1:ncol(y)
      colnames(D) <- 2:ncol(y)
    }
    
    w <- rowSums(y)
    yy <- c(t(y/w))
    weights <- rep(w*weights,each=ncol(y))
  }
  s <- rep(seq_len(nrow(X)),each=nrow(D))
  
  XD <- X%x%D

  fit <- mclogit.fit(y=yy,s=s,w=weights,X=XD,
                     control=control)
  
  coefficients <- fit$coefficients
  coefmat <- matrix(coefficients,nrow=ncol(D),
                         dimnames=list("Response categories"=colnames(D),
                                       "Predictors"=colnames(X)
                         ))
  names(coefficients) <- paste0(rep(colnames(D),ncol(X)),
                                "~",
                                rep(colnames(X),each=ncol(D)))
  
  
  fit$coefmat <- coefmat
  fit$coefficients <- coefficients
  
  fit <- c(fit,list(call = call, formula = formula,
                    terms = mt,
                    random = NULL,
                    data = data,
                    contrasts = contrasts,
                    xlevels = xlevels,
                    na.action = na.action,
                    prior.weights=prior.weights,
                    weights=weights,
                    model=mf,
                    D=D,
                    N=N))
  class(fit) <- c("mblogit","mclogit")
  fit
}


print.mblogit <- function(x,digits= max(3, getOption("digits") - 3), ...){
  cat("\nCall: ", deparse(x$call), "\n\n")
  
  D <- x$D
  
  categs <- colnames(D)
  basecat <- rownames(D)[!(rownames(D)%in%categs)]
  
  coefmat <- x$coefmat
  if(getOption("mblogit.show.basecat",TRUE)){
    
    rn <- paste0(rownames(coefmat), getOption("mblogit.basecat.sep","/"), basecat)
    rownames(coefmat) <- rn
  }
  
  if(length(coefmat)) {
    cat("Coefficients")
    if(is.character(co <- x$contrasts))
      cat("  [contrasts: ",
          apply(cbind(names(co),co), 1, paste, collapse="="), "]")
    cat(":\n")
    print.default(format(coefmat, digits=digits),
                  print.gap = 2, quote = FALSE)
  } else cat("No coefficients\n\n")
  
  cat("\nNull Deviance:    ",   format(signif(x$null.deviance, digits)),
      "\nResidual Deviance:", format(signif(x$deviance, digits)))
  if(nchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep="")
  else cat("\n")
  invisible(x)
}


summary.mblogit <- function(object,...){
  ans <- NextMethod()
  ans$D <- object$D
  class(ans) <- c("summary.mblogit","summary.mclogit")
  return(ans)
}
  
print.summary.mblogit <-
  function (x, digits = max(3, getOption("digits") - 3),
            symbolic.cor = x$symbolic.cor,
            signif.stars = getOption("show.signif.stars"), ...){
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep="\n", collapse="\n"), "\n\n", sep="")

    D <- x$D
    
    categs <- colnames(D)
    basecat <- rownames(D)[!(rownames(D)%in%categs)]
    
    coefs <- x$coefficients
    rn.coefs <- rownames(coefs)
    ncategs <- length(categs)
    
    for(i in 1:ncategs){
      cat <- categs[i]
      patn <- paste0(cat,"~")
      ii <- grep(patn,rn.coefs,fixed=TRUE)
      coefs.cat <- coefs[ii,,drop=FALSE]
      rownames(coefs.cat) <- gsub(patn,"",rownames(coefs.cat))
      if(i>1) cat("\n")
      cat("Equation for ",cat," vs ",basecat,":\n",sep="")
      printCoefmat(coefs.cat, digits=digits, signif.stars=signif.stars,
                   signif.legend=signif.stars && i==ncategs,
                   na.print="NA", ...)
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
    
    if(nchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep="")
    else cat("\n")
    invisible(x)
  }


fitted.mblogit <- function(object,type=c("probabilities","counts"),...){
  weights <- object$weights
  nobs <- length(weights)
  res <- object$fitted.values
  type <- match.arg(type)
  
  na.act <- object$na.action
  
  longfit <- switch(type,
                probabilities=res,
                counts=weights*res)
  
  if(!is.null(na.act))
    longfit <- napredict(na.act,longfit)
  
  ncat <- nrow(object$D)
  fit <- t(matrix(longfit,nrow=ncat))
  fit
}


predict.mblogit <- function(object, newdata=NULL,type=c("link","response"),se.fit=FALSE,...){
  
  type <- match.arg(type)
  rhs <- object$formula[-2]
  if(missing(newdata)){
    m <- model.frame(rhs,data=object$model)
    na.act <- object$na.action
  }
  else{
    m <- model.frame(rhs,data=newdata,na.action=na.exclude)
    na.act <- attr(m,"na.action")
  }
  X <- model.matrix(rhs,m,
                    contrasts.arg=object$contrasts,
                    xlev=object$xlevels
  )
  D <- object$D
  XD <- X%x%D
  rspmat <- function(x){
    y <- t(matrix(x,nrow=nrow(D)))
    colnames(y) <- rownames(D)
    y
  }
  
  eta <- c(XD %*% coef(object))
  eta <- rspmat(eta)
  if(se.fit){
    V <- vcov(object)
    stopifnot(ncol(XD)==ncol(V))
  }
  
  
  
  if(type=="response") {
    exp.eta <- exp(eta)
    sum.exp.eta <- rowSums(exp.eta)
    p <- exp.eta/sum.exp.eta
    
    if(se.fit){
      p.long <- as.vector(t(p))
      s <- rep(1:nrow(X),each=nrow(D))
      
      wX <- p.long*(XD - rowsum(p.long*XD,s)[s,,drop=FALSE])
      se.p.long <- sqrt(rowSums(wX * (wX %*% V)))
      se.p <- rspmat(se.p.long)
      
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
    se.eta <- sqrt(rowSums(XD * (XD %*% V)))
    se.eta <- rspmat(se.eta)
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


weights.mblogit <- function (object, ...) 
{
  res <- object$prior.weights
  if (is.null(object$na.action)) 
    res
  else naresid(object$na.action, res)
}