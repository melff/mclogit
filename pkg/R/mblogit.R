#' Baseline-Category Logit Models for Categorical and Multinomial Responses
#'
#' The function \code{mblogit} fits baseline-category logit models for categorical
#' and multinomial count responses with fixed alternatives. 
#'
#' @param formula the model formula. The response must be a factor or a matrix
#'     of counts.
#' @param data an optional data frame, list or environment (or object coercible
#'     by \code{\link{as.data.frame}} to a data frame) containing the variables
#'     in the model.  If not found in \code{data}, the variables are taken from
#'     \code{environment(formula)}, typically the environment from which
#'     \code{glm} is called.
#' @param random an optional formula or list of formulas that specify the
#'     random-effects structure or NULL.
#' @param subset an optional vector specifying a subset of observations to be
#'     used in the fitting process.
#' @param weights an optional vector of weights to be used in the fitting
#'     process.  Should be \code{NULL} or a numeric vector.
#' @param na.action a function which indicates what should happen when the data
#'     contain \code{NA}s.  The default is set by the \code{na.action} setting
#'     of \code{\link{options}}, and is \code{\link{na.fail}} if that is unset.
#'     The \sQuote{factory-fresh} default is \code{\link{na.omit}}.  Another
#'     possible value is \code{NULL}, no action.  Value \code{\link{na.exclude}}
#'     can be useful.
#' @param model a logical value indicating whether \emph{model frame} should be
#'     included as a component of the returned value.
#' @param x,y logical values indicating whether the response vector and model
#'     matrix used in the fitting process should be returned as components of
#'     the returned value.
#' @param contrasts an optional list. See the \code{contrasts.arg} of
#'     \code{model.matrix.default}.
#' @param method \code{NULL} or a character string, either "PQL" or "MQL",
#'     specifies the type of the quasilikelihood approximation to be used if a
#'     random-effects model is to be estimated.
#' @param estimator a character string; either "ML" or "REML", specifies which
#'     estimator is to be used/approximated.
#' @param dispersion a logical value or a character string; whether and how a
#'     dispersion parameter should be estimated. For details see
#'     \code{\link{dispersion}}.
#' @param from.table a logical value; do the data represent a contingency table,
#'     e.g. were created by applying \code{as.data.frame()} a the result of
#'     \code{table()} or \code{xtabs()}.  This relevant only if the response is
#'     a factor. This argument should be set to \code{TRUE} if the data do come
#'     from a contingency table. Correctly setting \code{from.table=TRUE} in
#'     this case, will lead to efficiency gains in computing, but more
#'     importantly overdispersion will correctly be computed if present.
#' @param groups an optional formula that specifies groups of observations
#'     relevant for the specification of overdispersed response counts.
#' @param control a list of parameters for the fitting process.  See
#'     \code{\link{mclogit.control}}
#' @param \dots arguments to be passed to \code{mclogit.control} or
#'     \code{mmclogit.control}
#'
#' @return \code{mblogit} returns an object of class "mblogit", which has almost
#'     the same structure as an object of class "\link[stats]{glm}". The
#'     difference are the components \code{coefficients}, \code{residuals},
#'     \code{fitted.values}, \code{linear.predictors}, and \code{y}, which are
#'     matrices with number of columns equal to the number of response
#'     categories minus one.
#' 
#' @details The function \code{mblogit} internally rearranges the data into a
#'     'long' format and uses \code{\link{mclogit.fit}} to compute
#'     estimates. Nevertheless, the 'user data' are unaffected.
#'
#' @seealso The function \code{\link[nnet]{multinom}} in package \pkg{nnet} also
#'     fits multinomial baseline-category logit models, but has a slightly less
#'     convenient output and does not support overdispersion or random
#'     effects. However, it provides some other options. Baseline-category logit
#'     models are also supported by the package \pkg{VGAM}, as well as some
#'     reduced-rank and (semi-parametric) additive generalisations.  The package
#'     \pkg{mnlogit} estimates logit models in a way optimized for large numbers
#'     of alternatives.
#' 
#' @example examples/mblogit-ex.R
#' 
#' @references
#'    Agresti, Alan. 2002.
#'    \emph{Categorical Data Analysis.} 2nd ed, Hoboken, NJ: Wiley.
#'    \url{https://doi.org/10.1002/0471249688}
#'
#'    Breslow, N.E. and D.G. Clayton. 1993.
#'    "Approximate Inference in Generalized Linear Mixed Models".
#'    \emph{Journal of the American Statistical Association} 88 (421): 9-25.
#'    \url{https://doi.org/10.1080/01621459.1993.10594284}
#'
#' 
#' @aliases print.mblogit summary.mblogit print.summary.mblogit fitted.mblogit
#'     weights.mblogit print.mmblogit summary.mmblogit print.summary.mmblogit
mblogit <- function(formula,
                    data=parent.frame(),
                    random=NULL,
                    subset,
                    weights=NULL,
                    na.action = getOption("na.action"),
                    model = TRUE, x = FALSE, y = TRUE,
                    contrasts=NULL,
                    method = NULL,
                    estimator=c("ML","REML"),
                    dispersion = FALSE,
                    from.table = FALSE,
                    groups = NULL,
                    control=if(length(random))
                                mmclogit.control(...)
                            else mclogit.control(...),
                    ...){
    
    call <- match.call(expand.dots = TRUE)
    
    if(missing(data)) data <- environment(formula)
    else if(is.table(data)){
        from.table <- TRUE
        data <- as.data.frame(data)
    }
    else 
        data <- as.data.frame(data)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "offset", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")

    if(length(random)){
        mf0 <- eval(mf, parent.frame())
        mt <- attr(mf0,"terms")
        if(inherits(random,"formula")){
            rf <- paste(c(".~.",all.vars(random)),collapse="+")
        }
        else if(inherits(random,"list")) {
            rf <- paste(c(".~.",unlist(lapply(random,all.vars))),collapse="+")
        }
        else
            stop("'random' argument must be either a formula or a list of formulae")
        rf <- as.formula(rf)
        if (typeof(mf$formula) == "symbol") {
          mff <- formula
        }
        else {
          mff <- structure(mf$formula,class="formula")
        }
        mff <- eval(mff, parent.frame())
        mf$formula <- update(mff,rf)
        mf <- eval(mf, parent.frame())
    }
    else if(length(groups)){
        mf0 <- eval(mf, parent.frame())
        mt <- attr(mf0,"terms")
        gf <- paste(c(".~.",all.vars(groups)),collapse="+")
        gf <- as.formula(gf)
        if (typeof(mf$formula) == "symbol") {
          mff <- formula
        }
        else {
          mff <- structure(mf$formula,class="formula")
        }
        mff <- eval(mff, parent.frame())
        mf$formula <- update(mff,gf)
        mf <- eval(mf, parent.frame())
    }
    else {
        mf <- eval(mf, parent.frame())
        mt <- attr(mf,"terms")
    }
    
    na.action <- attr(mf,"na.action")
    weights <- as.vector(model.weights(mf))
    offset <- as.vector(model.offset(mf))
    if(!is.null(weights) && !is.numeric(weights))
        stop("'weights' must be a numeric vector")
    
    Y <- model.response(mf, "any")
    X <- model.matrix(mt,mf,contrasts)
    contrasts <- attr(X, "contrasts")
    xlevels <- .getXlevels(mt,mf)
    
    if(is.null(weights))
        weights <- rep(1,nrow(X))
    N <- sum(weights)
    prior.weights <- weights

    if(is.factor(Y)){
        response.type <- "factor"
        if(from.table){
            # Create an appropriate response matrix if data
            # come from a table of frequencies
            tmf <- terms(mf)
            respix <- attr(tmf,"response")
            vars <- as.character(attr(tmf,"variables")[-1])
            respname <- vars[respix]
            respix <- match(respname,names(mf))
        
            wghix <- match("(weights)",names(mf))
            mf1 <- mf[-c(respix,wghix)]

            umf1 <- !duplicated(mf1)
            i <- cumsum(umf1)
            j <- as.integer(Y)
            attr(mf,"ij") <- cbind(i,j)
            attr(mf,"j==1") <- umf1
            
            levs <- levels(Y)
            m <- nlevels(Y)
            n <- i[length(i)]

            Y <- matrix(0,nrow=n,ncol=m)
            Y[cbind(i,j)] <- prior.weights
            w <- rowSums(Y)
            Y <- Y/w
            if(any(w==0)){
                Y[w==0,] <- 0
                N <- sum(weights[w>0])
                warning(sprintf("ignoring %d observerations with counts that sum to zero",
                                sum(w==0)),
                        call. = FALSE, immediate. = TRUE)
            }
            Y <- as.vector(t(Y))
            weights <- rep(w,each=m)
            D <- diag(m)[,-1, drop=FALSE]
            dimnames(D) <- list(levs,levs[-1])
            X <- X[umf1,,drop=FALSE]
        }
        else {
            weights <- rep(weights,each=nlevels(Y))
            D <- diag(nlevels(Y))[,-1, drop=FALSE]
            dimnames(D) <- list(levels(Y),levels(Y)[-1])
            I <- diag(nlevels(Y))
            dimnames(I) <- list(levels(Y),levels(Y))
            Y <- as.vector(I[,Y])
        }
    } else if(is.matrix(Y)){
        response.type <- "matrix"
        
        D <- diag(ncol(Y))[,-1, drop=FALSE]
        if(length(colnames(Y))){
            rownames(D) <- colnames(Y)
            colnames(D) <- colnames(Y)[-1]
        }
        else {
            rownames(D) <- 1:ncol(Y)
            colnames(D) <- 2:ncol(Y)
        }
        
        w <- rowSums(Y)
        Y <- Y/w
        if(any(w==0)){
            Y[w==0,] <- 0
            N <- sum(weights[w>0])
            warning(sprintf("ignoring %d observerations with counts that sum to zero",
                            sum(w==0)),
                    call. = FALSE, immediate. = TRUE)
        }
        weights <- rep(w*weights,each=ncol(Y))
        Y <- as.vector(t(Y))
    }
    else stop("response must either be a factor or a matrix of counts or dummies")
    
    s <- rep(seq_len(nrow(X)),each=nrow(D))
    
    XD <- X%x%D

    colnames(XD) <- paste0(rep(colnames(D),ncol(X)),
                                  "~",
                                  rep(colnames(X),each=ncol(D)))

    if(!length(random))
        fit <- mclogit.fit(y=Y,s=s,w=weights,X=XD,
                           dispersion=dispersion,
                           control=control)
    else { ## random effects

        if(!length(method)) method <- "PQL"

        if(inherits(random,"formula"))
            random <- list(random)

        random <- lapply(random,setupRandomFormula)
        rt <- lapply(random,"[[","formula")
        rt <- lapply(rt,terms)
        suppressWarnings(Z <- lapply(rt,model.matrix,mf,
                                     contrasts.arg=contrasts))
        # Use suppressWarnings() to stop complaining about unused contasts
        
        ZD <- lapply(Z,`%x%`,D)
        d <- sapply(ZD,ncol)

        nn <- length(ZD)
        for(k in 1:nn){
            colnames(ZD[[k]]) <- paste0(rep(colnames(D),ncol(Z[[k]])),
                                        "~",
                                        rep(colnames(Z[[k]]),each=ncol(D)))
            colnames(ZD[[k]]) <- gsub("(Intercept)","1",colnames(ZD[[k]]),fixed=TRUE)
        }

        randstruct <- lapply(1:nn,function(k){
            group.labels <- random[[k]]$groups
            groups <- mf[group.labels]
            groups <- lapply(groups,as.factor)
            nlev <- length(groups)
            if(nlev > 1){
                for(i in 2:nlev){
                    groups[[i]] <- interaction(groups[c(i-1,i)])
                    group.labels[i] <- paste(group.labels[i-1],group.labels[i],sep=":")
                }
            }
            groups <- lapply(groups,rep,each=nrow(D))
            
            VarCov.names.k <- rep(list(colnames(ZD[[k]])),nlev)
            ZD_k <- lapply(groups,mkZ,rX=ZD[[k]])
            d <- rep(d[k],nlev)
            names(groups) <- group.labels
            list(ZD_k,groups,d,VarCov.names.k)
        })
        ZD <- lapply(randstruct,`[[`,1)
        groups <- lapply(randstruct,`[[`,2)
        d <- lapply(randstruct,`[[`,3)
        VarCov.names <- lapply(randstruct,`[[`,4)
        ZD <- unlist(ZD,recursive=FALSE)
        groups <- unlist(groups,recursive=FALSE)
        VarCov.names <- unlist(VarCov.names,recursive=FALSE)
        d <- unlist(d)
        ZD <- blockMatrix(ZD,ncol=length(ZD))
        fit <- mmclogit.fitPQLMQL(y=Y,s=s,w=weights,
                                  X=XD,Z=ZD,d=d,
                                  method=method,
                                  estimator=estimator,
                                  control=control,
                                  offset = offset)
        nlev <- length(fit$VarCov)
        for(k in 1:nlev)
            dimnames(fit$VarCov[[k]]) <- list(VarCov.names[[k]],VarCov.names[[k]])
        names(fit$VarCov) <- names(groups)
    }
    
    coefficients <- fit$coefficients
    coefmat <- matrix(coefficients,nrow=ncol(D),
                      dimnames=list("Response categories"=colnames(D),
                                    "Predictors"=colnames(X)
                                    ))
    
    
    fit$coefmat <- coefmat
    fit$coefficients <- coefficients
    
    if(x) fit$x <- X
    if(x && length(random)) fit$z <- Z
    if(!y) {
        fit$y <- NULL
        fit$s <- NULL
    }
    fit <- c(fit,list(call = call, formula = formula,
                      terms = mt,
                      random = random,
                      groups = groups,
                      data = data,
                      contrasts = contrasts,
                      xlevels = xlevels,
                      na.action = na.action,
                      prior.weights=prior.weights,
                      weights=weights,
                      model=mf,
                      D=D,
                      N=N,
                      response.type=response.type,
                      from.table=from.table))

    if(length(random)){
        class(fit) <- c("mmblogit","mblogit","mmclogit","mclogit","lm")
    }
    else
        class(fit) <- c("mblogit","mclogit","lm")
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
  if(x$phi != 1)
      cat("\nDispersion: ",x$phi)
  
  cat("\nNull Deviance:    ",   format(signif(x$null.deviance, digits)),
      "\nResidual Deviance:", format(signif(x$deviance, digits)))
  if(!x$converged) cat("\n\nNote: Algorithm did not converge.\n")
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
    if(x$dispersion != 1)
        cat("\nDispersion: ",x$dispersion," on ",x$df.residual," degrees of freedom")

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
    
    if(!x$converged) cat("\n\nNote: Algorithm did not converge.\n")
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
  ncat <- nrow(object$D)
  fit <- t(matrix(longfit,nrow=ncat))
  
  if(!is.null(na.act))
    fit <- napredict(na.act,fit)
  fit
}


predict.mblogit <- function(object, newdata=NULL,type=c("link","response"),se.fit=FALSE,...){
  
  type <- match.arg(type)
  mt <- terms(object)
  rhs <- delete.response(mt)
  if(missing(newdata)){
    m <- object$model
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
  rn <- rownames(X)
  D <- object$D
  XD <- X%x%D
  rspmat <- function(x){
    y <- t(matrix(x,nrow=nrow(D)))
    colnames(y) <- rownames(D)
    y
  }
  
  eta <- c(XD %*% coef(object))
  eta <- rspmat(eta)
  rownames(eta) <- rn
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
      rownames(se.p) <- rownames(p)
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
    eta <- eta[,-1,drop=FALSE]
    se.eta <- se.eta[,-1,drop=FALSE]
    if(is.null(na.act))
        list(fit=eta,se.fit=se.eta) 
    else
      list(fit=napredict(na.act,eta),
           se.fit=napredict(na.act,se.eta))
  }
  else {
      eta <- eta[,-1,drop=FALSE]
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


print.mmblogit <- function(x,digits= max(3, getOption("digits") - 3), ...){
    cat(paste(deparse(x$call), sep="\n", collapse="\n"), "\n\n", sep="")
  
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

    cat("\n(Co-)Variances:\n")
    VarCov <- x$VarCov
    for(k in 1:length(VarCov)){
        if(k > 1) cat("\n")
        cat("Grouping level:",names(VarCov)[k],"\n")
        VarCov.k <- VarCov[[k]]
        VarCov.k[] <- format(VarCov.k, digits=digits)
        VarCov.k[upper.tri(VarCov.k)] <- ""
        print.default(VarCov.k, print.gap = 2, quote = FALSE)
    }
    
    cat("\nNull Deviance:    ",   format(signif(x$null.deviance, digits)),
        "\nResidual Deviance:", format(signif(x$deviance, digits)))
    if(!x$converged) cat("\n\nNote: Algorithm did not converge.\n")
    if(nchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep="")
    else cat("\n")
    invisible(x)
}

summary.mmblogit <- function(object,...){
  ans <- NextMethod()
  ans$D <- object$D
  class(ans) <- c("summary.mmblogit","summary.mmclogit")
  return(ans)
}

print.summary.mmblogit <-
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
 
    cat("\n(Co-)Variances:\n")
    VarCov <- x$VarCov
    se_VarCov <- x$se_VarCov
    for(k in 1:length(VarCov)){
        if(k > 1) cat("\n")
        cat("Grouping level:",names(VarCov)[k],"\n")
        VarCov.k <- VarCov[[k]]
        VarCov.k[] <- format(VarCov.k, digits=digits)
        VarCov.k[upper.tri(VarCov.k)] <- ""
        #print.default(VarCov.k, print.gap = 2, quote = FALSE)
        VarCov.k <- format_Mat(VarCov.k,title="Estimate")

        se_VarCov.k <- se_VarCov[[k]]
        se_VarCov.k[] <- format(se_VarCov.k, digits=digits)
        se_VarCov.k[upper.tri(se_VarCov.k)] <- ""
        se_VarCov.k <- format_Mat(se_VarCov.k,title="Std.Err.",rownames=" ")

        VarCov.k <- paste(VarCov.k,se_VarCov.k)
        writeLines(VarCov.k)
    }

    cat("\nNull Deviance:    ", format(signif(x$null.deviance, digits)),
        "\nResidual Deviance:", format(signif(x$deviance, digits)),
        "\nNumber of Fisher Scoring iterations: ", x$iter)

    cat("\nNumber of observations")
    for(i in seq_along(x$groups)){
        g <- nlevels(x$groups[[i]])
        nm.group <- names(x$groups)[i]
        cat("\n  Groups by",
            paste0(nm.group,": ",format(g)))
    }
    cat("\n  Individual observations: ",x$N)

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
    
    if(!x$converged) cat("\nNote: Algorithm did not converge.\n")
    if(nchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep="")
    else cat("\n")
    invisible(x)
  }

simulate.mblogit <- function(object, nsim = 1, seed = NULL, ...){

    if(object$phi > 1)
        stop("Simulating responses from models with oversdispersion is not supported yet")

    if(object$response.type=="matrix" || object$from.table){
        yy <- NextMethod()
        seed_attr <- attr(yy,"seed")
        nm <- nrow(yy)
        m <- nrow(object$D)
        n <- nm %/% m
        yy <- lapply(yy,array,
                     dim=c(m,n),
                     dimnames=list(rownames(object$D),
                                   NULL))
        yy <- lapply(yy,t)
        names(yy) <- paste0("sim_",1:nsim)
        
        if(object$response.type=="matrix"){
            class(yy) <- "data.frame"
            attr(yy,"row.names") <- rownames(object$model)
            attr(yy,"seed") <- seed_attr
            return(yy)
        }
        else {
            ij <- attr(object$model,"ij")
            n <- nrow(ij)
            yy <- lapply(yy,"[",ij)
            yy <- as.data.frame(yy)
            attr(yy,"seed") <- seed_attr
            return(yy)
        }
    }
    else { # response.type == "factor"
        probs <- object$fitted.values
        response <- model.response(object$model)
        nm <- length(probs)
        m <- nrow(object$D)
        n <- nm %/% m
        dim(probs) <- c(m,n)
        yy <- sample_factor(probs,nsim=nsim,seed=seed)
        seed_attr <- attr(yy,"seed")
        colnames(yy) <- paste0("sim_",1:nsim)
        rownames(yy) <- rownames(object$model)
        yy <- as.data.frame(yy)
        yy <- lapply(yy,factor,labels=levels(response))
        yy <- as.data.frame(yy)
        attr(yy,"seed") <- seed_attr
        return(yy)
    }
}

simulate.mmblogit <- function(object, nsim = 1, seed = NULL, ...)
    stop("Simulating responses from random-effects models is not supported yet")

sample_factor <- function(probs, nsim =1, seed = NULL, ...){
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
        runif(1)
    if (is.null(seed)) 
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    yy <- apply(probs,2,sample.int,size=nsim,n=nrow(probs),replace=TRUE)
    yy <- t(yy)
    attr(yy,"seed") <- RNGstate
    return(yy)
}

lenuniq <- function(x) length(unique(x))

predict.mmblogit <- function(object, newdata=NULL,type=c("link","response"),se.fit=FALSE,
                             conditional=TRUE, ...){
    
    type <- match.arg(type)
    rhs <- object$formula[-2]
    random <- object$random  
    if(missing(newdata)){
        mf <- object$model
        na.act <- object$na.action
        rmf <- mf
    }
    else{
        mf <- model.frame(rhs,data=newdata,na.action=na.exclude)
        rnd <- object$random
        for(i in seq_along(rnd)){
            rf_i <- random2formula(rnd[[i]])
            if(i == 1)
                rfo <- rf_i
            else
                rfo <- c_formulae(rfo,rf_i)
        }
        rmf <- model.frame(rfo,data=newdata,na.action=na.exclude)
        na.act <- attr(mf,"na.action")
    }
    X <- model.matrix(rhs,mf,
                      contrasts.arg=object$contrasts,
                      xlev=object$xlevels
                      )
    D <- object$D
    XD <- X%x%D
    eta <- c(XD %*% coef(object))

    if(object$method=="PQL" && conditional){

        rf <- lapply(random,"[[","formula")
        rt <- lapply(rf,terms)
        suppressWarnings(Z <- lapply(rt,model.matrix,rmf,
                                     contrasts.arg=object$contrasts,
                                     xlev=object$xlevels))
        
        ZD <- lapply(Z,`%x%`,D)
        d <- sapply(ZD,ncol)

        nn <- length(ZD)
        for(k in 1:nn){
            colnames(ZD[[k]]) <- paste0(rep(colnames(D),ncol(Z[[k]])),
                                        "~",
                                        rep(colnames(Z[[k]]),each=ncol(D)))
            colnames(ZD[[k]]) <- gsub("(Intercept)","1",colnames(ZD[[k]]),fixed=TRUE)
        }

        orig.groups <- object$groups
        olevels <- lapply(orig.groups,levels)
        randstruct <- lapply(1:nn,function(k){
            group.labels <- random[[k]]$groups
            groups <- mf[group.labels]
            groups <- lapply(groups,as.factor)
            nlev <- length(groups)
            if(nlev > 1){
                for(i in 2:nlev){
                    groups[[i]] <- interaction(groups[c(i-1,i)])
                    group.labels[i] <- paste(group.labels[i-1],group.labels[i],sep=":")
                }
            }
            groups <- lapply(groups,rep,each=nrow(D))
            olevels <- olevels[group.labels]
            groups <- Map(factor,x=groups,levels=olevels)
            
            VarCov.names.k <- rep(list(colnames(ZD[[k]])),nlev)
            ZD_k <- lapply(groups,mkZ,rX=ZD[[k]])
            d <- rep(d[k],nlev)
            names(groups) <- group.labels
            list(ZD_k,groups,d,VarCov.names.k)
        })

        ZD <- lapply(randstruct,`[[`,1)
        groups <- lapply(randstruct,`[[`,2)
        ZD <- unlist(ZD,recursive=FALSE)
        d <- lapply(randstruct,`[[`,3)
        groups <- unlist(groups,recursive=FALSE)
        d <- unlist(d)
        
        ZD <- blockMatrix(ZD)
        b <- object$random.effects
        nlev <- length(ZD)
        
        for(k in 1:nlev)
            eta <- eta +  as.vector(ZD[[k]]%*%b[[k]])
    }
    
    rspmat <- function(x){
        y <- t(matrix(x,nrow=nrow(D)))
        colnames(y) <- rownames(D)
        y
    }
    eta <- rspmat(eta)

    nvar <- ncol(X)
    nobs <- nrow(X)
    
    if(se.fit || type=="response"){
        exp.eta <- exp(eta)
        sum.exp.eta <- rowSums(exp.eta)
        p <- exp.eta/sum.exp.eta
    }
    if(se.fit){
        ncat <- ncol(p)
        W <- Matrix(0,nrow=nobs*ncat,ncol=nobs)
        i <- seq.int(ncat*nobs)
        j <- rep(1:nobs,each=ncat)
        pv <- as.vector(t(p))
        W[cbind(i,j)] <- pv
        W <- Diagonal(x=pv)-tcrossprod(W)
        WX <- W%*%XD
        if(object$method=="PQL"){
            WZ <- bMatProd(W,ZD)
            H <- object$info.fixed.random
            K <- solve(H)
        }
    }
    
    if(type=="response") {
        if(se.fit){
            if(object$method=="PQL" && conditional){
                WXZ <- structure(cbind(blockMatrix(WX),WZ),class="blockMatrix")
                var.p <- bMatProd(WXZ,K)
                var.p <- Map(`*`,WXZ,var.p)
                var.p <- lapply(var.p,rowSums)
                var.p <- Reduce(`+`,var.p)
            }
            else {
                vcov.coef <- vcov(object)
                var.p <- rowSums(WX*(WX%*%vcov.coef))
            }
            se.p <- sqrt(var.p)
            se.p <- rspmat(se.p)
            if(is.null(na.act))
                list(fit=p,se.fit=se.p) 
            else
                list(fit=napredict(na.act,p),
                     se.fit=napredict(na.act,se.p))
        }
        else{
            if(is.null(na.act)) p
            else napredict(na.act,p)
        }
    }
    else {
        eta <- eta[,-1,drop=FALSE]
        if(se.fit){
            if(object$method=="PQL" && conditional){
                XZ <- structure(cbind(blockMatrix(XD),ZD),class="blockMatrix")
                var.eta <- bMatProd(XZ,K)
                var.eta <- Map(`*`,XZ,var.eta)
                var.eta <- lapply(var.eta,rowSums)
                var.eta <- Reduce(`+`,var.eta)
            }
            else {
                vcov.coef <- vcov(object)
                var.eta <- rowSums(XD*(XD%*%vcov.coef))
            }
            se.eta <- sqrt(var.eta)
            se.eta <- rspmat(se.eta)
            se.eta <- se.eta[,-1,drop=FALSE]
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
}
