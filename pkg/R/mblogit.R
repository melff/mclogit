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
#' @param catCov a character string that specifies optional restrictions
#'     on the covariances of random effects between the logit equations.
#'     "free" means no restrictions, "diagonal" means that random effects
#'     pertinent to different categories are uncorrelated, while "single" means
#'     that the random effect variances pertinent to all categories are identical.
#' @param subset an optional vector specifying a subset of observations to be
#'     used in the fitting process.
#' @param weights an optional vector of weights to be used in the fitting
#'     process.  Should be \code{NULL} or a numeric vector.
#' @param offset an optional model offset. If not NULL, must be a matrix
#'     if as many columns as the response has categories or one less.
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
#' @param start an optional matrix of starting values (with as many rows
#'     as logit equations). If the model has random effects, the matrix
#'     should have a "VarCov" attribute wtih starting values for
#'     the random effects (co-)variances. If the random effects model
#'     is estimated with the "PQL" method, the starting values matrix
#'     should also have a "random.effects" attribute, which should have
#'     the same structure as the "random.effects" component of an object
#'     returned by \code{mblogit()}.
#' @param aggregate a logical value; whether to aggregate responses by
#'     covariate classes and groups before estimating the model
#'     if the response variable is a factor. 
#'     
#'     This will not affect the estimates, but the dispersion and the
#'     residual degrees of freedom. If \code{aggregate=TRUE}, the 
#'     dispersion will be relative to a saturated model; it will be much 
#'     smaller than with \code{aggregate=TRUE}. In particular, with only
#'     a single covariate and no grouping, the deviance will be close to
#'     zero. If \code{dispersion} is not \code{FALSE}, then the
#'     default value of \code{aggregate} will be \code{TRUE}. For details see
#'     \code{\link{dispersion}}.
#'     
#'     This argument has consequences only if the response in \code{formula}
#'     is a factor.
#' @param groups an optional formula that specifies groups of observations
#'     relevant for the estimation of overdispersion. For details see
#'     \code{\link{dispersion}}.
#' @param from.table a logical value; should be FALSE. This argument
#'     only exists for the sake of compatibility and will be removed
#'     in the next relase.
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
#'    \doi{10.1002/0471249688}
#'
#'    Breslow, N.E. and D.G. Clayton. 1993.
#'    "Approximate Inference in Generalized Linear Mixed Models".
#'    \emph{Journal of the American Statistical Association} 88 (421): 9-25.
#'    \doi{10.1080/01621459.1993.10594284}
#'
#' 
#' @aliases print.mblogit summary.mblogit print.summary.mblogit fitted.mblogit
#'     weights.mblogit print.mmblogit summary.mmblogit print.summary.mmblogit
mblogit <- function(formula,
                    data=parent.frame(),
                    random=NULL,
                    catCov=c("free","diagonal","single"),
                    subset,
                    weights=NULL,
                    offset=NULL,
                    na.action = getOption("na.action"),
                    model = TRUE, x = FALSE, y = TRUE,
                    contrasts=NULL,
                    method = NULL,
                    estimator=c("ML","REML"),
                    dispersion = FALSE,
                    start = NULL,
                    aggregate = FALSE,
                    groups = NULL,
                    from.table = FALSE,
                    control=if(length(random))
                                mmclogit.control(...)
                            else mclogit.control(...),
                    ...){

    call <- match.call(expand.dots = TRUE)
    
    if(!missing(from.table)) {
        warning("Argument 'from.table' is deprecated. Use 'aggregate=TRUE' instead.")
        if(missing(aggregate))
            aggregate <- from.table
    }
    if(!aggregate) {
        if(length(groups))
            warning("Argument 'groups' is inconsequential unless aggregate=TRUE")
    }

    if(missing(data)) data <- environment(formula)
    else if(is.table(data)){
        from.table <- TRUE
        data <- as.data.frame(data)
    }
    else 
        data <- as.data.frame(data)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0)
    if("offset" %in% names(mf)) {
        offset <- eval(mf$offset,data,environment(formula))
    }
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
        check.names(control,
                    "epsilon","maxit",
                    "trace","trace.inner",
                    "avoid.increase",
                    "break.on.increase",
                    "break.on.infinite",
                    "break.on.negative")
        catCov <- match.arg(catCov)
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
        groups <- all.vars(groups)
        groups <- mf[groups]
        # if(length(groups) > 1) stop("Multiple groups not supported")
        check.names(control,
                    "epsilon","maxit",
                    "trace"
                    )
    }
    else {
        mf <- eval(mf, parent.frame())
        mt <- attr(mf,"terms")
        check.names(control,
                    "epsilon","maxit",
                    "trace")
    }
    
    na.action <- attr(mf,"na.action")
    weights <- as.vector(model.weights(mf))
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

    if(is.factor(Y)) {
        n.categs <- nlevels(Y)
        n.obs <- length(Y)
    } else if(is.matrix(Y)) {
        n.categs <- ncol(Y)
        n.obs <- nrow(Y)
    } else {
        stop("Response must be either a factor or a matrix of counts")
    }
    if(length(offset)) {
        if(!is.matrix(offset)) {
            if(length(offset) != n.obs) stop("'offset' has wrong length")
            offset <- matrix(offset,ncol=n.categs-1) 
            offset <- cbind(0, offset)
        } else {
            if(nrow(offset) != n.obs) 
                stop("'offset' has wrong number of rows")
            if(ncol(offset) != n.categs) {
                if(ncol(offset) != n.categs - 1)
                    stop(sprintf("'offset' must either have %d or %d columns",
                                n.categs-1, n.categs))
                offset <- cbind(0, offset)
            }
        }
    }

    if(is.factor(Y)){
        response.type <- "factor"
        if(!isFALSE(dispersion)) {
            aggregate <- TRUE
        }
        if(aggregate && !length(random)) {
            if(!length(random)) {
                D <- structure(diag(n.categs),
                            dimnames=rep(list(levels(Y)),2))[,-1, drop=FALSE]
                tmf <- terms(mf)
                respix <- attr(tmf,"response")
                vars <- as.character(attr(tmf,"variables")[-1])
                respname <- vars[respix]
                respix <- match(respname,names(mf),nomatch=0L)
            
                wghix <- match("(weights)",names(mf),nomatch=0L)
                mf1 <- mf[-c(respix,wghix)]

                strata <- quickInteraction(mf1)
                
                weights.tab <- rowsum(weights,
                                    quickInteraction(list(Y,strata)))
                dim(weights.tab) <- c(n.categs,attr(strata,"n"))
                w <- colSums(weights.tab)
                weights <- rep(w,each=n.categs)
                Y <- as.vector(weights.tab/weights)
                keep <- !duplicated(strata)
                X <- X[keep,,drop=FALSE]
                if(is.matrix(offset))
                    offset <- offset[keep,,drop=FALSE]
            } else {
                stop("'aggregate' is not yet supported with random effects")
            }
        } else {
            weights <- rep(weights,each=nlevels(Y))
            D <- diag(nlevels(Y))[,-1, drop=FALSE]
            dimnames(D) <- list(levels(Y),levels(Y)[-1])
            I <- diag(nlevels(Y))
            dimnames(I) <- list(levels(Y),levels(Y))
            Y <- as.vector(I[,Y])
        }
    } else if(is.matrix(Y)){
        response.type <- "matrix"
        n.categs <- ncol(Y)
        n.obs <- nrow(Y)

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

    start.VarCov <- NULL
    start.randeff <- NULL
    if(length(start)){
        start.VarCov <- attr(start,"VarCov")
        start.randeff <- attr(start,"random.effects")
        if(nrow(start)!=ncol(D))
            stop("Rows of 'start' argument do not match dependent variable.")
        start.names <- colnames(start)
        X.names <- colnames(X)
        if(length(start.names))
            start <- start[,X.names,drop=FALSE]
        if(ncol(start)!=ncol(X))
             stop("Columns of 'start' argument do not match independent variables.")
        start <- as.vector(start)
    }
    
    s <- rep(seq_len(nrow(X)),each=nrow(D))
    
    XD <- X%x%D

    colnames(XD) <- paste0(rep(colnames(D),ncol(X)),
                                  "~",
                                  rep(colnames(X),each=ncol(D)))

    if(is.matrix(offset)){
        Offset <- offset
        offset <- as.vector(t(offset))
    } else {
        Offset <- NULL
    }

    if(!length(random)){
        fit <- mclogit.fit(y=Y,s=s,w=weights,X=XD,
                           dispersion=dispersion,
                           control=control,
                           start=start,
                           offset = offset)
    }
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

        if(catCov == "free"){
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
        } else if(catCov =="single"){
            cc <- rep(1:n.categs,n.obs)
            stopifnot(length(Y)==length(cc))
            d <- sapply(Z,ncol)

            nn <- length(Z)

            for(k in 1:nn){
                colnames(Z[[k]]) <- paste0("~",colnames(Z[[k]]))
                colnames(Z[[k]]) <- gsub("(Intercept)","1",colnames(Z[[k]]),fixed=TRUE)
            }

            randstruct <- lapply(1:nn,function(k){
                group.labels <- random[[k]]$groups
                groups <- mf[group.labels]
                groups <- lapply(groups,as.factor)
                nlev <- length(groups)
                groups[[1]] <- interaction(cc,groups[[1]])
                if(nlev > 1){
                    for(i in 2:nlev){
                        groups[[i]] <- interaction(groups[c(i-1,i)])
                        group.labels[i] <- paste(group.labels[i-1],group.labels[i],sep=":")
                    }
                }
                
                VarCov.names.k <- rep(list(colnames(Z[[k]])),nlev)
                ZD_k <- lapply(groups,mkZ,rX=Z[[k]])
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
        } else { # catCov == "diagonal"
            categs <- 1:n.categs
            cc <- rep(categs,n.obs)
            stopifnot(length(Y)==length(cc))
            randstruct <- list()
            for(categ in categs){
                u <- as.integer(categ==categs)

                ZD <- lapply(Z,`%x%`,u)
                d <- sapply(ZD,ncol)

                nn <- length(ZD)

                for(k in 1:nn){
                    colnames(ZD[[k]]) <- paste0(rownames(D)[categ],"~",colnames(Z[[k]]))
                    colnames(ZD[[k]]) <- gsub("(Intercept)","1",colnames(ZD[[k]]),fixed=TRUE)
                }

                randstruct_c <- lapply(1:nn,function(k){
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
                randstruct <- c(randstruct,randstruct_c)
            }
            ZD <- lapply(randstruct,`[[`,1)
            groups <- lapply(randstruct,`[[`,2)
            d <- lapply(randstruct,`[[`,3)
            VarCov.names <- lapply(randstruct,`[[`,4)
            ZD <- unlist(ZD,recursive=FALSE)
            groups <- unlist(groups,recursive=FALSE)
            VarCov.names <- unlist(VarCov.names,recursive=FALSE)
            d <- unlist(d)
            ZD <- blockMatrix(ZD,ncol=length(ZD))
        }

        fit <- mmclogit.fitPQLMQL(y=Y,s=s,w=weights,
                                  X=XD,Z=ZD,d=d,
                                  start=start,
                                  start.Phi=start.VarCov,
                                  start.b=start.randeff,
                                  method=method,
                                  estimator=estimator,
                                  control=control,
                                  offset = offset)
        nlev <- length(fit$VarCov)
        for(k in 1:nlev)
            dimnames(fit$VarCov[[k]]) <- list(VarCov.names[[k]],VarCov.names[[k]])
        names(fit$VarCov) <- names(groups)
    }
    fit$offset <- Offset

    coefficients <- fit$coefficients
    coefmat <- matrix(coefficients,nrow=ncol(D),
                      dimnames=list("Logit eqn."=colnames(D),
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
                      start = start,
                      prior.weights=prior.weights,
                      weights=weights,
                      model=mf,
                      D=D,
                      N=N,
                      response.type=response.type,
                      aggregated = aggregate,
                      catCov = catCov))

    if(length(random)){
        class(fit) <- c("mmblogit","mblogit","mmclogit","mclogit","lm")
    }
    else
        class(fit) <- c("mblogit","mclogit","lm")
    fit
}


print.mblogit <- function(x,digits= max(3, getOption("digits") - 3), ...){
  cat("\nCall: ",paste(deparse(x$call), sep="\n", collapse="\n"), "\n\n", sep="")
  
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
      ii <- which(startsWith(rn.coefs,patn))
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

    cat("\nApproximate residual Deviance:", format(signif(x$deviance, digits)),
        "\nNumber of Fisher scoring iterations: ", x$iter,
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
    offset <- object$offset
  }
  else{
    m <- model.frame(rhs,data=newdata,na.action=na.exclude)
    na.act <- attr(m,"na.action")
    offset <- model.offset(m)
    offset_in_call <- object$call$offset
    if(!is.null(offset_in_call)){
        offset_in_call <- eval(offset_in_call,newdata,
                               environment(terms(object)))
        if(length(offset))
            offset <- offset + offset_in_call
        else
            offset <- offset_in_call
    }
  }
  X <- model.matrix(rhs,m,
                    contrasts.arg=object$contrasts,
                    xlev=object$xlevels
  )
  rn <- rownames(X)
  D <- object$D
  n.obs <- nrow(X)
  n.categs <- nrow(D)
  XD <- X%x%D
  eta <- c(XD %*% coef(object))
  
  if(length(offset)) {
      if(!is.matrix(offset)) {
          if(length(offset) != n.obs) stop("'offset' has wrong length")
          offset <- matrix(offset,ncol=n.categs-1) 
          offset <- cbind(0, offset)
      } else {
          if(nrow(offset) != n.obs) 
              stop("'offset' has wrong number of rows")
          if(ncol(offset) != n.categs) {
              if(ncol(offset) != n.categs - 1)
                  stop(sprintf("'offset' must either have %d or %d columns",
                              n.categs-1, n.categs))
              offset <- cbind(0, offset)
          }
      }
      offset <- as.vector(t(offset))
      eta <- eta + offset
  }

  rspmat <- function(x){
    y <- t(matrix(x,nrow=nrow(D)))
    colnames(y) <- rownames(D)
    y
  }
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

format_VarCov <- function(x, digits = 3){
    x <- format(x, digits = digits)
    x[upper.tri(x)] <- ""
    return(x)
}

rcomb <- function(x){
    total.cnames <- unique(unlist(lapply(x,colnames)))
    total.ncol <- length(total.cnames)
    res <- matrix(nrow=0,ncol=total.ncol)
    for(i in seq_along(x)){
        x.i <- x[[i]]
        res.i <- matrix("",nrow=nrow(x.i),ncol=total.ncol)
        res.i[,match(colnames(x.i),total.cnames)] <- x.i
        res <- rbind(res,res.i)
    }
    total.rnames <- unlist(lapply(x,rownames))
    colnames(res) <- total.cnames
    rownames(res) <- total.rnames
    return(res)
}

VC_colnames_drop_lhs <- function(x){
    coln <- colnames(x)
    coln <- strsplit(coln,"~",fixed=TRUE)
    coln <-  unlist(lapply(coln,"[",2))
    colnames(x) <- paste0(".~",coln)
    return(x)
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
    nVC <- names(VarCov)
    unVC <- unique(nVC)
    for(nm in unVC){
        cat("\nGrouping level:",nm,"\n")
        k <- which(nVC==nm)
        VarCov.nm <- VarCov[k]
        if(length(VarCov.nm) == 1){
            VarCov.nm <- format_VarCov(VarCov.nm[[1]], digits = digits)
            print.default(VarCov.nm, print.gap = 2, quote = FALSE)
        }
        else {
            VarCov.nm <- lapply(VarCov.nm, format_VarCov, digits = digits)
            VarCov.nm <- lapply(VarCov.nm,VC_colnames_drop_lhs)
            VarCov.nm <- rcomb(VarCov.nm)
            print.default(VarCov.nm, print.gap = 2, quote = FALSE)
        }
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
      patn <- paste0("^",cat,"~")
      ii <- grep(patn,rn.coefs)
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
    nVC <- names(VarCov)
    unVC <- unique(nVC)
    for(nm in unVC){
        cat("\nGrouping level:",nm,"\n")
        k <- which(nVC==nm)
        VarCov.nm <- VarCov[k]
        if(length(VarCov.nm) == 1){
            VarCov.k <- format_VarCov(VarCov[[k]], digits = digits)
            VarCov.k <- format_Mat(VarCov.k,title="Estimate")
            se_VarCov.k <- se_VarCov[[k]]
            se_VarCov.k <- format_VarCov(se_VarCov[[k]], digits = digits)
            se_VarCov.k <- format_Mat(se_VarCov.k,title="Std.Err.")
            VarCov.k <- paste(VarCov.k,se_VarCov.k)
            writeLines(VarCov.k)
        }
        else {
            VarCov.nm <- lapply(VarCov.nm, format_VarCov, digits = digits)
            VarCov.nm <- lapply(VarCov.nm,VC_colnames_drop_lhs)
            VarCov.nm <- rcomb(VarCov.nm)
            VarCov.nm <- format_Mat(VarCov.nm,title="Estimate")
            se_VarCov.nm <- se_VarCov[k]
            se_VarCov.nm <- lapply(se_VarCov.nm, format_VarCov, digits = digits)
            se_VarCov.nm <- lapply(se_VarCov.nm,VC_colnames_drop_lhs)
            se_VarCov.nm <- rcomb(se_VarCov.nm)
            se_VarCov.nm <- format_Mat(se_VarCov.nm,title="Std.Err.")
            VarCov.nm <- paste(VarCov.nm,se_VarCov.nm)
            writeLines(VarCov.nm)
        }
    }

    cat("\nApproximate residual deviance:", format(signif(x$deviance, digits)),
        "\nNumber of Fisher scoring iterations: ", x$iter)

    cat("\nNumber of observations")
    nm_grps <- names(x$groups)
    unm_grps <- unique(nm_grps)
    for(nm in unm_grps){
        k <- which(nm_grps == nm)
        grps_k <- x$groups[k]
        g <- nlevels(grps_k[[1]])
        cat("\n  Groups by",
            paste0(nm,": ",format(g)))
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

    if(object$response.type=="matrix" || object$aggregated){
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
    mt <- terms(object)
    rhs <- delete.response(mt)
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
        
        catCov <- object$catCov
        if(!length(catCov)) catCov <- "free"
        n.categs <- nrow(D)

        orig.groups <- object$groups
        olevels <- lapply(orig.groups,levels)

        if(catCov == "free"){
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
                groups <- rmf[group.labels]
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
            d <- lapply(randstruct,`[[`,3)
            VarCov.names <- lapply(randstruct,`[[`,4)
            ZD <- unlist(ZD,recursive=FALSE)
            groups <- unlist(groups,recursive=FALSE)
            VarCov.names <- unlist(VarCov.names,recursive=FALSE)
            d <- unlist(d)
            ZD <- blockMatrix(ZD,ncol=length(ZD))
        } else if(catCov =="single"){
            n.obs <- nrow(X)
            cc <- rep(1:n.categs,n.obs)
            d <- sapply(Z,ncol)

            nn <- length(Z)

            for(k in 1:nn){
                colnames(Z[[k]]) <- paste0("~",colnames(Z[[k]]))
                colnames(Z[[k]]) <- gsub("(Intercept)","1",colnames(Z[[k]]),fixed=TRUE)
            }

            randstruct <- lapply(1:nn,function(k){
                group.labels <- random[[k]]$groups
                groups <- mf[group.labels]
                groups <- lapply(groups,as.factor)
                nlev <- length(groups)
                groups[[1]] <- interaction(cc,groups[[1]])
                if(nlev > 1){
                    for(i in 2:nlev){
                        groups[[i]] <- interaction(groups[c(i-1,i)])
                        group.labels[i] <- paste(group.labels[i-1],group.labels[i],sep=":")
                    }
                }
                groups <- lapply(groups,rep,each=nrow(D))
                olevels <- olevels[group.labels]
                groups <- Map(factor,x=groups,levels=olevels)
                
                VarCov.names.k <- rep(list(colnames(Z[[k]])),nlev)
                ZD_k <- lapply(groups,mkZ,rX=Z[[k]])
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
        } else { # catCov == "diagonal"
            categs <- 1:n.categs
            randstruct <- list()
            for(categ in categs){
                u <- as.integer(categ==categs)

                ZD <- lapply(Z,`%x%`,u)
                d <- sapply(ZD,ncol)

                nn <- length(ZD)

                for(k in 1:nn){
                    colnames(ZD[[k]]) <- paste0(rownames(D)[categ],"~",colnames(Z[[k]]))
                    colnames(ZD[[k]]) <- gsub("(Intercept)","1",colnames(ZD[[k]]),fixed=TRUE)
                }

                randstruct_c <- lapply(1:nn,function(k){
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
                randstruct <- c(randstruct,randstruct_c)
            }
            ZD <- lapply(randstruct,`[[`,1)
            groups <- lapply(randstruct,`[[`,2)
            d <- lapply(randstruct,`[[`,3)
            VarCov.names <- lapply(randstruct,`[[`,4)
            ZD <- unlist(ZD,recursive=FALSE)
            groups <- unlist(groups,recursive=FALSE)
            VarCov.names <- unlist(VarCov.names,recursive=FALSE)
            d <- unlist(d)
            ZD <- blockMatrix(ZD,ncol=length(ZD))
        }

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
            H <- object$info.fixed.random
            K <- solve(H)
        }
    }
    
    if(type=="response") {
        if(se.fit){
            if(object$method=="PQL" && conditional){
                WZ <- bMatProd(W,ZD)
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
