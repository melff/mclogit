quickInteraction <- function(by){
  if(is.list(by)){
    n.arg <- length(by)
    f <- 0L
    uf <- 0L
    for(i in rev(1:n.arg)){
      y <- by[[i]]
      y <- as.numeric(y)
      uy <- unique(na.omit(y))
      y <- match(y,uy,NA)
      l <- length(uy)
      f <- f*l + y - 1
      uf <- unique(na.omit(f))
      f <- match(f,uf,NA)
      uf <- seq(length(uf))
    }
  }
  else {
    by <- as.numeric(by)
    uf <- unique(na.omit(by))
    f <- match(by,uf,NA)
    uf <- seq(length(uf))
  }
  return(structure(f,unique=uf))
}

matConstInSets <- function(X,sets){
    ans <- logical(ncol(X))
    for(i in 1:ncol(X)){
        v <- tapply(X[,i],sets,varies)
        ans[i] <- !any(v) 
    }
    ans
}

listConstInSets <- function(X,sets){
    ans <- logical(length(X))
    for(i in 1:length(X)){
        v <- tapply(X[[i]],sets,varies)
        ans[i] <- !any(v) 
    }
    ans
}

groupConstInSets <- function(X,sets){
    ans <- logical(length(X))
    for(i in 1:length(X)){
        v <- tapply(X[[i]],sets,varies)
        ans[i] <- !any(v) 
    }
    ans
}

varies <- function(x)
    !all(duplicated(x)[-1L])


mclogit <- function(
                formula,
                data=parent.frame(),
                random=NULL,
                subset,
                weights=NULL,
                offset=NULL,
                na.action = getOption("na.action"),
                model = TRUE, x = FALSE, y = TRUE,
                contrasts=NULL,
                method = NULL,
                estimator=c("ML","REML"),
                dispersion = FALSE,
                start=NULL,
                control=if(length(random))
                            mmclogit.control(...)
                        else mclogit.control(...),
                ...
            ){
# Assumptions:
#   left hand side of formula: cbind(counts,  choice set index)
#   right hand side of the formula: attributes
#   intercepts are removed!

    call <- match.call(expand.dots = TRUE)

    if(missing(data)) data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "offset", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    if(as.character(formula[[2]][[1]])=="|")
        mf$formula[[2]][[1]] <- as.name("cbind")

    if(length(random)){
        mf0 <- eval(mf, parent.frame())
        mt <- attr(mf0,"terms")
        if(is_formula(random)){
            rf <- paste(c(".~.",all.vars(random)),collapse="+")
        }
        else if(is.list(random)) {
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
    else {
        mf <- eval(mf, parent.frame())
        mt <- attr(mf,"terms")
    }
    
    na.action <- attr(mf,"na.action")
    weights <- as.vector(model.weights(mf))
    offset <- as.vector(model.offset(mf))
    if(!is.null(weights) && !is.numeric(weights))
        stop("'weights' must be a numeric vector")

    Y <- as.matrix(model.response(mf, "any"))
    if(!is.numeric(Y)) stop("The response matrix has to be numeric.")
    if(ncol(Y)<2) stop("need response counts and choice set indicators")
    sets <- Y[,2]
    sets <- match(sets,unique(sets))
    Y <- Y[,1]
    if (is.null(weights)){
        prior.weights <- rep(1,length(Y))
        N <- rowsum(Y,sets,na.rm=TRUE)
        weights <- N[sets]
        }
    else{
        prior.weights <- weights
        N <- rowsum(weights*Y,sets,na.rm=TRUE)
        weights <- N[sets]
        }
    N <- sum(N)
    Y <- Y/weights
    Y[weights==0] <- 0
    
    X <- model.matrix(mt,mf,contrasts)
    contrasts <- attr(X, "contrasts")
    xlevels <- .getXlevels(mt,mf)
    icpt <- match("(Intercept)",colnames(X),nomatch=0)
    if(icpt) X <- X[,-icpt,drop=FALSE]
    const <- matConstInSets(X,sets)
    if(any(const)){
        warning("removing ",
                gsub("(Intercept)","intercept",paste(colnames(X)[const],collapse=","),fixed=TRUE),
                " from model due to insufficient within-choice set variance")
        X <- X[,!const,drop=FALSE]
    }
    if(!length(start)){
      drop.coefs <- check.mclogit.drop.coefs(Y,sets,weights,X,
                                             offset = offset)
      if(any(drop.coefs)){
        warning("removing ",paste(colnames(X)[drop.coefs],collapse=",")," from model")
        X <- X[,!drop.coefs,drop=FALSE]
      }
    }
    if(ncol(X)<1)
        stop("No predictor variable remains in model")
    
    if(!length(random)){
        fit <- mclogit.fit(y=Y,s=sets,w=weights,X=X,
                       dispersion=dispersion,
                       control=control,
                       start = start,
                       offset = offset)
        groups <- NULL
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

        nn <- length(Z)
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
            Z_k <- Z[[k]]
            gconst <- groupConstInSets(groups,sets) # Is grouping factor constant within choice sets?
            if(any(gconst)){
                # If grouping factor is constant within choice sets, remove covariates that
                # are constants within choice sets
                rconst <- matConstInSets(Z_k,sets)
                if(any(rconst)){
                    cat("\n")
                    warning("removing ",
                            gsub("(Intercept)","intercept",paste(colnames(Z_k)[rconst],collapse=","),fixed=TRUE),
                            " from random part of the model\n because of insufficient within-choice set variance")
                    Z_k <- Z_k[,!rconst,drop=FALSE]
                }
                if(ncol(Z_k)<1)
                    stop("No predictor variable remains in random part of the model.\nPlease reconsider your model specification.")
            }
            d <- ncol(Z_k)
            colnames(Z_k) <- gsub("(Intercept)","(Const.)",colnames(Z_k),fixed=TRUE)
            VarCov.names.k <- rep(list(colnames(Z_k)),nlev)
            Z_k <- lapply(groups,mkZ,rX=Z_k)
            d <- rep(d,nlev)
            names(groups) <- group.labels
            list(Z_k,groups,d,VarCov.names.k)
        })
        Z <- lapply(randstruct,`[[`,1)
        groups <- lapply(randstruct,`[[`,2)
        d <- lapply(randstruct,`[[`,3)
        VarCov.names <- lapply(randstruct,`[[`,4)
        Z <- unlist(Z,recursive=FALSE)
        groups <- unlist(groups,recursive=FALSE)
        VarCov.names <- unlist(VarCov.names,recursive=FALSE)
        d <- unlist(d)
        Z <- blockMatrix(Z,ncol=length(Z))
        fit <- mmclogit.fitPQLMQL(Y,sets,weights,X,Z,
                                  d=d,
                                  method = method,
                                  estimator=estimator,
                                  control=control,
                                  offset = offset)
        nlev <- length(fit$VarCov)
        for(k in 1:nlev)
            dimnames(fit$VarCov[[k]]) <- list(VarCov.names[[k]],VarCov.names[[k]])
        names(fit$VarCov) <- names(groups)
    }
    
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
                      N=N))
    if(length(random))
        class(fit) <- c("mmclogit","mclogit","lm")
    else
        class(fit) <- c("mclogit","lm")
    fit
}


check.mclogit.drop.coefs <- function(y,
                                     s,
                                     w,
                                     X,
                                     offset){
  nvar <- ncol(X)
  nobs <- length(y)
  if(!length(offset))
    offset <- rep.int(0, nobs)
  eta <- mclogitLinkInv(y,s,w)
  pi <- mclogitP(eta,s)
  y.star <- eta - offset + (y-pi)/pi
  yP.star <- y.star - rowsum(pi*y.star,s)[s]
  XP <- X - rowsum(pi*X,s)[s,,drop=FALSE]
  ww <- w*pi
  good <- ww > 0
  wlsFit <- lm.wfit(x=XP[good,,drop=FALSE],y=yP.star[good],w=ww[good])  
  is.na(wlsFit$coef)
}

setupRandomFormula <- function(formula){

    trms <- terms(formula)
    fo <- delete.response(trms)
    
    attributes(fo) <- NULL
    if(length(fo[[2]]) < 2 || as.character(fo[[2]][1])!="|")
        stop("missing '|' operator")

    groups <- fo
    fo[2] <- fo[[2]][2]
    groups[2] <- groups[[2]][3]
    checkRandomFormula(groups[[2]])
    list(
        formula=structure(fo,class="formula"),
        groups=all.vars(groups)
    )
}

checkRandomFormula <- function(x){
    l <- as.list(x)
    if(length(l) < 3) return(NULL)
    if(!as.character(l[[1]])=="/") stop("Invalid random formula",call.=FALSE)
    x <- x[[2]]
    if(length(x)>1) Recall(x)
}



print.mclogit <- function(x,digits= max(3, getOption("digits") - 3), ...){
    cat(paste(deparse(x$call), sep="\n", collapse="\n"), "\n\n", sep="")
    if(length(coef(x))) {
        cat("Coefficients")
        if(is.character(co <- x$contrasts))
            cat("  [contrasts: ",
                apply(cbind(names(co),co), 1, paste, collapse="="), "]")
        cat(":\n")
        print.default(format(x$coefficients, digits=digits),
                      print.gap = 2, quote = FALSE)
    } else cat("No coefficients\n\n")
    if(x$phi != 1)
        cat("\nDispersion: ",x$phi)
    cat("\nNull Deviance:    ",   format(signif(x$null.deviance, digits)),
        "\nResidual Deviance:", format(signif(x$deviance, digits)))
    if(!x$converged) cat("\nNote: Algorithm did not converge.\n")
    if(nchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep="")
    else cat("\n")
    invisible(x)
}

vcov.mclogit <- function(object,...){
    phi <- object$phi
    if(!length(phi)) phi <- 1
    return(object$covmat * phi)
}

weights.mclogit <- function(object,...){
  return(object$weights)
}

deviance.mclogit <- function(object,...){
  return(object$deviance)
}

summary.mclogit <- function(object,dispersion=NULL,correlation = FALSE, symbolic.cor = FALSE,...){

    ## calculate coef table

    coef <- object$coefficients

    if(is.null(dispersion))
        dispersion <- object$phi
    
    covmat.scaled <- object$covmat * dispersion
    
    var.cf <- diag(covmat.scaled)
    s.err <- sqrt(var.cf)
    zvalue <- coef/s.err

    if(dispersion == 1)
        pvalue <- 2*pnorm(-abs(zvalue))
    else
        pvalue <- 2*pt(-abs(zvalue),df=object$df.residual)

    coef.table <- array(NA,dim=c(length(coef),4))
    rownames(coef.table) <- names(coef)
    if(dispersion == 1)
        colnames(coef.table) <- c("Estimate", "Std. Error","z value","Pr(>|z|)")
    else
        colnames(coef.table) <- c("Estimate", "Std. Error","t value","Pr(>|t|)")
    coef.table[,1] <- coef
    coef.table[,2] <- s.err
    coef.table[,3] <- zvalue
    coef.table[,4] <- pvalue

    ans <- c(object[c("call","terms","deviance","contrasts",
                      "null.deviance","iter","na.action","model.df",
                      "df.residual","N","converged")],
              list(coefficients = coef.table,
                   cov.coef=covmat.scaled,
                   dispersion = dispersion
                   ))
    p <- length(coef)
    if(correlation && p > 0) {
        dd <- sqrt(diag(ans$cov.coef))
        ans$correlation <-
            ans$cov.coef/outer(dd,dd)
        ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.mclogit"
    return(ans)
}

print.summary.mclogit <-
    function (x, digits = max(3, getOption("digits") - 3),
              symbolic.cor = x$symbolic.cor,
              signif.stars = getOption("show.signif.stars"), ...){
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep="\n", collapse="\n"), "\n\n", sep="")

    coefs <- x$coefficients
    printCoefmat(coefs, digits=digits, signif.stars=signif.stars,
                     na.print="NA", ...)
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
    if(nchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n\n", sep="")
    else cat("\n\n")
    invisible(x)
}




fitted.mclogit <- function(object,type=c("probabilities","counts"),...){
  weights <- object$weights
  
  res <- object$fitted.values
  type <- match.arg(type)
  
  na.act <- object$na.action
  
  res <- switch(type,
          probabilities=res,
          counts=weights*res)
          
  if(is.null(na.act))
    res
  else 
    napredict(na.act,res)
}

predict.mclogit <- function(object, newdata=NULL,type=c("link","response"),se.fit=FALSE,...){

    type <- match.arg(type)
    fo <- object$formula
    if(as.character(fo[[2]][[1]])=="|")
        fo[[2]][[1]] <- as.name("cbind")
    lhs <- fo[[2]]
    rhs <- fo[-2]
    if(length(lhs)==3)
        sets <- lhs[[3]]
    else stop("no way to determine choice set ids")
    if(missing(newdata)){
        m <- model.frame(fo,data=object$data)
        set <- m[[1]][,2]
        na.act <- object$na.action
    }
    else{
        lhs <- lhs[[3]]
        fo[[2]] <- lhs
        m <- model.frame(fo,data=newdata)
        set <- m[[1]]
        na.act <- attr(m,"na.action")
    }
    X <- model.matrix(rhs,m,
                      contasts.arg=object$contrasts,
                      xlev=object$xlevels
                      )

    cf <- coef(object)
    X <- X[,names(cf), drop=FALSE]
    
    eta <- c(X %*% cf)
    if(se.fit){
        V <- vcov(object)
        stopifnot(ncol(X)==ncol(V))
    }

    
    if(type=="response") {
        
        set <- match(set,unique(set))
        exp.eta <- exp(eta)
        sum.exp.eta <- rowsum(exp.eta,set)
        p <- exp.eta/sum.exp.eta[set]
        if(se.fit){
            wX <- p*(X - rowsum(p*X,set)[set,,drop=FALSE])
            se.p <- sqrt(rowSums(wX * (wX %*% V)))
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
        se.eta <- sqrt(rowSums(X * (X %*% V)))
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

logLik.mclogit <- function(object,...){
    if (length(list(...)))
        warning("extra arguments discarded")
    val <- if(length(object$ll))
            object$ll
           else NA
    attr(val, "nobs") <- object$N
    attr(val, "df") <- object$model.df
    class(val) <- "logLik"
    return(val)
}

residuals.mclogit <- 
 function(object, type = c("deviance", "pearson", "working",
                           "response", "partial"), ...){
   type <- match.arg(type)
   
   resid <- switch(type,
                 deviance=mclogit.dev.resids(object),
                 pearson=stop("not yet implemented"),
                 working=object$working.residuals,
                 response=object$response.residuals,
                 partial=stop("not yet implemented")
                 )
   naresid(object$na.action,resid)
 }
 
mclogit.dev.resids <- function(obj){
 y <- obj$y
 s <- obj$s
 w <- obj$weights
 pi <- obj$fitted.values
 
 n <- w*y+0.5
 f <- n/(rowsum(n,s)[s])
 #sign(y-p)*sqrt(2*abs(log(f)-log(y)))
 r <- 2*(f*log(f/pi))
 r - ave(r,s)
}


nobs.mclogit <- function(object,...) object$N

extractAIC.mclogit <- function(fit, scale = 0, k = 2, ...) 
{
  N <- fit$N
  edf <- N - fit$df.residual
  aic <- AIC(fit)
  c(edf, aic + (k - 2) * edf)
}

weights.mclogit <- function(object, type = c("prior", "working"),...) {
  type <- match.arg(type)
  res <- if (type == "prior") 
    object$prior.weights
  else object$weights
  if (is.null(object$na.action)) 
    res
  else naresid(object$na.action, res)
}

print.mmclogit <- function(x,digits= max(3, getOption("digits") - 3), ...){
    cat(paste(deparse(x$call), sep="\n", collapse="\n"), "\n\n", sep="")
    if(length(coef(x))) {
        cat("Coefficients")
        if(is.character(co <- x$contrasts))
            cat("  [contrasts: ",
                apply(cbind(names(co),co), 1, paste, collapse="="), "]")
        cat(":\n")
        print.default(format(x$coefficients, digits=digits),
                      print.gap = 2, quote = FALSE)
    } else cat("No coefficients\n\n")

    cat("\n(Co-)Variances:\n")
    VarCov <- x$VarCov
    names(VarCov) <- names(x$groups)
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

vcov.mmclogit <- function(object,...){
    info.coef <- object$info.coef
    vcov.cf <- solve(info.coef)
    return(vcov.cf)
}



summary.mmclogit <- function(object,dispersion=NULL,correlation = FALSE, symbolic.cor = FALSE,...){

    ## calculate coef table

    coef <- object$coefficients
    info.coef <- object$info.coef
    vcov.cf <- solve2(info.coef)
    var.cf <- diag(vcov.cf)
    s.err <- sqrt(var.cf)
    zvalue <- coef/s.err
    pvalue <- 2*pnorm(-abs(zvalue))

    coef.table <- array(NA,dim=c(length(coef),4))
    dimnames(coef.table) <- list(names(coef),
            c("Estimate", "Std. Error","z value","Pr(>|z|)"))
    coef.table[,1] <- coef
    coef.table[,2] <- s.err
    coef.table[,3] <- zvalue
    coef.table[,4] <- pvalue

    VarCov <- object$VarCov
    info.lambda <- object$info.lambda
    se_VarCov <- se_Phi(VarCov,info.lambda)

    names(VarCov) <- names(object$groups)
    names(se_VarCov) <- names(VarCov)
    
    ans <- c(object[c("call","terms","deviance","contrasts",
                      "null.deviance","iter","na.action","model.df",
                      "df.residual","groups","N","converged")],
              list(coefficients = coef.table,
                   vcov.coef = vcov.cf,
                   VarCov    = VarCov,
                   se_VarCov = se_VarCov))
    p <- length(coef)
    if(correlation && p > 0) {
        dd <- sqrt(diag(ans$cov.coef))
        ans$correlation <-
            ans$cov.coef/outer(dd,dd)
        ans$symbolic.cor <- symbolic.cor
    }

    ans$ngrps <- sapply(object$groups,nlevels)
    
    class(ans) <- "summary.mmclogit"
    return(ans)
}


print.summary.mmclogit <-
    function (x, digits = max(3, getOption("digits") - 3),
              symbolic.cor = x$symbolic.cor,
              signif.stars = getOption("show.signif.stars"), ...){
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep="\n", collapse="\n"), "\n\n", sep="")

    coefs <- x$coefficients
    cat("Coefficents:\n")
    printCoefmat(coefs, digits=digits, signif.stars=signif.stars,
                     na.print="NA", ...)

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
    if(nchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n\n", sep="")
    else cat("\n\n")
    invisible(x)
}

predict.mmclogit <- function(object, newdata=NULL,type=c("link","response"),se.fit=FALSE,
                             conditional=TRUE, ...){
    
    type <- match.arg(type)
    fo <- object$formula
    if(as.character(fo[[2]][[1]])=="|")
        fo[[2]][[1]] <- as.name("cbind")
    lhs <- fo[[2]]
    rhs <- fo[-2]
    random <- object$random  
    if(length(lhs)==3)
        sets <- lhs[[3]]
    else stop("no way to determine choice set ids")
    if(missing(newdata)){
        mf <- object$model
        sets <- mf[[1]][,2]
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
        sets <- eval(sets,newdata)
        na.act <- attr(mf,"na.action")
    }
    X <- model.matrix(rhs,mf,
                      contrasts.arg=object$contrasts,
                      xlev=object$xlevels
                      )
    cf <- coef(object)
    X <- X[,names(cf), drop=FALSE]
    eta <- c(X %*% cf)

    if(object$method=="PQL" && conditional){
        
        rf <- lapply(random,"[[","formula")
        rt <- lapply(rf,terms)
        suppressWarnings(Z <- lapply(rt,model.matrix,rmf,
                                     contrasts.arg=object$contrasts,
                                     xlev=object$xlevels))
        d <- sapply(Z,ncol)
        nn <- length(Z)

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
            olevels <- olevels[group.labels]
            groups <- Map(factor,x=groups,levels=olevels)
            
            VarCov.names.k <- rep(list(colnames(Z[[k]])),nlev)
            Z_k <- lapply(groups,mkZ,rX=Z[[k]])
            d <- rep(d[k],nlev)
            names(groups) <- group.labels
            list(Z_k,groups,d,VarCov.names.k)
        })
        Z <- lapply(randstruct,`[[`,1)
        groups <- lapply(randstruct,`[[`,2)
        Z <- unlist(Z,recursive=FALSE)
        d <- lapply(randstruct,`[[`,3)
        groups <- unlist(groups,recursive=FALSE)
        d <- unlist(d)
        
        Z <- blockMatrix(Z)
        b <- object$random.effects
        nlev <- length(Z)

        for(k in 1:nlev)
            eta <- eta +  as.vector(Z[[k]]%*%b[[k]])
    }
    
    nvar <- ncol(X)
    nobs <- nrow(X)
    
    if(type=="response" || object$method=="PQL" && conditional ){
        j <- match(sets,unique(sets))
        exp.eta <- exp(eta)
        sum.exp.eta <- rowsum(exp.eta,j)
        p <- exp.eta/sum.exp.eta[j]
    }
    if(se.fit && (type=="response" || object$method=="PQL" && conditional)){
        nsets <- j[length(j)]
        W <- Matrix(0,nrow=nobs,ncol=nsets)
        i <- 1:nobs
        W[cbind(i,j)] <- p
        W <- Diagonal(x=p)-tcrossprod(W)
        WX <- W%*%X
        if(object$method=="PQL" && conditional){
            WZ <- bMatProd(W,Z)
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
                se.p <- sqrt(var.p)
            }
            else {
                vcov.coef <- vcov(object)
                se.p <- sqrt(rowSums(WX*(WX%*%vcov.coef)))
            }
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
        if(se.fit){
            if(object$method=="PQL" && conditional){
                XZ <- structure(cbind(blockMatrix(X),Z),class="blockMatrix")
                var.eta <- bMatProd(XZ,K)
                var.eta <- Map(`*`,XZ,var.eta)
                var.eta <- lapply(var.eta,rowSums)
                var.eta <- Reduce(`+`,var.eta)
            }
            else {
                vcov.coef <- vcov(object)
                var.eta <- rowSums(X*(X%*%vcov.coef))
            }
            se.eta <- sqrt(var.eta)
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



tr <- function(x) sum(diag(x))

mkZ <- function(groups,rX){

    n <- length(groups)
    m <- nlevels(groups)
    p <- ncol(rX)
    
    Z <- Matrix(0,nrow=n,ncol=m*p)

    i <- 1:n
    k <- 1:p
    j <- as.integer(groups)
    
    i <- rep(i,p)
    jk <- rep((j-1)*p,p)+rep(k,each=n)
    i.jk <- cbind(i,jk)

    Z[i.jk] <- rX
    Z
}

mkZ2 <- function(all.groups,
                groups,
                rX){
    n <- length(groups)
    ug <- unique(all.groups)
    m <- length(ug)
    p <- ncol(rX)
    
    Z <- Matrix(0,nrow=n,ncol=m*p)

    i <- 1:n
    k <- 1:p
    j <- groups
    
    i <- rep(i,p)
    jk <- rep((j-1)*p,p)+rep(k,each=n)
    i.jk <- cbind(i,jk)

    Z[i.jk] <- rX
    Z
}



mkG <- function(rX){

    p <- ncol(rX)
    nms <- colnames(rX)
    
    G <- matrix(0,p,p)
    ltT <- lower.tri(G,diag=TRUE)
    ltF <- lower.tri(G,diag=FALSE)

    n <- p*(p+1)/2
    m <- p*(p-1)/2
    
    diag(G) <- 1:p
    G[ltF] <- p + 1:m
    G <- lwr2sym(G)
    rownames(G) <- colnames(G) <- nms
    
    lapply(1:n,mkG1,G)
}

mkG1 <- function(i,G) Matrix(array(as.integer(i==G),
                                   dim=dim(G),
                                   dimnames=dimnames(G)
                                   ))

fillG <- function(G,theta){

    Phi <- Map(`*`,theta,G)

    if(length(Phi)>1){
        for(i in 2:length(Phi))
            Phi[[1]] <- Phi[[1]] + Phi[[i]]
    }
    Phi[[1]]
}


lunq <- function(x)length(attr(x,"unique"))

G.star1 <- function(I,G)Map(`%x%`,list(I),G)
quadform <- function(A,x) as.numeric(crossprod(x,A%*%x))
tr.crossprod <- function(A,B) sum(A*B)

lwr2sym <- function(X){

    lwrX <- lower.tri(X)
    x.lwr <- X[lwrX]
    Y <- t(X)
    Y[lwrX] <- x.lwr
    Y
}


fuseMat <- function(x){
    if(ncol(x)>1){
        y <- lapply(1:nrow(x),
                    fuseCols,x=x)
    }
    else
        y <- x
    
    do.call(rbind,y)
}

cbindList <- function(x) do.call(cbind,x)
fuseCols <- function(x,i) do.call(cbind,x[i,]) 

format_Mat <- function(x,title="",rownames=NULL){
    if(length(rownames))
        rn <- format(c("",rownames))
    else 
        rn <- format(c("",rownames(x)))
    x <- format(x)
    x <- apply(x,1,paste,collapse=" ")
    x <- format(c(title,x))
    paste(rn,x)
}

update.mclogit <-  function(object, formula., dispersion, ...) {
    if(!inherits(object,"mmclogit") &&
       (missing(formula.) || formula. == object$formula)
       && !missing(dispersion))
        update_mclogit_dispersion(object,dispersion)
    else NextMethod()
}


getFirst <- function(x) x[1]

simulate.mclogit <- function(object, nsim = 1, seed = NULL, ...){

    if(object$phi > 1)
        stop("Simulating responses from models with oversdispersion is not supported yet")
    
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
    
    weights <- object$weights
    probs <- object$fitted.values
    set <- object$s
    i <- 1:length(probs)
    
    probs <- split(probs,set)
    weights <- split(weights,set)
    i <- split(i,set)
    weights <- sapply(weights,getFirst)

    yy <- mapply(rmultinom,size=weights,prob=probs,
                 MoreArgs=list(n=nsim),SIMPLIFY=FALSE)
    yy <- do.call(rbind,yy)

    i <- unlist(i)
    yy[i,] <- yy
    rownames(yy) <- names(object$working.residuals)
    colnames(yy) <- paste0("sim_",1:nsim)
    yy <- as.data.frame(yy)
    attr(yy,"seed") <- RNGstate
    yy
}

simulate.mmclogit <- function(object, nsim = 1, seed = NULL, ...)
    stop("Simulating responses from random-effects models is not supported yet")

eigen.solve <- function(x){
    ev <- eigen(x)
    d <- ev$values
    V <- ev$vectors
    id <- 1/d
    V %*% (id*t(V))
}

solve2 <- function(x){
    ix <- try(solve(x))
    if(inherits(ix,"try-error"))
        return(eigen.solve(x))
    else return(ix)
}
