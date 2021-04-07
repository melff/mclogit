### emmeans support for mblogit and mmblogit models

recover_data.mblogit <- function(object, ...) {
    fcall <- object$call
    emmeans::recover_data(fcall, delete.response(terms(object)), object$na.action, ...)
}

emm_basis.mblogit <- function(object, trms, xlev, grid, 
                              mode = c("prob", "latent"), ...) {
    mode <- match.arg(mode)
    bhat <- t(object$coefmat)
    V <- emmeans::.my.vcov(object, ...)
    # we have to rearrange the vcov elements in row-major order...
    perm <- matrix(seq_along(as.numeric(object$coefmat)), 
                   ncol <- ncol(object$coefmat))
    perm <- as.numeric(t(perm))
    V <- V[perm, perm]
    k <- ifelse(is.matrix(object$coefmat), ncol(bhat), 1)
    m <- model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X <- model.matrix(trms, m, contrasts.arg = object$contrasts)
    # recenter for latent predictions
    pat <- (rbind(0, diag(k + 1, k)) - 1) / (k + 1)
    X <- kronecker(pat, X)
    nbasis <- estimability::all.estble
    nbasis <- kronecker(rep(1,k), nbasis)
    misc <- list(tran = "log", inv.lbl = "e^y")
    dfargs <- list(df = Inf)
    dffun <- function(k, dfargs) dfargs$df
    ylevs <- list(class = levels(object$model[[1]]))
    if (is.null(ylevs)) 
        ylevs <- list(class = seq_len(k))
    names(ylevs) <- as.character.default(eval(object$call$formula, environment(trms))[[2]])
    misc$ylevs <- ylevs
    if (mode == "prob")
        misc$postGridHook <- .multinom.postGrid
    list(X = X, bhat = as.numeric(bhat), nbasis = nbasis, V = V, 
         dffun = dffun, dfargs = dfargs, misc = misc)
}


### post-processing when mode = "prob" (copied from emmeans support for "multinom" objects)
.multinom.postGrid <- function(object, N.sim, ...) {
    linfct <- object@linfct
    misc <- object@misc
    # grid will have multresp as slowest-varying factor...
    idx <- matrix(seq_along(linfct[, 1]), 
                 ncol = length(object@levels[[object@roles$multresp]]))
    bhat <- as.numeric(idx) # right length, contents will be replaced
    if(sim <- !missing(N.sim)) {
        message("Simulating a sample of size ", N.sim)
        bsamp <- mvtnorm::rmvnorm(N.sim, object@bhat, object@V)
        postb <- matrix(0, nrow = N.sim, ncol = length(bhat))
    }
    for (i in 1:nrow(idx)) {
        rows <- idx[i, ]
        exp.psi <- exp(linfct[rows, , drop = FALSE] %*% object@bhat)
        p <- as.numeric(exp.psi / sum(exp.psi))
        bhat[rows] <- p
        if (sim) {
            ex <- exp(linfct[rows, , drop = FALSE] %*% t(bsamp))  # p x N
            px <- t(apply(ex, 2, function(x) x / sum(x)))
            postb[, rows] <- px
        }
        A <- emmeans::.diag(p) - outer(p, p)    # partial derivs
        linfct[rows, ] <- A %*% linfct[rows, ]
    }
    misc$postGridHook <- misc$tran <- misc$inv.lbl <- NULL
    misc$estName <- "prob"
    
    object@bhat <- bhat
    object@V <- linfct %*% tcrossprod(object@V, linfct)
    object@linfct <- diag(1, length(bhat))
    object@misc <- misc
    if (sim)
        object@post.beta <- postb
    object
}

### Documentation & example: See emmeans-support.Rd

