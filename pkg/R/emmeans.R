## Added by Russel V. Lenth
### emmeans support for mblogit and mmblogit models

recover_data.mblogit <- function(object, ...) {
    rd <- get("recover_data.multinom", asNamespace("emmeans"))
    rd(object, ...)
}

emm_basis.mblogit <- function(object, trms, xlev, grid, 
                              mode = c("prob", "latent"), vcov., ...) {
    object$coefficients <- object$coefmat
    object$lev <- levels(object$model[[1]])
    object$edf <- Inf
    # we have to rearrange the vcov elements in row-major order
    if(missing(vcov.))
        vcov. <- vcov(object)
    perm <- matrix(seq_along(as.numeric(object$coefmat)), 
                  ncol = ncol(object$coefmat))
    perm <- as.numeric(t(perm))
    vcov. <- vcov.[perm, perm]
    emb <- get("emm_basis.multinom", asNamespace("emmeans"))
    emb(object, trms = trms, xlev = xlev, grid = grid, mode = mode, vcov. = vcov., ...)
}


