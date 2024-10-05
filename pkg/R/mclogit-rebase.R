
get_categs <- function(object){
    D <- object$D
    rownames(D)
}

get_baseline_cat <- function(object){
    D <- object$D
    j <- which(!rownames(D)%in%colnames(D))
    rownames(D)[j]
}

rebase_mat <- function(categs,from,to){
    m <- length(categs)
    j <- match(from,categs)
    k <- match(to,categs)
    res <- diag(nrow=m)
    dimnames(res) <- list(categs,categs)
    res[,k] <- -1
    res <- res[,-j] 
    res <- res[-k,]
    res
}


#' Change baseline category of multinomial logit or similar model
#'
#' `rebase` returns an model object that is equivalent to the one
#' given as argument but differs in parameterization
#'
#' @param object a statistical model object
#' @param to usually, a string; the baseline category
#' @param ... other arguments, currently ignored
rebase <- function(object,to,...) UseMethod("rebase")

#' @rdname rebase
rebase.mblogit <- function(object,to,...){
    categs <- get_categs(object)
    m <- length(categs)
    from <- get_baseline_cat(object)
    TMat <- rebase_mat(categs,from=from,to=to)
    coefmat <- object$coefmat
    p <- ncol(coefmat)
    coefmat.new <- TMat%*%coefmat
    coefficients.new <- as.vector(coefmat.new)
    coefficients.new.names <- outer(rownames(coefmat.new),colnames(coefmat.new),paste,sep="~")
    coefficients.new.names <- as.vector(coefficients.new.names)
    names(coefficients.new) <- coefficients.new.names
    iTMat <- rebase_mat(categs,from=to,to=from)
    iMMat <- diag(p)%x%t(iTMat)
    info.matrix <- object$information.matrix
    info.matrix.new <- iMMat%*%info.matrix%*%t(iMMat)
    dimnames(info.matrix.new) <- list(coefficients.new.names,
                                      coefficients.new.names)
    D.new <- diag(m)
    dimnames(D.new) <- list(categs,categs)
    D.new <- D.new[,-match(to,categs)]
    object.new <- object
    object.new$coefmat <- coefmat.new
    object.new$coefficients <- coefficients.new
    object.new$information.matrix <- info.matrix.new
    object.new$D <- D.new
    object.new
}
