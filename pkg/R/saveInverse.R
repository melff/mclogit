safeInverse <- function(x,tol=1e-7){
    tryCatch(solve(x),
             error=function(e){
                 warning(e$message,call.=FALSE,immediate.=TRUE)
                 warning("saveInverse: Using Moore-Penrose inverse",call.=FALSE,immediate.=TRUE)
                 moore.penrose(x,tol=tol)
             })
}

mach.eps <- .Machine$double.eps

moore.penrose <- function(x,tol=mach.eps*max(dim(x))*max(abs(d))){
    svd.x <- svd(x)
    d <- svd.x$d
    u <- svd.x$u
    v <- svd.x$v
    good <- abs(d) > tol
    id <- 1/d
    id[!good] <- 0
    v %*% diag(id,nrow=length(id)) %*% t(u)
}
