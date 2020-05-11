blockMatrix <- function(x=list(),nrow=1,ncol=1){
    if(!is.list(x)) x <- list(x)
    y <- matrix(x,nrow=nrow,ncol=ncol)
    structure(y,class="blockMatrix")
}

Ops.blockMatrix <- function(e1, e2){
    if(!inherits(e1,"blockMatrix")) e1 <- blockMatrix(e1)
    if(!inherits(e2,"blockMatrix")) e2 <- blockMatrix(e2)
    stopifnot(dim(e1)==dim(e2))
    d <- dim(e1)
    if(!(.Generic%in% c("+","-","*","==")))
        stop(sQuote(.Generic)," not implemented for block matrices")
    res <- switch(.Generic,
                  `+`= mapply(`+`,e1,e2,SIMPLIFY=FALSE),
                  `-`= mapply(`-`,e1,e2,SIMPLIFY=FALSE),
                  `*`= mapply(`*`,e1,e2,SIMPLIFY=FALSE),
                  `==`= all(Reduce(`&`,mapply(`==`,e1,e2)))
                  )
    if(is.list(res)){
        dim(res) <- d
        structure(res,                  
                  class=class(e1))
    }
    else res
}

bMatProd <- function(x,y){
    if(!inherits(x,"blockMatrix")) x <- blockMatrix(x)
    if(!inherits(y,"blockMatrix")) y <- blockMatrix(y)
    dim.x <- dim(x)
    dim.y <- dim(y)
    stopifnot(dim.x[2]==dim.y[1])
    m <- dim.x[1]
    n <- dim.y[2]
    q <- dim.x[2]
    res <- blockMatrix(nrow=m,ncol=n)
    for(i in 1:m)
        for(j in 1:n){
            res[[i,j]] <- inner_p(x[i,],y[,j])
        }
    res
}

bMatCrsProd <- function(x,y=NULL){
    if(missing(y))
        y <- x
    if(!inherits(x,"blockMatrix")) x <- blockMatrix(x)
    if(!inherits(y,"blockMatrix")) y <- blockMatrix(y)
    dim.x <- dim(x)
    dim.y <- dim(y)
    stopifnot(dim.x[1]==dim.y[1])
    m <- dim.x[2]
    n <- dim.y[2]
    q <- dim.x[1]
    res <- blockMatrix(nrow=m,ncol=n)
    for(i in 1:m)
        for(j in 1:n){
            res[[i,j]] <- inner_crsp(x[,i],y[,j])
        }
    res
}

bMatTCrsProd <- function(x,y=NULL){
    if(missing(y))
        y <- x
    if(!inherits(x,"blockMatrix")) x <- blockMatrix(x)
    if(!inherits(y,"blockMatrix")) y <- blockMatrix(y)
    dim.x <- dim(x)
    dim.y <- dim(y)
    stopifnot(dim.x[2]==dim.y[2])
    m <- dim.x[1]
    n <- dim.y[1]
    q <- dim.x[2]
    res <- blockMatrix(nrow=m,ncol=n)
    for(i in 1:m)
        for(j in 1:n){
            res[[i,j]] <- inner_tcrsp(x[i,],y[j,])
        }
    res
}

bMatTrns <- function(x){
    m <- nrow(x)
    n <- ncol(x)
    res <- blockMatrix(nrow=n,ncol=m)
    for(i in 1:m)
        for(j in 1:n){
            res[[i,j]] <- t(res[[i,j]])
        }
    res
}

inner_p <- function(x,y){
    xy <- mapply(`%*%`,x,y,SIMPLIFY=FALSE)
    Reduce(`+`,xy)
}

inner_crsp <- function(x,y){
    xy <- mapply(crossprod,x,y,SIMPLIFY=FALSE)
    Reduce(`+`,xy)
}

inner_tcrsp <- function(x,y){
    xy <- mapply(tcrossprod,x,y,SIMPLIFY=FALSE)
    Reduce(`+`,xy)
}



matprod1 <- function(x,y){
    if(!length(x) || !length(y)) NULL
    else x %*% y
}

blockDiag <- function(x,n=length(x)){
    y <- blockMatrix(nrow=n,ncol=n)
    i <- 1:n
    y[cbind(i,i)] <- x
    bM_fill(y)
}

bM_check <- function(x){
    nnrow <- sapply(x,NROW)
    nncol <- sapply(x,NCOL)
    dim(nnrow) <- dim(x)
    dim(nncol) <- dim(x)
    lunq.cols <- apply(nncol,2,lunq)
    lunq.rows <- apply(nnrow,1,lunq)
    ok <- all(lunq.cols==1) && all(lunq.cols)
    return(ok)
}

bM_nrow <- function(x) sapply(x[,1],nrow)

bM_ncol <- function(x) sapply(x[1,],ncol)

to_bM <- function(x,nnrow,nncol){
    nnrow1 <- cumsum(c(0,nnrow[-length(nnrow)])) + 1
    nncol1 <- cumsum(c(0,nncol[-length(nncol)])) + 1
    rows <- mapply(seq.int,from=nnrow1,length.out=nnrow,SIMPLIFY=FALSE)
    cols <- mapply(seq.int,from=nncol1,length.out=nncol,SIMPLIFY=FALSE)
    m <- length(nnrow)
    n <- length(nncol)
    y <- blockMatrix(nrow=m,ncol=n)
    for(i in 1:m)
        for(j in 1:n)
            y[i,j] <- list(Matrix(x[rows[[i]],cols[[j]]]))
    return(y)
}

bM_fill <- function(x){
    nnrow <- Sapply(x,NROW)
    nncol <- Sapply(x,NCOL)
    dim(nnrow) <- dim(x)
    dim(nncol) <- dim(x)
    nnrow <- apply(nnrow,1,max)
    nncol <- apply(nncol,2,max)
    m <- nrow(x)
    n <- ncol(x)
    for(i in 1:m)
        for(j in 1:n){
            if(is.null(x[[i,j]])){
                x[[i,j]] <- Matrix(0,nnrow[i],nncol[j])
            }
        }
    return(x)
}

solve.blockMatrix <- function(a,b,...){
    nnrow.a <- bM_nrow(a)
    nncol.a <- bM_ncol(a)
    a <- fuseMat(a)
    if(missing(b)){
        x <- solve(a)
        return(to_bM(x,nnrow=nnrow.a,nncol=nncol.a))
    }
    else {
        nnrow.b <- bM_nrow(b)
        nncol.b <- bM_ncol(b)
        b <- fuseMat(b)
        x <- solve(a,b)
        return(to_bM(x,nnrow=nnrow.a,nncol=nncol.b))
    }
}

format_dims <- function(x){
    sprintf("<%d x %d>",nrow(x),ncol(x))
}

print.blockMatrix <- function(x,quote=FALSE,...){
    cat(sprintf("Block matrix with %d x %d blocks\n\n",nrow(x),ncol(x)))
    y <- sapply(x,format_dims)
    dim(y) <- dim(x)
    print.default(y,quote=quote,...)
    invisible(x)
}

sum_blockDiag <- function(x,n){
    i <- rep(1:n,n)
    j <- rep(1:n,each=n)
    nblks <- nrow(x) %/% n
    offs <- rep(seq.int(from=0,to=nblks-1),each=n*n)
    i <- rep(i,nblks) + offs
    j <- rep(j,nblks) + offs
    y <- x[cbind(i,j)]
    dim(y) <- c(n*n,nblks)
    y <- rowSums(y)
    dim(y) <- c(n,n)
    Matrix(y)
}

v_bCrossprod <- function(x,d){
    n <- length(x)%/%d
    dim(x) <- c(d,n)
    tcrossprod(x)
}

v_bQuadfm <- function(x,W){
    d <- nrow(W)
    n <- length(x)%/%d
    dim(x) <- c(d,n)
    colSums((W%*%x)*x)
}

set_blockDiag <- function(x,v){
    n <- ncol(v)
    i <- rep(1:n,n)
    j <- rep(1:n,each=n)
    nblks <- ncol(x) %/% n
    offs <- rep(seq.int(from=0,to=nblks-1)*n,each=n*n)
    i <- rep(i,nblks) + offs
    j <- rep(j,nblks) + offs
    x[cbind(i,j)] <- v
    return(x)
}

logDet_blockMatrix <- function(x){
    d <- determinant(fuseMat(x),logarithm=TRUE)
    d$modulus
}

chol_blockMatrix <- function(x,resplit=TRUE){
    y <- chol(fuseMat(x))
    if(resplit){
        nnrow <- bM_nrow(x)
        nncol <- bM_ncol(x)
        return(to_bM(y,nnrow=nnrow,nncol=nncol))
    }
    else return(y)
}

kron_bM <- function(x,y){
    m1 <- nrow(x)
    m2 <- nrow(y)
    n1 <- ncol(x)
    n2 <- ncol(y)
    attributes(x) <- NULL
    attributes(y) <- NULL
    lx <- length(x)
    ly <- length(y)
    x <- rep(x,each=ly)
    y <- rep(y,lx)
    xy <- mapply(`%x%`,x,y,SIMPLIFY=FALSE)
    blockMatrix(xy,m1*m2,n1*n2)
}
