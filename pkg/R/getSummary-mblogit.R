rbind_list <- function(x) do.call(rbind, x)

getSummary.mblogit <- function(obj,
                               alpha = .05,
                               ...) {
  smry <- summary(obj)
  N <- obj$N
  coef <- smry$coefficients

  lower.cf <- qnorm(p = alpha / 2, mean = coef[, 1], sd = coef[, 2])
  upper.cf <- qnorm(p = 1 - alpha / 2, mean = coef[, 1], sd = coef[, 2])
  coef <- cbind(coef, lower.cf, upper.cf)
  ttl <- c("est", "se", "stat", "p", "lwr", "upr")
  colnames(coef) <- ttl

  modcat <- colnames(obj$D)
  basecat <- rownames(obj$D)[rownames(obj$D) %nin% modcat]

  eqs <- paste0(modcat, "~")

  rn.coef <- rownames(coef)
  coef.grps <- lapply(eqs, function(eq) {
    ii <- grep(eq, rn.coef, fixed = TRUE)
    coef.grp <- coef[ii, , drop = FALSE]
    rownames(coef.grp) <- gsub(eq, "", rownames(coef.grp), fixed = TRUE)
    coef.grp
  })

  if (getOption("mblogit.show.basecat", TRUE)) {
    grp.titles <- paste(modcat, basecat, sep = getOption("mblogit.basecat.sep", "/"))
  } else {
    grp.titles <- modcat
  }

  names(coef.grps) <- grp.titles
  coef <- do.call(memisc::collect, coef.grps)

  VarPar <- NULL
  VarCov <- smry$VarCov
  se_VarCov <- smry$se_VarCov

  n.eq <- length(eqs)

  for (i in seq_along(VarCov)) {
    lv.i <- names(VarCov)[i]
    vc.i <- VarCov[[i]]
    se_vc.i <- se_VarCov[[i]]
    vp.i <- array(NA, c(
      nrow(vc.i),
      ncol(vc.i),
      6
    ))
    vp.i[, , 1] <- vc.i
    vp.i[, , 2] <- se_vc.i
    m.i <- ncol(vc.i) %/% n.eq
    d <- c(n.eq, m.i)
    dim(vp.i) <- c(d, d, 6)
    vn.i <- colnames(vc.i)
    vn.i <- strsplit(vn.i, "~")
    vn.i <- unique(sapply(vn.i, "[", 2))
    dn <- list(eqs, vn.i)
    dimnames(vp.i) <- c(dn, dn, list(ttl))
    vp.i.arr <- aperm(vp.i, c(4, 2, 3, 1, 5))

    # vp.i <- lapply(eqs,function(eq){
    #   ii <- grep(eq,dn.4,fixed=TRUE)
    #   browser()
    #   vp.i.grp <- vp.i[,,,ii,,drop=FALSE]
    #   nr.i.g <- nrow(vp.i.grp)
    #   nc.i.g <- ncol(vp.i.grp)
    #   dn1.i.grp <- dimnames(vp.i.grp)[[1]]
    #   dn2.i.grp <- dimnames(vp.i.grp)[[2]]
    #   dn2.i.grp <- gsub(eq,"~",dn2.i.grp,fixed=TRUE)
    #   dn3.i.grp <- dimnames(vp.i.grp)[[3]]
    #   dim(vp.i.grp) <- c(nr.i.g*nc.i.g,6)
    #   rn.i.g.1 <- rep(dn1.i.grp,nc.i.g)
    #   rn.i.g.2 <- rep(dn2.i.grp,each=nr.i.g)
    #   #rn.i.g <- ifelse(dn1.i.grp == dn2.i.grp,"Var","Cov")
    #   rn.i.g <- paste0(rn.i.g.1,",",rn.i.g.2)
    #   rownames(vp.i.grp) <- rn.i.g
    #   colnames(vp.i.grp) <- dn3.i.grp
    #   vp.i.grp
    # })
    vp.i_ <- matrix(list(NULL), n.eq, n.eq)
    for (j in 1:n.eq) {
      for (k in 1:n.eq) {
        vp.ijk <- vp.i.arr[, , j, k, ]
        dim(vp.ijk) <- c(m.i^2, 6)
        rn.i.1 <- rep(vn.i, m.i)
        rn.i.2 <- rep(vn.i, each = m.i)
        jk.1 <- rep(1:m.i, m.i)
        jk.2 <- rep(1:m.i, each = m.i)
        rownames(vp.ijk) <- paste0("VCov(~", rn.i.1, ",", "~", rn.i.2, ")")
        rownames(vp.ijk)[1] <- paste0(grp.titles[j], ": ", rownames(vp.ijk)[1])
        rownames(vp.ijk) <- format(rownames(vp.ijk), justify = "right")
        colnames(vp.ijk) <- ttl
        ii <- c(which(jk.1 == jk.2), which(jk.1 < jk.2))
        ii <- which(jk.1 <= jk.2)
        vp.ijk <- vp.ijk[ii, , drop = FALSE]
        vp.i_[[j, k]] <- vp.ijk
      }
    }
    vp.i_ <- lapply(1:n.eq, function(j) do.call(rbind, vp.i_[, j]))

    vp.i <- list()
    # vp.i <- array(NA,c(dim(vp.i_[[1]]),n.eq),dimnames=c(dimnames(vp.i_[[1]]),list(grp.titles)))
    vp.i <- array(NA, c(dim(vp.i_[[1]]), n.eq), dimnames = c(dimnames(vp.i_[[1]]), list(NULL)))
    for (j in 1:n.eq) {
      vp.i[, , j] <- vp.i_[[j]]
    }
    VarPar <- c(VarPar, structure(list(vp.i), names = lv.i))
  }


  phi <- smry$phi
  LR <- smry$null.deviance - smry$deviance
  df <- obj$model.df
  deviance <- deviance(obj)

  if (df > 0) {
    p <- pchisq(LR, df, lower.tail = FALSE)
    L0.pwr <- exp(-smry$null.deviance / N)
    LM.pwr <- exp(-smry$deviance / N)

    McFadden <- 1 - smry$deviance / smry$null.deviance
    Cox.Snell <- 1 - exp(-LR / N)
    Nagelkerke <- Cox.Snell / (1 - L0.pwr)
  } else {
    LR <- NA
    df <- NA
    p <- NA
    McFadden <- NA
    Cox.Snell <- NA
    Nagelkerke <- NA
  }

  ll <- obj$ll
  AIC <- AIC(obj)
  BIC <- AIC(obj, k = log(N))
  sumstat <- c(
    phi = phi,
    LR = LR,
    df = df,
    # p             = p,
    logLik = ll,
    deviance = deviance,
    McFadden = McFadden,
    Cox.Snell = Cox.Snell,
    Nagelkerke = Nagelkerke,
    AIC = AIC,
    BIC = BIC,
    N = N
  )

  ans <- list(coef = coef)
  ans <- c(ans, VarPar)
  parameter.types <- c("coef", names(VarPar))

  if(length(smry$ngrps)){
      G <-as.integer(smry$ngrps)
      names(G) <- names(smry$ngrps)
      names(G) <- paste("Groups by",names(G))
      G <- c(G,"Total obs."=N)

      sumstat <- list(sumstat,N=G)
      
      c(ans,
        list(sumstat=sumstat,
             parameter.types=parameter.types,
             call=obj$call,
             contrasts = obj$contrasts,
             xlevels = obj$xlevels))     
  }
  else {

    sumstat <- c(sumstat,N=N)
    c(ans,
    list(sumstat=sumstat,
         call=obj$call,
         contrasts = obj$contrasts,
         xlevels = obj$xlevels))     
  }  

}

getSummary.mmblogit <- getSummary.mblogit
