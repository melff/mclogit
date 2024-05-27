getSummary.mclogit <- function(obj,
                               alpha=.05,
                               rearrange=NULL,
                               ...){
  
  smry <- summary(obj)
  N <- obj$N
  coef <- smry$coefficients
  varPar <- smry$varPar
  
  lower.cf <- qnorm(p=alpha/2,mean=coef[,1],sd=coef[,2])
  upper.cf <- qnorm(p=1-alpha/2,mean=coef[,1],sd=coef[,2])
  coef <- cbind(coef,lower.cf,upper.cf)
  colnames(coef) <- c("est","se","stat","p","lwr","upr")
  if(length(rearrange)){
    coef.grps <- lapply(rearrange,function(ii){
      if(is.character(ii) && !all(ii %in% rownames(coef)))
        stop("coefficient(s) ",dQuote(unname(ii[!(ii %in% rownames(coef))]))," do not exist")
      structure(coef[ii,],
                dimnames=list(names(ii),dimnames(coef)[[2]])
      )
    })
    grp.titles <- names(rearrange)
    coef.grps <- do.call(memisc::collect,coef.grps)
    coef <- array(NA,dim=c(
      dim(coef.grps)[1] + NROW(varPar),
      dim(coef.grps)[2],
      dim(coef.grps)[3]
    ))
    coef[seq(dim(coef.grps)[1]),,] <- coef.grps
    if(length(varPar))
      coef[dim(coef.grps)[1]+seq(nrow(varPar)),,1] <- varPar
    dimnames(coef) <- list(
      c(dimnames(coef.grps)[[1]],rownames(varPar)),
      dimnames(coef.grps)[[2]],
      grp.titles
    )
  }
  
  VarPar <- NULL
  VarCov <- smry$VarCov
  se_VarCov <- smry$se_VarCov

  for(i in seq_along(VarCov)){
    lv.i <- names(VarCov)[i]
    vc.i <- VarCov[[i]]
    vr.i <- diag(vc.i)
    cv.i <- vc.i[lower.tri(vc.i)]
    se_vc.i <- se_VarCov[[i]]
    se_vr.i <- diag(se_vc.i)
    se_cv.i <- se_vc.i[lower.tri(se_vc.i)]   
    nms.i <- rownames(vc.i)
    nms.i <- gsub("(Intercept)","1",nms.i,fixed=TRUE)
    vrnames.i <- paste0("Var(~",nms.i,"|",lv.i,")")
    cvnames.i <- t(outer(nms.i,nms.i,FUN=paste,sep=":"))
    cvnames.i <- cvnames.i[lower.tri(cvnames.i)]
    if(length(cvnames.i))
      cvnames.i <- paste0("Cov(~",cvnames.i,"|",lv.i,")")
    vp.i <- matrix(NA,nrow=length(vr.i)+length(cv.i),ncol=6)
    vp.i[,1] <- c(vr.i,cv.i)
    vp.i[,2] <- c(se_vr.i,se_cv.i)
    dimnames(vp.i) <- list(c(vrnames.i,cvnames.i),
                           c("est","se","stat","p","lwr","upr"))
    VarPar <- c(VarPar,structure(list(vp.i),names=lv.i))
  }
   
  phi <- smry$phi
  LR <- smry$null.deviance - smry$deviance
  df <- obj$model.df
  deviance <- deviance(obj)
  
  
  if(df > 0){
    p <- pchisq(LR,df,lower.tail=FALSE)
    L0.pwr <- exp(-smry$null.deviance/N)
    LM.pwr <- exp(-smry$deviance/N)
    
    McFadden <- 1- smry$deviance/smry$null.deviance
    Cox.Snell <- 1 - exp(-LR/N)
    Nagelkerke <- Cox.Snell/(1-L0.pwr)
  }
  else {
    LR <- NA
    df <- NA
    p <- NA
    McFadden <- NA
    Cox.Snell <- NA
    Nagelkerke <- NA
  }
  
  ll <- obj$ll
  AIC <- AIC(obj)
  BIC <- AIC(obj,k=log(N))
  sumstat <- c(
    phi         = phi,
    LR             = LR,
    df         = df,
    #p             = p,
    logLik        = ll,
    deviance      = deviance,
    McFadden      = McFadden,
    Cox.Snell       = Cox.Snell,
    Nagelkerke    = Nagelkerke,
    AIC           = AIC,
    BIC           = BIC,
    N = N
  )

  ans <- list(coef= coef)
  ans <- c(ans,VarPar)
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

getSummary.mmclogit <- getSummary.mclogit
