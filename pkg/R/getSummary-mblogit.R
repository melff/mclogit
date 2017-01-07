getSummary.mblogit <- function(obj,
                               alpha=.05,
                               ...){
  
  smry <- summary(obj)
  N <- obj$N
  coef <- smry$coefficients
  
  lower.cf <- qnorm(p=alpha/2,mean=coef[,1],sd=coef[,2])
  upper.cf <- qnorm(p=1-alpha/2,mean=coef[,1],sd=coef[,2])
  coef <- cbind(coef,lower.cf,upper.cf)
  colnames(coef) <- c("est","se","stat","p","lwr","upr")

  modcat <- colnames(obj$D)
  basecat <- rownames(obj$D)[rownames(obj$D)%nin%modcat]
  
  eqs <- paste0(modcat,"~")
  
  coef.grps <- lapply(eqs,function(eq){
    
    ii <- grep(eq,rownames(coef),fixed=TRUE)
    coef.grp <- coef[ii,]
    rownames(coef.grp) <- gsub(eq,"",rownames(coef.grp),fixed=TRUE)
    coef.grp
  })
  
  if(getOption("mblogit.show.basecat",TRUE))
    grp.titles <- paste(modcat,basecat,sep=getOption("mblogit.basecat.sep","/"))
  else
    grp.titles <- modcat
  
  names(coef.grps) <- grp.titles
  coef <- do.call(memisc::collect,coef.grps)
  
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
    N             = N
  )
  list(call=obj$call,
       coef=coef,
       contrasts = obj$contrasts,
       xlevels = obj$xlevels,
       sumstat=sumstat)
}

getSummary.mmblogit <- getSummary.mblogit
