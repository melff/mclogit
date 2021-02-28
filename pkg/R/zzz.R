.onLoad <- function(lib,pkg){
  if(requireNamespace("memisc",quietly = TRUE)){
    
    memisc::setSummaryTemplate(
      mclogit = c(
        "Likelihood-ratio" = "($LR:f1#)",
        #p             = "($p:#)",
        "Log-likelihood" = "($logLik:f1#)",
        Deviance      = "($deviance:f1#)",
        AIC           = "($AIC:f1#)",
        BIC           = "($BIC:f1#)",
        N             = "($N:d)"
      ),
      mmclogit = c(
        #"Likelihood-ratio" = "($LR:f1#)",
        #p             = "($p:#)",
        #"Log-likelihood" = "($logLik:f1#)",
        Deviance      = "($deviance:f1#)",
        #AIC           = "($AIC:f1#)",
        #BIC           = "($BIC:f1#)",
        N             = "($N:d)"
      ),
      mblogit = c(
        "Log-likelihood" = "($logLik:f1#)",
        Deviance      = "($deviance:f1#)",
        AIC           = "($AIC:f1#)",
        BIC           = "($BIC:f1#)",
        N             = "($N:d)"
      )
    )    
  }
  
  options(mblogit.basecat.sep="/")
  options(mblogit.show.basecat=TRUE)
  options(summary.stats.mclogit=c("Deviance","N"))
  options(summary.stats.mmclogit=c("Deviance","N"))
}

