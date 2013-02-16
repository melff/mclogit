.onLoad <- function(lib,pkg){
    setSummaryTemplate(
      mclogit = c(
              "McFadden R-sq." = "($McFadden:f#)",
              "Cox-Snell R-sq." = "($Cox.Snell:f#)",
              "Nagelkerke R-sq."  = "($Nagelkerke:f#)",
              Dispersion         = "($phi:#)",
              "Likelihood-ratio" = "($LR:f1#)",
              #p             = "($p:#)",
              "Log-likelihood" = "($logLik:f1#)",
              Deviance      = "($deviance:f1#)",
              AIC           = "($AIC:f1#)",
              BIC           = "($BIC:f1#)",
              N             = "($N:d)"
      ),
      mclogitRandeff = c(
              #"McFadden R-sq." = "($McFadden:f#)",
              #"Cox-Snell R-sq." = "($Cox.Snell:f#)",
              #"Nagelkerke R-sq."  = "($Nagelkerke:f#)",
              Dispersion         = "($phi:#)",
              #"Likelihood-ratio" = "($LR:f1#)",
              #p             = "($p:#)",
              #"Log-likelihood" = "($logLik:f1#)",
              Deviance      = "($deviance:f1#)",
              #AIC           = "($AIC:f1#)",
              #BIC           = "($BIC:f1#)",
              N             = "($N:d)"
      )
    )
}

