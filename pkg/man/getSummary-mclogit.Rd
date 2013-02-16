\name{getSummary.mclogit}
\alias{getSummary.mclogit}
\title{getSummary.mclogit}
\description{
  A \code{\link[memisc]{getSummary}} method for use by \code{\link[memisc]{mtable}}
}
\usage{
\method{getSummary}{mclogit}(obj,
            alpha=.05,
            rearrange=NULL,
            \dots)
}
\arguments{
  \item{obj}{an object returned by \code{\link{mclogit}}}
  \item{alpha}{level of the confidence intervals; their coverage should
    be 1-alpha/2 }
  \item{rearrange}{an optional named list of character vectors.
      Each element of the list designates a column in the table of estimates,
      and each element of a character vector refers to a coefficient.
      Names of list elements become column heads and names
      of the character vector elements become coefficient labels.
      }
  \item{\dots}{further arguments; ignored.}
}
\examples{
\dontrun{
summary(classd.model <- mclogit(cbind(Freq,choice.set)~
                   (econdim1.sq+nonmatdim1.sq+nonmatdim2.sq)+
                   (econdim1+nonmatdim1+nonmatdim2)+
                   (econdim1+nonmatdim1+nonmatdim2):classd,
                  data=mvoteint.classd,random=~1|mvoteint/eb,
                  subset=classd!="Farmers"))
myGetSummary.classd <- function(x)getSummary.mclogit(x,rearrange=list(
        "Econ. Left/Right"=c(
                    "Squared effect"="econdim1.sq",
                    "Linear effect"="econdim1",
                    " x Intermediate/Manual worker"="econdim1:classdIntermediate",
                    " x Service class/Manual worker"="econdim1:classdService class",
                    " x Self-employed/Manual worker"="econdim1:classdSelf-employed"
                    ),
        "Lib./Auth."=c(
                    "Squared effect"="nonmatdim1.sq",
                    "Linear effect"="nonmatdim1",
                    " x Intermediate/Manual worker"="nonmatdim1:classdIntermediate",
                    " x Service class/Manual worker"="nonmatdim1:classdService class",
                    " x Self-employed/Manual worker"="nonmatdim1:classdSelf-employed"
                    ),
        "Mod./Trad."=c(
                    "Squared effect"="nonmatdim2.sq",
                    "Linear effect"="nonmatdim2",
                    " x Intermediate/Manual worker"="nonmatdim2:classdIntermediate",
                    " x Service class/Manual worker"="nonmatdim2:classdService class",
                    " x Self-employed/Manual worker"="nonmatdim2:classdSelf-employed"
                    )
        ))

mtable(classd.model,getSummary=myGetSummary.classd)
# Output would look like so:
# ==================================================================================
#                                 Econ. Left/Right    Lib./Auth.       Mod./Trad.
# ----------------------------------------------------------------------------------
# Squared effect                      0.030            0.008           -0.129**
#                                    (0.081)          (0.041)          (0.047)
# Linear effect                      -0.583***        -0.038            0.137**
#                                    (0.063)          (0.041)          (0.045)
#  x Intermediate/Manual worker       0.632***        -0.029           -0.015
#                                    (0.026)          (0.020)          (0.019)
#  x Service class/Manual worker      1.158***         0.084**          0.000
#                                    (0.040)          (0.032)          (0.030)
#  x Self-employed/Manual worker      1.140***         0.363***         0.112***
#                                    (0.035)          (0.027)          (0.026)
# Var(mvoteint)                       1.080***
#                                    (0.000)
# Var(mvoteint x eb)                  0.118***
#                                    (0.000)
# ----------------------------------------------------------------------------------
# Dispersion                              1.561
# Deviance                            15007.0
# N                                  173445
# ==================================================================================
}
}