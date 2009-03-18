ACF <- function (x, type="spearman", cov=FALSE, modified=TRUE) 
{
######################################################
##
##  R-function 'ACF'
##  computes the autocorrelation function/sequence 
##  using Spearman's rank correlation coefficient
##  or highly robust autocovariance function proposed
##  by Ma & Genton (2000), 
##  (non-parametric alternative to the autocorrelation
##  coefficient)
##
##  author: Bernhard Spangl
##  date: 2008-06-01, 2007-04-06
##
##  x ... time series (vector)
##  type ... which correlation coefficient should be used,
##           Spearman's rho ("spearman") or Qn proposed by 
##           Rousseeuw & Croux (1993) ("Qn"),
##           default = "spearman"
##  cov ... should the autocovariance function also
##          be computed in case of using 'type="Qn"'
##  modified ... should the modified autocovariance function
##               be computed in case of using 'type="Qn"'
##
######################################################

    n <- length(x)
    mlag <- n-2
    acovf <- rep(NA, mlag+1)
    acf <- rep(NA, mlag+1)

##    cor <- switch(type, spearman=get("Spearman", mode="function"),
##                  Qn=get("CovQn", mode="function"))

    for (i in 0:mlag) {
        x1 <- x[(1+i):n]
        x2 <- x[1:(n-i)]
        if (type == "spearman") {    # Spearman's rho
            acf[i+1] <- Spearman(x1, x2, n)
        }
        if (type == "Qn") {    # highly robust scale estimator 'Qn'
            if (modified) {
                res <- CovQn(x1, x2, cov, n)
                acf[i+1] <- res$acf.m
                if (cov) {
                    acovf[i+1] <- res$acovf.m
                }
            } else {
                res <- CovQn(x1, x2, cov)
                acf[i+1] <- res$acf
                if (cov) {
                    acovf[i+1] <- res$acovf
                }
            }
        }
    }
    return(list(ACF=c(acf, 0), ACovF=c(acovf, 0)))
}
