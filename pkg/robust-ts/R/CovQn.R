CovQn <- function (u, v, cov=FALSE, N) 
{
######################################################
##
##  R-function 'CovQn'
##  computes the covariance and correlation coefficient
##  using the highly robust scale estimator 'Qn', to be
##  applied to time series (as non-parametric alternative
##  to the autocovariance and autocorrelation coefficient)
##  (cf. Ma and Genton, 2000)
##
##  IMPORTANT: requires R-package 'robustbase' 
##
##  author: Bernhard Spangl
##  date: 2008-06-01
##
##  u, v ... appropriate vectors
##  cov ... boolean, should the autocovariance also be
##          calculated? 
##  N ... length of time series for modified estimates
##
######################################################

    acf <- acovf <- NA
    acf.m <- acovf.m <- NA
    n <- length(u)
    if (missing(N)) {
        N <- n
    }
    mod.fact <- n/N

    Q1 <- Qn(u+v)^2
    Q2 <- Qn(u-v)^2

    acf <- (Q1 - Q2)/(Q1 + Q2)
    acf.m <- mod.fact*acf
    if (cov) {
        acovf <- (Q1 - Q2)/4
        acovf.m <- mod.fact*acovf
    }

    return(list(acf=acf, acf.m=acf.m, acovf=acovf, acovf.m=acovf.m))

}
