Spearman <- function (x, y, N) 
{
######################################################
##
##  R-function 'Spearman'
##  computes Spearman's (modified) rank correlation 
##  coefficient, to be applied to time series (as 
##  non-parametric alternative to the autocorrelation
##  coefficient)
##  (cf. Ahdesmäki et al., 2005, Equation 12)
##
##  author: Bernhard Spangl
##  date: 2007-04-06
##
##  x, y ... vectors
##  N ... interger, should biased or unbiased version be used?
##
######################################################

    if (missing(N)) {
        N <- length(x)
    }
    n <- length(x)
    rx <- rank(x)
    ry <- rank(y)
##  rho <- cor(x, y, method="spear")*n/N    # build-in version
    rho <- 12/(N*(n^2 - 1))*((rx - (n + 1)/2)%*%(ry - (n + 1)/2))
    return(drop(rho))
}
