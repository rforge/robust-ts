AcfPacf <- function (x, mx, type="spearman")
{
######################################################
##
##  R-function 'AcfPacf':
##  computes the autocorrelation cofficients recursively 
##  via the partial autocorrelation cofficients
##  (cf. Möttönen et al., 1999)
##
##  author: Bernhard Spangl
##  date: 2007-04-09
##
##  x ... time series (vector)
##  mx ... mean of time series x
##  type ... which correlation coefficient should be used,
##           Pearson's or Spearman's? (default = spearman)
##
######################################################

    if (missing(mx)) {
        mx <- mean(x)
    }
    x <- x - mx    # centering 
    n <- length(x)
    mlag <- n - 2     # maximum lag using this recursive algorithm
    rho <- rep(NA, n)
    A <- matrix(NA, n-2, n-2)

    ##  Doesn't work for 'pearson', because the biased estimator for
    ##  the Pearson's correlation coefficient has to be used (i.e.,
    ##  should always be divided by length of series)!
    cor <- switch(type, pearson=get("cor", mode="function"), 
                  spearman=get("Spearman", mode="function"))
    
    rho[1] <- 1     # lag zero
    
    A[1, 1] <- cor(x[1:(n-1)], x[2:n], n)
##  A[1, 1] <- pacf(x, lag.max=n-1, plot=F)$acf[1, 1, 1]    # built-in
    rho[2] <- A[1, 1]
    
    for (h in 2:mlag) {
        X <- embed(x, h+1)
        u <- X[, 1] - X[, 2:h]%*%t(A[h-1, 1:(h-1), drop=F])
        v <- X[, h+1] - X[, h:2]%*%t(A[h-1, 1:(h-1), drop=F])
        A[h, h] <- cor(u, v, n)
##      A[h, h] <- pacf(x, lag.max=n-1, plot=F)$acf[h, 1, 1]    # built-in
        for (i in 1:(h-1)) {
            A[h, i] <- A[h-1, i] - A[h, h]*A[h-1, h-i]
        }
        rho[h+1] <- A[h-1, 1:(h-1)]%*%rev(rho[2:h]) +
                    A[h, h]*(1 - A[h-1, 1:(h-1)]%*%rho[2:h])
    }

    rho[n] <- 0 
    return(list(ACF=rho, A=A))
    
}
