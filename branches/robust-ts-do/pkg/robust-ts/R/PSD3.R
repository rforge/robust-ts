PSD3 <- function (r, f,  finv, delta=0.05, ...) 
{
######################################################
##
##  R-function 'PSD3' - non-linear shrinking:
##  computes a positiv semidefinite autocorrelation function 
##  via the shrinking method (cf. Rousseeuw and Molenberghs, 1993),
##  using implementations (function 'is.positive.definite')
##  of R package 'corpcor' (cf. Higham, 1988),
##  cf. also Devlin et al. (1975)
##
##  author: Bernhard Spangl
##  date: 2008-06-11
##
##  r ... originial autocorrelation function
##  f ... monotone increasing continuous function
##  finv ... inverse of 'f'
##  delta ... shifting parameter (e.g.: delta=0.05)
##  ... ... tolerance (e.g.: tol=1e-08)
##
######################################################

    r <- as.vector(r)
    R <- toeplitz(r)

    was.psd <- is.positive.definite(R, method="eigen", ...)
    is.psd <- was.psd

    if (!was.psd) {
        limit <- finv(delta)
        while (!is.psd) {
            r <- sapply(r[-1], Trans, f=f, finv=finv, delta=delta, limit=limit)
            r <- c(1, r)
            R.psd <- toeplitz(r)
            is.psd <- is.positive.definite(R.psd, method="eigen", ...)
        }
        res <- list(r.psd=R.psd[1, ], R.psd=R.psd, R=R,
                    was.psd=was.psd, is.psd=is.psd)
    } else {
        res <- list(r.psd=r, R.psd=R, R=R,
                    was.psd=was.psd, is.psd=is.psd)
    }

    return(res)

}

Trans <- function (r, f, finv, delta, limit)
{
    if (r < -limit) {
        finv(f(r) + delta)
    } else {
        if (r <= limit) {
            0
        } else {
            finv(f(r) - delta)
        }
    }
}
