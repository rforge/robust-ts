\name{PSD3}
\Rdversion{1.1}
\alias{PSD3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PSD3(r, f, finv, delta = 0.05, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{f}{
%%     ~~Describe \code{f} here~~
}
  \item{finv}{
%%     ~~Describe \code{finv} here~~
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (r, f, finv, delta = 0.05, ...) 
{
    r <- as.vector(r)
    R <- toeplitz(r)
    was.psd <- is.positive.definite(R, method = "eigen", ...)
    is.psd <- was.psd
    if (!was.psd) {
        limit <- finv(delta)
        while (!is.psd) {
            r <- sapply(r[-1], Trans, f = f, finv = finv, delta = delta, 
                limit = limit)
            r <- c(1, r)
            R.psd <- toeplitz(r)
            is.psd <- is.positive.definite(R.psd, method = "eigen", 
                ...)
        }
        res <- list(r.psd = R.psd[1, ], R.psd = R.psd, R = R, 
            was.psd = was.psd, is.psd = is.psd)
    }
    else {
        res <- list(r.psd = r, R.psd = R, R = R, was.psd = was.psd, 
            is.psd = is.psd)
    }
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
