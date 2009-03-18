\name{ACF}
\Rdversion{1.1}
\alias{ACF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ACF(x, type = "spearman", cov = FALSE, modified = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{cov}{
%%     ~~Describe \code{cov} here~~
}
  \item{modified}{
%%     ~~Describe \code{modified} here~~
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
function (x, type = "spearman", cov = FALSE, modified = TRUE) 
{
    n <- length(x)
    mlag <- n - 2
    acovf <- rep(NA, mlag + 1)
    acf <- rep(NA, mlag + 1)
    for (i in 0:mlag) {
        x1 <- x[(1 + i):n]
        x2 <- x[1:(n - i)]
        if (type == "spearman") {
            acf[i + 1] <- Spearman(x1, x2, n)
        }
        if (type == "Qn") {
            if (modified) {
                res <- CovQn(x1, x2, cov, n)
                acf[i + 1] <- res$acf.m
                if (cov) {
                  acovf[i + 1] <- res$acovf.m
                }
            }
            else {
                res <- CovQn(x1, x2, cov)
                acf[i + 1] <- res$acf
                if (cov) {
                  acovf[i + 1] <- res$acovf
                }
            }
        }
    }
    return(list(ACF = c(acf, 0), ACovF = c(acovf, 0)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
