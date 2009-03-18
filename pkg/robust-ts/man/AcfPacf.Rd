\name{AcfPacf}
\Rdversion{1.1}
\alias{AcfPacf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
AcfPacf(x, mx, type = "spearman")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{mx}{
%%     ~~Describe \code{mx} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
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
function (x, mx, type = "spearman") 
{
    if (missing(mx)) {
        mx <- mean(x)
    }
    x <- x - mx
    n <- length(x)
    mlag <- n - 2
    rho <- rep(NA, n)
    A <- matrix(NA, n - 2, n - 2)
    cor <- switch(type, pearson = get("cor", mode = "function"), 
        spearman = get("Spearman", mode = "function"))
    rho[1] <- 1
    A[1, 1] <- cor(x[1:(n - 1)], x[2:n], n)
    rho[2] <- A[1, 1]
    for (h in 2:mlag) {
        X <- embed(x, h + 1)
        u <- X[, 1] - X[, 2:h] \%*\% t(A[h - 1, 1:(h - 1), drop = F])
        v <- X[, h + 1] - X[, h:2] \%*\% t(A[h - 1, 1:(h - 1), 
            drop = F])
        A[h, h] <- cor(u, v, n)
        for (i in 1:(h - 1)) {
            A[h, i] <- A[h - 1, i] - A[h, h] * A[h - 1, h - i]
        }
        rho[h + 1] <- A[h - 1, 1:(h - 1)] \%*\% rev(rho[2:h]) + 
            A[h, h] * (1 - A[h - 1, 1:(h - 1)] \%*\% rho[2:h])
    }
    rho[n] <- 0
    return(list(ACF = rho, A = A))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
