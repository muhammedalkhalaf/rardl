\name{print.rardl}
\alias{print.rardl}
\alias{summary.rardl}
\title{Print and Summary Methods for rardl Objects}
\usage{
\method{print}{rardl}(x, digits = 4L, ...)

\method{summary}{rardl}(object, ...)
}
\arguments{
\item{x}{An object of class \code{"rardl"}.}

\item{digits}{Integer. Significant digits.}

\item{...}{Further arguments (ignored).}

\item{object}{An object of class \code{"rardl"}.}
}
\value{
Invisibly returns the input object.
}
\description{
Print and summary methods for objects of class \code{"rardl"} returned
by \code{\link{rardl}}.
}
\examples{
set.seed(1)
n  <- 60
df <- data.frame(y = cumsum(rnorm(n)), x1 = cumsum(rnorm(n)))
\donttest{
res <- rardl(df, type = "recursive", depvar = "y",
             indepvars = "x1", init_obs = 25, max_lag = 1,
             nsim = 50)
print(res)
summary(res)
}
}
