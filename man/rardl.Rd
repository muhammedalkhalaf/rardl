\name{rardl}
\alias{rardl}
\title{Rolling-Window and Recursive ARDL Cointegration Analysis}
\usage{
rardl(data = NULL, type,
      depvar = NULL, indepvars = NULL,
      window_size = NULL, init_obs = 30L,
      max_lag = 4L, ic = "bic",
      case = 3L, level = 5,
      nsim = 1000L, seed = -1L,
      tobs = c(60L, 120L, 250L, 500L),
      nregs = 1L,
      adf_case = "c",
      transform = NULL)
}
\arguments{
\item{data}{A data frame. First column is the dependent variable for
\code{"rolling"}, \code{"recursive"}, and \code{"rgranger"} types.}

\item{type}{Character. One of \code{"rolling"}, \code{"recursive"},
\code{"radf"}, \code{"rgranger"}, \code{"simulate"}.}

\item{depvar}{Character. Name of the dependent variable column.}

\item{indepvars}{Character vector. Names of independent variable columns.}

\item{window_size}{Integer or \code{NULL}. Rolling window width.}

\item{init_obs}{Integer. Initial observations for recursive methods.
Default \code{30}.}

\item{max_lag}{Non-negative integer. Maximum lag for selection.
Default \code{4}.}

\item{ic}{Character. Information criterion: \code{"bic"} (default) or
\code{"aic"}.}

\item{case}{Integer 1--5. PSS deterministic case. Default \code{3}.}

\item{level}{Numeric. Significance level: \code{1}, \code{5} (default),
or \code{10}.}

\item{nsim}{Positive integer. Monte Carlo replications. Default \code{1000}.}

\item{seed}{Integer. Random seed (\code{-1} for no fixed seed).}

\item{tobs}{Integer vector. Sample sizes for \code{"simulate"}.}

\item{nregs}{Positive integer. Number of regressors for \code{"simulate"}.}

\item{adf_case}{Character. ADF model: \code{"c"} (default), \code{"ct"},
or \code{"n"}.}

\item{transform}{Character or \code{NULL}. Pre-transformation:
\code{"log"} or \code{"diff"}. Default \code{NULL}.}
}
\value{
An object of class \code{"rardl"} with components:
\describe{
  \item{type}{The analysis type.}
  \item{results}{Data frame of window/period-specific test statistics.}
  \item{critical_values}{Named numeric vector of simulated critical values.}
  \item{depvar}{Dependent variable name.}
  \item{indepvars}{Independent variable names.}
  \item{nobs}{Total observations.}
  \item{window_size}{Window size (rolling type only).}
  \item{init_obs}{Initial observations (recursive types only).}
}
}
\description{
Provides rolling-window and recursive implementations of the ARDL bounds
cointegration test. Five analysis types detect structural instability in
cointegrating relationships over time.
}
\examples{
set.seed(42)
n  <- 80
x1 <- cumsum(rnorm(n))
y  <- 0.5 * x1 + rnorm(n, sd = 0.5)
df <- data.frame(y = y, x1 = x1)
\donttest{
res <- rardl(df, type = "recursive", depvar = "y",
             indepvars = "x1", init_obs = 30, max_lag = 2,
             nsim = 100)
print(res)
}
}
\references{
Shahbaz, M., Khan, A. and Mubarak, M. S. (2023). Do Geopolitical Risks
Obstruct the Chinese Economic Progress? Evidence from Rolling-Window ARDL.
\emph{Energy Economics}, 116, 106642.
\doi{10.1016/j.eneco.2022.106642}

Khan, A., Shahbaz, M. and Napari, A. (2023). Recursive ARDL: Detecting
Unstable Cointegration. \emph{Energy Economics}, 117, 106683.
\doi{10.1016/j.eneco.2022.106683}

Pesaran, M. H., Shin, Y. and Smith, R. J. (2001). Bounds testing
approaches to the analysis of level relationships.
\emph{Journal of Applied Econometrics}, 16(3), 289--326.
\doi{10.1002/jae.616}
}
