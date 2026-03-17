#' Rolling-Window and Recursive ARDL Cointegration Analysis
#'
#' Implements five analysis types for time-varying ARDL cointegration and
#' unit root testing:
#' \describe{
#'   \item{\code{"rolling"}}{Rolling-window ARDL bounds test.}
#'   \item{\code{"recursive"}}{Recursive (expanding-window) ARDL bounds test.}
#'   \item{\code{"radf"}}{Recursive ADF unit root test.}
#'   \item{\code{"rgranger"}}{Recursive Granger causality test.}
#'   \item{\code{"simulate"}}{Monte Carlo simulation of critical values.}
#' }
#'
#' @param data A data frame. The first column is treated as the dependent
#'   variable for \code{"rolling"}, \code{"recursive"}, and
#'   \code{"rgranger"} types. For \code{"radf"}, only one column is needed.
#'   For \code{"simulate"}, \code{data} may be \code{NULL}.
#' @param type Character. One of \code{"rolling"}, \code{"recursive"},
#'   \code{"radf"}, \code{"rgranger"}, \code{"simulate"}.
#' @param depvar Character. Name of the dependent variable column in
#'   \code{data}.
#' @param indepvars Character vector. Names of independent variable columns.
#'   Required for \code{"rolling"} and \code{"recursive"}.
#' @param window_size Integer or \code{NULL}. Rolling window width for
#'   \code{"rolling"}. Defaults to approximately \code{floor(0.3 * T)}.
#' @param init_obs Integer. Initial number of observations for recursive
#'   methods. Default \code{30}.
#' @param max_lag Non-negative integer. Maximum lag for ARDL / ADF selection.
#'   Default \code{4}.
#' @param ic Character. Information criterion: \code{"bic"} (default) or
#'   \code{"aic"}.
#' @param case Integer 1--5. PSS deterministic case. Default \code{3}.
#' @param level Numeric (0, 100). Significance level for critical value
#'   display. One of \code{1}, \code{5} (default), or \code{10}.
#' @param nsim Positive integer. Monte Carlo simulation replications.
#'   Default \code{1000}.
#' @param seed Integer. Random seed (\code{-1} for no fixed seed). Default
#'   \code{-1}.
#' @param tobs Integer vector. Sample sizes for \code{"simulate"}. Defaults
#'   to \code{c(60, 120, 250, 500)}.
#' @param nregs Positive integer. Number of regressors for \code{"simulate"}.
#'   Default \code{1}.
#' @param adf_case Character. ADF model for \code{"radf"}: \code{"c"}
#'   (constant, default), \code{"ct"} (constant + trend), or \code{"n"}
#'   (none).
#' @param transform Character or \code{NULL}. Optional pre-transformation:
#'   \code{"log"} or \code{"diff"}. Default \code{NULL}.
#'
#' @return An object of class \code{"rardl"} with components depending on
#'   \code{type}:
#'   \describe{
#'     \item{type}{The analysis type used.}
#'     \item{results}{A data frame of window/period-specific test statistics.}
#'     \item{critical_values}{Named numeric vector of simulated or tabulated
#'       critical values.}
#'     \item{depvar}{Dependent variable name.}
#'     \item{indepvars}{Independent variable names.}
#'     \item{nobs}{Total number of observations in the full sample.}
#'     \item{window_size}{Window size used (rolling type only).}
#'     \item{init_obs}{Initial observation count (recursive types only).}
#'   }
#'
#' @references
#' Shahbaz, M., Khan, A. and Mubarak, M. S. (2023). Do Geopolitical Risks
#' Obstruct the Chinese Economic Progress? Evidence from Rolling-Window ARDL.
#' \emph{Energy Economics}, 116, 106642.
#' \doi{10.1016/j.eneco.2022.106642}
#'
#' Khan, A., Shahbaz, M. and Napari, A. (2023). Recursive ARDL: Detecting
#' Unstable Cointegration. \emph{Energy Economics}, 117, 106683.
#' \doi{10.1016/j.eneco.2022.106683}
#'
#' Pesaran, M. H., Shin, Y. and Smith, R. J. (2001). Bounds testing approaches
#' to the analysis of level relationships. \emph{Journal of Applied Econometrics},
#' 16(3), 289--326. \doi{10.1002/jae.616}
#'
#' @examples
#' set.seed(42)
#' n  <- 80
#' x1 <- cumsum(rnorm(n))
#' y  <- 0.5 * x1 + rnorm(n, sd = 0.5)
#' df <- data.frame(y = y, x1 = x1)
#' \donttest{
#' res <- rardl(df, type = "recursive", depvar = "y",
#'              indepvars = "x1", init_obs = 30, max_lag = 2)
#' print(res)
#' }
#'
#' @importFrom stats lm coef vcov residuals logLik pf pt rnorm sd var
#'   AIC BIC quantile setNames
#' @export
rardl <- function(data = NULL, type,
                  depvar = NULL, indepvars = NULL,
                  window_size = NULL, init_obs = 30L,
                  max_lag = 4L, ic = "bic",
                  case = 3L, level = 5,
                  nsim = 1000L, seed = -1L,
                  tobs = c(60L, 120L, 250L, 500L),
                  nregs = 1L,
                  adf_case = "c",
                  transform = NULL) {

  ## --- Input validation ---
  type <- tolower(as.character(type))
  if (!type %in% c("rolling", "recursive", "radf", "rgranger", "simulate")) {
    stop("'type' must be one of: \"rolling\", \"recursive\", \"radf\",",
         " \"rgranger\", \"simulate\".", call. = FALSE)
  }
  ic       <- tolower(as.character(ic))
  adf_case <- tolower(as.character(adf_case))
  if (!ic %in% c("aic", "bic")) {
    stop("'ic' must be \"aic\" or \"bic\".", call. = FALSE)
  }
  if (!adf_case %in% c("c", "ct", "n")) {
    stop("'adf_case' must be \"c\", \"ct\", or \"n\".", call. = FALSE)
  }

  if (seed >= 0L) set.seed(seed)

  ## --- Dispatch ---
  if (type == "simulate") {
    return(.rardl_simulate(tobs = tobs, nregs = nregs, nsim = nsim,
                           max_lag = max_lag, ic = ic, case = case,
                           level = level))
  }

  ## For all non-simulate types, data is required
  if (is.null(data)) {
    stop("'data' must be provided for type = \"", type, "\".", call. = FALSE)
  }
  if (!is.data.frame(data)) data <- as.data.frame(data)

  ## Identify variables
  if (is.null(depvar)) {
    depvar <- names(data)[1L]
    message("'depvar' not specified. Using first column: ", depvar)
  }
  if (!depvar %in% names(data)) {
    stop(sprintf("'depvar' = \"%s\" not found in data.", depvar), call. = FALSE)
  }

  ## Optional transform
  if (!is.null(transform)) {
    data <- .rardl_transform(data, transform)
  }

  y_full <- as.numeric(data[[depvar]])
  n_full <- length(y_full)

  if (type == "radf") {
    return(.rardl_radf(y = y_full, init_obs = init_obs, max_lag = max_lag,
                       ic = ic, adf_case = adf_case, level = level,
                       nsim = nsim, depvar = depvar))
  }

  if (is.null(indepvars)) {
    indepvars <- setdiff(names(data), depvar)
    if (length(indepvars) == 0L) {
      stop("No independent variables found.", call. = FALSE)
    }
    message("'indepvars' not specified. Using: ",
            paste(indepvars, collapse = ", "))
  }
  for (iv in indepvars) {
    if (!iv %in% names(data)) {
      stop(sprintf("'indepvars' variable \"%s\" not found in data.", iv),
           call. = FALSE)
    }
  }
  x_full <- as.matrix(data[, indepvars, drop = FALSE])

  if (type == "rgranger") {
    if (length(indepvars) != 1L) {
      stop("'rgranger' type requires exactly one independent variable.",
           call. = FALSE)
    }
    return(.rardl_rgranger(y = y_full, x = x_full[, 1L],
                           init_obs = init_obs, max_lag = max_lag,
                           ic = ic, level = level, nsim = nsim,
                           depvar = depvar, indepvar = indepvars[1L]))
  }

  if (type == "rolling") {
    if (is.null(window_size)) {
      window_size <- max(floor(0.3 * n_full), init_obs)
    }
    window_size <- as.integer(window_size)
    if (window_size < 20L || window_size >= n_full) {
      stop("'window_size' must be between 20 and T - 1.", call. = FALSE)
    }
    return(.rardl_rolling(y = y_full, xmat = x_full,
                           window_size = window_size,
                           max_lag = max_lag, ic = ic, case = case,
                           level = level, nsim = nsim,
                           depvar = depvar, indepvars = indepvars))
  }

  if (type == "recursive") {
    return(.rardl_recursive(y = y_full, xmat = x_full,
                             init_obs = init_obs, max_lag = max_lag,
                             ic = ic, case = case, level = level,
                             nsim = nsim, depvar = depvar,
                             indepvars = indepvars))
  }
}


## ===========================================================================
## INTERNAL IMPLEMENTATIONS
## ===========================================================================

#' @keywords internal
.rardl_transform <- function(data, transform) {
  if (transform == "log") {
    for (nm in names(data)) {
      if (is.numeric(data[[nm]]) && all(data[[nm]] > 0, na.rm = TRUE)) {
        data[[nm]] <- log(data[[nm]])
      }
    }
  } else if (transform == "diff") {
    for (nm in names(data)) {
      if (is.numeric(data[[nm]])) {
        data[[nm]] <- c(NA_real_, diff(data[[nm]]))
      }
    }
    data <- data[!is.na(data[[1L]]), , drop = FALSE]
  }
  data
}


#' Compute ARDL ECM test statistics for a given window
#' @keywords internal
.rardl_ecm_stats <- function(y, xmat, max_lag, ic) {
  n   <- length(y)
  n_x <- ncol(xmat)

  best_ic <- Inf
  opt_p   <- 1L
  opt_q   <- rep(0L, n_x)

  for (p in 1L:min(max_lag, floor(n / 5L))) {
    for (q in 0L:min(max_lag, floor(n / 5L))) {
      reg <- .rardl_build_ecm_simple(y, xmat, p, rep(q, n_x))
      if (length(reg$y_dep) < ncol(reg$X_mat) + 5L) next
      fit <- tryCatch(lm(reg$y_dep ~ reg$X_mat - 1),
                      error = function(e) NULL)
      if (is.null(fit)) next
      n_u  <- length(reg$y_dep)
      k_u  <- length(coef(fit))
      ll_u <- as.numeric(logLik(fit))
      ic_v <- if (ic == "aic") -2 * ll_u + 2 * k_u else
        -2 * ll_u + k_u * log(n_u)
      if (ic_v < best_ic) {
        best_ic <- ic_v
        opt_p   <- p
        opt_q   <- rep(q, n_x)
      }
    }
  }

  reg_final <- .rardl_build_ecm_simple(y, xmat, opt_p, opt_q)
  if (length(reg_final$y_dep) < ncol(reg_final$X_mat) + 3L) {
    return(list(F_pss = NA_real_, t_pss = NA_real_, F_ind = NA_real_))
  }
  fit <- tryCatch(lm(reg_final$y_dep ~ reg_final$X_mat - 1),
                  error = function(e) NULL)
  if (is.null(fit)) {
    return(list(F_pss = NA_real_, t_pss = NA_real_, F_ind = NA_real_))
  }

  b   <- coef(fit)
  V   <- vcov(fit)
  nm  <- colnames(reg_final$X_mat)
  lv  <- grep("^(y_lag1|x_lag)", nm)
  ylv <- which(nm == "y_lag1")
  xlv <- grep("^x_lag", nm)

  F_pss <- tryCatch({
    R <- matrix(0, nrow = length(lv), ncol = length(b))
    for (i in seq_along(lv)) R[i, lv[i]] <- 1
    Rb <- R %*% b
    as.numeric(t(Rb) %*% solve(R %*% V %*% t(R)) %*% Rb) / length(lv)
  }, error = function(e) NA_real_)

  t_pss <- if (length(ylv) == 1L) {
    b[[ylv]] / sqrt(V[ylv, ylv])
  } else NA_real_

  F_ind <- tryCatch({
    if (length(xlv) == 0L) return(NA_real_)
    R <- matrix(0, nrow = length(xlv), ncol = length(b))
    for (i in seq_along(xlv)) R[i, xlv[i]] <- 1
    Rb <- R %*% b
    as.numeric(t(Rb) %*% solve(R %*% V %*% t(R)) %*% Rb) / length(xlv)
  }, error = function(e) NA_real_)

  list(F_pss = F_pss, t_pss = t_pss, F_ind = F_ind)
}


#' Simple ECM builder
#' @keywords internal
.rardl_build_ecm_simple <- function(y, xmat, p, q_vec) {
  n    <- length(y)
  n_x  <- ncol(xmat)
  ml   <- max(p, max(q_vec, 0L)) + 1L
  if (n <= ml + 5L) {
    return(list(y_dep = numeric(0),
                X_mat = matrix(nrow = 0L, ncol = 1L)))
  }
  dy   <- diff(y)
  idx  <- (ml + 1L):(n - 1L)
  if (length(idx) < 5L) {
    return(list(y_dep = numeric(0),
                X_mat = matrix(nrow = 0L, ncol = 1L)))
  }
  n_used <- length(idx)

  y_dep   <- dy[idx]
  y_lag1  <- y[idx]
  intcpt  <- rep(1, n_used)

  x_lag_mat <- xmat[idx, , drop = FALSE]
  colnames(x_lag_mat) <- paste0("x_lag", seq_len(n_x))

  ly_diff_mat <- matrix(NA_real_, nrow = n_used, ncol = p)
  for (j in seq_len(p)) ly_diff_mat[, j] <- dy[idx - j]
  colnames(ly_diff_mat) <- paste0("Ldy_", seq_len(p))

  x_diff_list <- list()
  for (i in seq_len(n_x)) {
    dx_i <- diff(xmat[, i])
    qi   <- q_vec[i]
    for (j in 0L:qi) {
      nm <- if (j == 0L) paste0("Dx_", i) else paste0("Ldx_", i, "_", j)
      x_diff_list[[nm]] <- dx_i[idx - j]
    }
  }
  x_diff_mat <- if (length(x_diff_list) > 0L) {
    do.call(cbind, x_diff_list)
  } else {
    matrix(nrow = n_used, ncol = 0L)
  }

  X_mat <- cbind(
    intercept = intcpt,
    y_lag1    = y_lag1,
    x_lag_mat,
    ly_diff_mat,
    x_diff_mat
  )
  list(y_dep = y_dep, X_mat = X_mat)
}


#' Rolling-window ARDL bounds test
#' @keywords internal
.rardl_rolling <- function(y, xmat, window_size, max_lag, ic, case,
                            level, nsim, depvar, indepvars) {
  n  <- length(y)
  n_wins <- n - window_size + 1L

  results <- data.frame(
    window_end = seq(window_size, n),
    F_pss      = NA_real_,
    t_pss      = NA_real_,
    F_ind      = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n_wins)) {
    idx <- i:(i + window_size - 1L)
    yw  <- y[idx]
    xw  <- xmat[idx, , drop = FALSE]
    st  <- .rardl_ecm_stats(yw, xw, max_lag = max_lag, ic = ic)
    results$F_pss[i] <- st$F_pss
    results$t_pss[i] <- st$t_pss
    results$F_ind[i] <- st$F_ind
  }

  ## Simulated critical values for window_size
  cvs <- .rardl_sim_cvs(n_obs = window_size, n_x = ncol(xmat),
                         nsim = nsim, max_lag = max_lag, ic = ic,
                         level = level)

  structure(
    list(type            = "rolling",
         results         = results,
         critical_values = cvs,
         depvar          = depvar,
         indepvars       = indepvars,
         nobs            = n,
         window_size     = window_size),
    class = "rardl"
  )
}


#' Recursive ARDL bounds test
#' @keywords internal
.rardl_recursive <- function(y, xmat, init_obs, max_lag, ic, case,
                              level, nsim, depvar, indepvars) {
  n <- length(y)

  results <- data.frame(
    obs    = seq(init_obs, n),
    F_pss  = NA_real_,
    t_pss  = NA_real_,
    F_ind  = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(results$obs)) {
    t_end <- results$obs[i]
    yr    <- y[1L:t_end]
    xr    <- xmat[1L:t_end, , drop = FALSE]
    st    <- .rardl_ecm_stats(yr, xr, max_lag = max_lag, ic = ic)
    results$F_pss[i] <- st$F_pss
    results$t_pss[i] <- st$t_pss
    results$F_ind[i] <- st$F_ind
  }

  ## Simulated critical values at full sample
  cvs <- .rardl_sim_cvs(n_obs = n, n_x = ncol(xmat),
                         nsim = nsim, max_lag = max_lag, ic = ic,
                         level = level)

  structure(
    list(type            = "recursive",
         results         = results,
         critical_values = cvs,
         depvar          = depvar,
         indepvars       = indepvars,
         nobs            = n,
         init_obs        = init_obs),
    class = "rardl"
  )
}


#' Recursive ADF unit root test
#' @keywords internal
.rardl_radf <- function(y, init_obs, max_lag, ic, adf_case, level,
                         nsim, depvar) {
  n <- length(y)

  results <- data.frame(
    obs     = seq(init_obs, n),
    adf_stat = NA_real_,
    opt_lag  = NA_integer_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(results$obs)) {
    t_end   <- results$obs[i]
    yr      <- y[1L:t_end]
    adf_res <- .rardl_adf_one(yr, max_lag = max_lag, ic = ic,
                               adf_case = adf_case)
    results$adf_stat[i] <- adf_res$stat
    results$opt_lag[i]  <- adf_res$opt_lag
  }

  ## Bootstrap critical values for ADF (standard tables, MacKinnon approx)
  cv_vals <- .rardl_adf_cvs(adf_case = adf_case, level = level)

  structure(
    list(type            = "radf",
         results         = results,
         critical_values = cv_vals,
         depvar          = depvar,
         nobs            = n,
         init_obs        = init_obs),
    class = "rardl"
  )
}


#' ADF test for one sample
#' @keywords internal
.rardl_adf_one <- function(y, max_lag, ic, adf_case) {
  n   <- length(y)
  dy  <- diff(y)
  n_d <- length(dy)

  best_ic  <- Inf
  opt_lag  <- 0L
  opt_stat <- NA_real_

  for (p in 0L:min(max_lag, floor(n / 4L))) {
    ml <- p + 1L
    if (n_d - ml < 5L) next
    idx    <- (ml + 1L):n_d
    y_dep  <- dy[idx]
    y_lag1 <- y[idx]    # y_{t-1}

    Xlist <- list(intercept = rep(1, length(idx)))
    if (adf_case == "ct") Xlist$trend <- seq_along(idx)
    Xlist$y_lag1 <- y_lag1

    if (p > 0L) {
      for (j in seq_len(p)) {
        Xlist[[paste0("Ldy", j)]] <- dy[idx - j]
      }
    }
    X_mat <- do.call(cbind, Xlist)
    if (adf_case == "n") {
      X_mat <- X_mat[, -1L, drop = FALSE]  # remove intercept
    }
    fit <- tryCatch(lm(y_dep ~ X_mat - 1), error = function(e) NULL)
    if (is.null(fit)) next

    ll_v <- as.numeric(logLik(fit))
    k_v  <- length(coef(fit))
    n_u  <- length(y_dep)
    ic_v <- if (ic == "aic") -2 * ll_v + 2 * k_v else
      -2 * ll_v + k_v * log(n_u)
    if (ic_v < best_ic) {
      best_ic  <- ic_v
      opt_lag  <- p
      b        <- coef(fit)
      V        <- vcov(fit)
      ycol     <- which(grepl("y_lag1", names(b)))
      if (length(ycol) == 1L) {
        opt_stat <- b[[ycol]] / sqrt(V[ycol, ycol])
      }
    }
  }
  list(stat = opt_stat, opt_lag = opt_lag)
}


#' MacKinnon (1994) approximate ADF critical values
#' @keywords internal
.rardl_adf_cvs <- function(adf_case, level) {
  ## Critical values from MacKinnon (1994) for n -> infinity
  cv_c  <- c("1%" = -3.43, "5%" = -2.86, "10%" = -2.57)
  cv_ct <- c("1%" = -3.96, "5%" = -3.41, "10%" = -3.13)
  cv_n  <- c("1%" = -2.58, "5%" = -1.95, "10%" = -1.62)
  switch(adf_case,
         c  = cv_c,
         ct = cv_ct,
         n  = cv_n)
}


#' Recursive Granger causality test
#' @keywords internal
.rardl_rgranger <- function(y, x, init_obs, max_lag, ic, level, nsim,
                             depvar, indepvar) {
  n <- length(y)

  results <- data.frame(
    obs     = seq(init_obs, n),
    F_granger = NA_real_,
    p_granger = NA_real_,
    opt_lag   = NA_integer_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(results$obs)) {
    t_end  <- results$obs[i]
    yr     <- y[1L:t_end]
    xr     <- x[1L:t_end]
    gc_res <- .rardl_granger_one(yr, xr, max_lag = max_lag, ic = ic)
    results$F_granger[i] <- gc_res$F
    results$p_granger[i] <- gc_res$p
    results$opt_lag[i]   <- gc_res$opt_lag
  }

  structure(
    list(type      = "rgranger",
         results   = results,
         depvar    = depvar,
         indepvars = indepvar,
         nobs      = n,
         init_obs  = init_obs),
    class = "rardl"
  )
}


#' Granger causality test for one sample
#' @keywords internal
.rardl_granger_one <- function(y, x, max_lag, ic) {
  n   <- length(y)
  best_ic  <- Inf
  opt_lag  <- 1L
  opt_F    <- NA_real_
  opt_p    <- NA_real_

  for (p in 1L:min(max_lag, floor(n / 4L))) {
    ml <- p
    if (n - ml < 5L) next
    idx <- (ml + 1L):n

    ## Unrestricted: y ~ lags(y) + lags(x)
    Xu <- cbind(1, do.call(cbind, lapply(seq_len(p), function(j) y[idx - j])),
                do.call(cbind, lapply(seq_len(p), function(j) x[idx - j])))
    ## Restricted: y ~ lags(y)
    Xr <- Xu[, 1L:(p + 1L), drop = FALSE]
    y_dep <- y[idx]

    fit_u <- tryCatch(lm(y_dep ~ Xu - 1), error = function(e) NULL)
    fit_r <- tryCatch(lm(y_dep ~ Xr - 1), error = function(e) NULL)
    if (is.null(fit_u) || is.null(fit_r)) next

    n_u  <- length(y_dep)
    k_u  <- ncol(Xu)
    ll_u <- as.numeric(logLik(fit_u))
    ic_v <- if (ic == "aic") -2 * ll_u + 2 * k_u else
      -2 * ll_u + k_u * log(n_u)

    if (ic_v < best_ic) {
      best_ic <- ic_v
      opt_lag  <- p
      rss_u <- sum(residuals(fit_u)^2)
      rss_r <- sum(residuals(fit_r)^2)
      df_r  <- n_u - k_u
      F_val <- ((rss_r - rss_u) / p) / (rss_u / df_r)
      opt_F <- F_val
      opt_p <- pf(F_val, df1 = p, df2 = df_r, lower.tail = FALSE)
    }
  }
  list(F = opt_F, p = opt_p, opt_lag = opt_lag)
}


#' Monte Carlo simulation of critical values
#' @keywords internal
.rardl_simulate <- function(tobs, nregs, nsim, max_lag, ic, case, level) {
  alpha_map <- c("1" = 0.01, "5" = 0.05, "10" = 0.10)
  alpha     <- alpha_map[as.character(level)]
  if (is.na(alpha)) alpha <- 0.05

  cv_list <- list()
  for (nn in tobs) {
    F_vec <- numeric(nsim)
    for (s in seq_len(nsim)) {
      y   <- cumsum(rnorm(nn))
      xmat <- matrix(cumsum(rnorm(nn * nregs)), nrow = nn, ncol = nregs)
      st   <- tryCatch(
        .rardl_ecm_stats(y, xmat, max_lag = max_lag, ic = ic),
        error = function(e) list(F_pss = NA_real_)
      )
      F_vec[s] <- if (!is.null(st$F_pss)) st$F_pss else NA_real_
    }
    q_val <- quantile(F_vec, probs = 1 - alpha, na.rm = TRUE)
    cv_list[[as.character(nn)]] <- q_val
  }

  results_df <- data.frame(
    n_obs  = tobs,
    cv_val = unlist(cv_list),
    level  = level,
    nregs  = nregs,
    stringsAsFactors = FALSE
  )

  structure(
    list(type            = "simulate",
         results         = results_df,
         critical_values = setNames(unlist(cv_list),
                                    paste0("n=", tobs)),
         nregs           = nregs,
         nsim            = nsim),
    class = "rardl"
  )
}


#' Simulate critical values for a specific window size
#' @keywords internal
.rardl_sim_cvs <- function(n_obs, n_x, nsim, max_lag, ic, level) {
  ## Approximate using reduced nsim for within-function calls
  nsim_fast <- min(nsim, 500L)
  alpha_map <- c("1" = 0.01, "5" = 0.05, "10" = 0.10)
  alpha     <- alpha_map[as.character(level)]
  if (is.na(alpha)) alpha <- 0.05

  F_vec <- numeric(nsim_fast)
  t_vec <- numeric(nsim_fast)
  for (s in seq_len(nsim_fast)) {
    y    <- cumsum(rnorm(n_obs))
    xmat <- matrix(replicate(n_x, cumsum(rnorm(n_obs))), nrow = n_obs)
    st   <- tryCatch(
      .rardl_ecm_stats(y, xmat, max_lag = max_lag, ic = ic),
      error = function(e) list(F_pss = NA_real_, t_pss = NA_real_)
    )
    F_vec[s] <- if (!is.null(st$F_pss)) st$F_pss else NA_real_
    t_vec[s] <- if (!is.null(st$t_pss)) st$t_pss else NA_real_
  }
  c(
    F_cv = quantile(F_vec, probs = 1 - alpha, na.rm = TRUE),
    t_cv = quantile(t_vec, probs = alpha, na.rm = TRUE)
  )
}
