#' Print method for rardl objects
#'
#' @param x An object of class \code{"rardl"}.
#' @param digits Integer. Significant digits.
#' @param ... Further arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.rardl <- function(x, digits = 4L, ...) {
  message(strrep("-", 62))
  message("Rolling/Recursive ARDL Analysis (rardl)")
  message(sprintf("Type: %s", x$type))
  if (!is.null(x$depvar)) {
    message(sprintf("Dependent variable : %s", x$depvar))
  }
  if (!is.null(x$indepvars)) {
    message(sprintf("Independent vars   : %s",
                    paste(x$indepvars, collapse = ", ")))
  }
  if (!is.null(x$nobs)) {
    message(sprintf("Total observations : %d", x$nobs))
  }
  if (!is.null(x$window_size)) {
    message(sprintf("Window size        : %d", x$window_size))
  }
  if (!is.null(x$init_obs)) {
    message(sprintf("Initial obs        : %d", x$init_obs))
  }
  message(strrep("-", 62))

  if (!is.null(x$critical_values)) {
    message("Critical values:")
    cv <- x$critical_values
    for (nm in names(cv)) {
      message(sprintf("  %s : %.*f", nm, digits, cv[[nm]]))
    }
    message(strrep("-", 62))
  }

  if (!is.null(x$results) && nrow(x$results) > 0L) {
    if (x$type %in% c("rolling", "recursive")) {
      message("Test statistics (first 10 windows):")
      show_n <- min(10L, nrow(x$results))
      res_show <- x$results[seq_len(show_n), , drop = FALSE]
      for (i in seq_len(nrow(res_show))) {
        row <- res_show[i, ]
        obs_col <- if ("window_end" %in% names(row)) row$window_end else row$obs
        message(sprintf("  [%4d]  F_pss=%.*f  t_pss=%.*f  F_ind=%.*f",
                        obs_col,
                        digits, row$F_pss,
                        digits, row$t_pss,
                        digits, row$F_ind))
      }
      if (nrow(x$results) > 10L) {
        message(sprintf("  ... (%d more rows)", nrow(x$results) - 10L))
      }
    } else if (x$type == "radf") {
      message("Recursive ADF statistics (first 10 obs):")
      show_n <- min(10L, nrow(x$results))
      for (i in seq_len(show_n)) {
        row <- x$results[i, ]
        message(sprintf("  [%4d]  ADF=%.*f  (lag=%d)",
                        row$obs, digits, row$adf_stat, row$opt_lag))
      }
    } else if (x$type == "rgranger") {
      message("Recursive Granger statistics (first 10 obs):")
      show_n <- min(10L, nrow(x$results))
      for (i in seq_len(show_n)) {
        row <- x$results[i, ]
        message(sprintf("  [%4d]  F=%.*f  p=%.*f  (lag=%d)",
                        row$obs, digits, row$F_granger,
                        digits, row$p_granger, row$opt_lag))
      }
    } else if (x$type == "simulate") {
      message("Simulated critical values:")
      for (i in seq_len(nrow(x$results))) {
        row <- x$results[i, ]
        message(sprintf("  n=%d  CV(%d%%) = %.*f",
                        row$n_obs, row$level, digits, row$cv_val))
      }
    }
  }
  message(strrep("-", 62))
  invisible(x)
}


#' Summary method for rardl objects
#'
#' @param object An object of class \code{"rardl"}.
#' @param ... Further arguments passed to \code{print.rardl}.
#' @return Invisibly returns \code{object}.
#' @export
summary.rardl <- function(object, ...) {
  print(object, ...)
}
