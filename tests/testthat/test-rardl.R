test_that("rardl recursive returns correct class", {
  set.seed(1)
  n  <- 60
  df <- data.frame(y = cumsum(rnorm(n)), x1 = cumsum(rnorm(n)))
  res <- rardl(df, type = "recursive", depvar = "y",
               indepvars = "x1", init_obs = 25, max_lag = 1,
               nsim = 50)
  expect_s3_class(res, "rardl")
  expect_equal(res$type, "recursive")
})

test_that("rardl rolling returns correct class", {
  set.seed(2)
  n  <- 60
  df <- data.frame(y = cumsum(rnorm(n)), x1 = cumsum(rnorm(n)))
  res <- rardl(df, type = "rolling", depvar = "y",
               indepvars = "x1", window_size = 25, max_lag = 1,
               nsim = 50)
  expect_s3_class(res, "rardl")
  expect_equal(res$type, "rolling")
  expect_equal(res$window_size, 25L)
})

test_that("rardl radf returns correct class", {
  set.seed(3)
  n  <- 60
  df <- data.frame(y = cumsum(rnorm(n)))
  res <- rardl(df, type = "radf", depvar = "y",
               init_obs = 25, max_lag = 2)
  expect_s3_class(res, "rardl")
  expect_equal(res$type, "radf")
})

test_that("rardl rgranger requires exactly one indepvar", {
  set.seed(4)
  n  <- 60
  df <- data.frame(y = cumsum(rnorm(n)), x1 = cumsum(rnorm(n)),
                   x2 = cumsum(rnorm(n)))
  expect_error(
    rardl(df, type = "rgranger", depvar = "y",
          indepvars = c("x1", "x2"), init_obs = 25),
    "exactly one independent variable"
  )
})

test_that("rardl simulate returns correct class", {
  set.seed(5)
  res <- rardl(type = "simulate", tobs = c(30L, 60L),
               nregs = 1L, nsim = 50L)
  expect_s3_class(res, "rardl")
  expect_equal(res$type, "simulate")
})

test_that("rardl validates type", {
  set.seed(1)
  df <- data.frame(y = rnorm(40))
  expect_error(rardl(df, type = "unknown"), "'type' must be one of")
})

test_that("print.rardl uses message not cat", {
  set.seed(1)
  n  <- 50
  df <- data.frame(y = cumsum(rnorm(n)), x1 = cumsum(rnorm(n)))
  res <- rardl(df, type = "recursive", depvar = "y",
               indepvars = "x1", init_obs = 20, max_lag = 1,
               nsim = 50)
  expect_message(print(res))
  out <- capture.output(print(res))
  expect_equal(length(out), 0L)
})
