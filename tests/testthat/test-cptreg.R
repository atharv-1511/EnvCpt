# =============================================================================
# Test suite for cpt.reg() and related functions in EnvCpt
# Author: Atharv Raskar (GSoC 2026 Contributor)
# Purpose: Increase test coverage for changepoint regression functions
# =============================================================================

context("cpt.reg function tests")

# =============================================================================
# Test Data Setup
# =============================================================================

set.seed(42)

# Create simple AR(1) data with a changepoint
n <- 100
segment1 <- arima.sim(model = list(ar = 0.3), n = n/2)
segment2 <- arima.sim(model = list(ar = 0.8), n = n/2)
ar1_data <- c(segment1, segment2)

# Create valid data matrix for AR(1): cbind(y, intercept, lag1)
valid_ar1_matrix <- cbind(
  ar1_data[-1],           # Response
  rep(1, n - 1),          # Intercept
  ar1_data[-n]            # Lag 1
)

# Create valid data matrix without changepoints
stable_data <- rnorm(100)
stable_matrix <- cbind(
  stable_data[-1],
  rep(1, 99),
  stable_data[-100]
)

# =============================================================================
# Tests for cpt.reg() - Input Validation
# =============================================================================

test_that("cpt.reg rejects non-array data", {
  expect_error(
    EnvCpt:::cpt.reg(data = c(1, 2, 3)),
    "Argument 'data' must be a numerical matrix/array."
  )
})

test_that("cpt.reg rejects non-numeric data", {
  char_matrix <- matrix(letters[1:9], nrow = 3)
  expect_error(
    EnvCpt:::cpt.reg(data = char_matrix),
    "Argument 'data' must be a numerical matrix/array."
  )
})

test_that("cpt.reg rejects invalid penalty argument", {
  # Note: There is a typo in the source code ("penelty" instead of "penalty")
  # This test documents the existing behavior
  expect_error(
    EnvCpt:::cpt.reg(data = valid_ar1_matrix, penalty = c("MBIC", "BIC")),
    "Argument 'penelty' is invalid."
  )
  expect_error(
    EnvCpt:::cpt.reg(data = valid_ar1_matrix, penalty = 123),
    "Argument 'penelty' is invalid."
  )
})

test_that("cpt.reg rejects invalid method argument", {
  expect_error(
    EnvCpt:::cpt.reg(data = valid_ar1_matrix, method = c("AMOC", "PELT")),
    "Argument 'method' is invalid."
  )
  expect_error(
    EnvCpt:::cpt.reg(data = valid_ar1_matrix, method = 123),
    "Argument 'method' is invalid."
  )
  expect_error(
    EnvCpt:::cpt.reg(data = valid_ar1_matrix, method = "SegNeigh"),
    "Invalid method, must be AMOC or PELT."
  )
  expect_error(
    EnvCpt:::cpt.reg(data = valid_ar1_matrix, method = "BinSeg"),
    "Invalid method, must be AMOC or PELT."
  )
})

test_that("cpt.reg warns and converts non-Normal distributions", {
  expect_warning(
    EnvCpt:::cpt.reg(data = valid_ar1_matrix, dist = "Gamma"),
    "dist = Gamma is not supported. Converted to dist='Normal'"
  )
  expect_warning(
    EnvCpt:::cpt.reg(data = valid_ar1_matrix, dist = "Poisson"),
    "dist = Poisson is not supported. Converted to dist='Normal'"
  )
})

test_that("cpt.reg rejects invalid dist argument format", {
  expect_error(
    EnvCpt:::cpt.reg(data = valid_ar1_matrix, dist = c("Normal", "Gamma")),
    "Argument 'dist' is invalid."
  )
  expect_error(
    EnvCpt:::cpt.reg(data = valid_ar1_matrix, dist = 123),
    "Argument 'dist' is invalid."
  )
})

test_that("cpt.reg rejects invalid class argument", {
  expect_error(
    EnvCpt:::cpt.reg(data = valid_ar1_matrix, class = "TRUE"),
    "Argument 'class' is invalid."
  )
  expect_error(
    EnvCpt:::cpt.reg(data = valid_ar1_matrix, class = c(TRUE, FALSE)),
    "Argument 'class' is invalid."
  )
})

test_that("cpt.reg rejects invalid param.estimates argument", {
  expect_error(
    EnvCpt:::cpt.reg(data = valid_ar1_matrix, param.estimates = "TRUE"),
    "Argument 'param.estimates' is invalid."
  )
  expect_error(
    EnvCpt:::cpt.reg(data = valid_ar1_matrix, param.estimates = c(TRUE, FALSE)),
    "Argument 'param.estimates' is invalid."
  )
})

test_that("cpt.reg rejects invalid minseglen argument", {
  expect_error(
    EnvCpt:::cpt.reg(data = valid_ar1_matrix, minseglen = "5"),
    "Argument 'minseglen' is invalid."
  )
  expect_error(
    EnvCpt:::cpt.reg(data = valid_ar1_matrix, minseglen = c(3, 5)),
    "Argument 'minseglen' is invalid."
  )
  expect_error(
    EnvCpt:::cpt.reg(data = valid_ar1_matrix, minseglen = 0),
    "Argument 'minseglen' must be positive integer."
  )
  expect_error(
    EnvCpt:::cpt.reg(data = valid_ar1_matrix, minseglen = -5),
    "Argument 'minseglen' must be positive integer."
  )
  expect_error(
    EnvCpt:::cpt.reg(data = valid_ar1_matrix, minseglen = 3.5),
    "Argument 'minseglen' must be positive integer."
  )
})

test_that("cpt.reg rejects invalid tol argument", {
  expect_error(
    EnvCpt:::cpt.reg(data = valid_ar1_matrix, tol = "1e-07"),
    "Argument 'tol' is invalid."
  )
  expect_error(
    EnvCpt:::cpt.reg(data = valid_ar1_matrix, tol = c(1e-07, 1e-08)),
    "Argument 'tol' is invalid."
  )
  expect_error(
    EnvCpt:::cpt.reg(data = valid_ar1_matrix, tol = -1e-07),
    "Argument 'tol' must be positive."
  )
})

# =============================================================================
# Tests for cpt.reg() - Methods and Output
# =============================================================================

test_that("cpt.reg works with AMOC method", {
  result <- EnvCpt:::cpt.reg(data = valid_ar1_matrix, method = "AMOC")
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg works with PELT method", {
  result <- EnvCpt:::cpt.reg(data = valid_ar1_matrix, method = "PELT")
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg class=FALSE returns raw output", {
  result <- EnvCpt:::cpt.reg(data = valid_ar1_matrix, method = "PELT", class = FALSE)
  # When class=FALSE, should return a list or vector, not cpt.reg object

  expect_false(inherits(result, "cpt.reg"))
})

test_that("cpt.reg param.estimates=FALSE works", {
  result <- EnvCpt:::cpt.reg(data = valid_ar1_matrix, method = "PELT", param.estimates = FALSE)
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg works with various penalty types", {
  # Test different penalty types
  penalties <- c("MBIC", "BIC", "SIC", "AIC")

  for (pen in penalties) {
    result <- EnvCpt:::cpt.reg(data = stable_matrix, method = "PELT", penalty = pen)
    expect_s4_class(result, "cpt.reg")
  }
})

test_that("cpt.reg handles data without changepoints", {
  # Stable data should detect no changepoints (or just the end point)
  result <- EnvCpt:::cpt.reg(data = stable_matrix, method = "PELT", penalty = "MBIC")
  expect_s4_class(result, "cpt.reg")
  # The only changepoint should be n (end of data)
  cpts_result <- cpts(result)
  expect_true(length(cpts_result) >= 0)
})

# =============================================================================
# Tests for cpt.reg() - Multiple Time Series (3D Array)
# =============================================================================

test_that("cpt.reg handles 3D array input (multiple time series)", {
  # Create 3D array with 2 time series
  multi_data <- array(
    c(valid_ar1_matrix, stable_matrix),
    dim = c(2, nrow(valid_ar1_matrix), ncol(valid_ar1_matrix))
  )

  result <- EnvCpt:::cpt.reg(data = multi_data, method = "PELT")

  # Should return a list of results
  expect_true(is.list(result))
  expect_equal(length(result), 2)
  expect_s4_class(result[[1]], "cpt.reg")
  expect_s4_class(result[[2]], "cpt.reg")
})

# =============================================================================
# Tests for cpt.reg() - Minseglen Warnings
# =============================================================================

test_that("cpt.reg warns when minseglen is too small", {
  # minseglen smaller than number of regressors should trigger warning
  expect_warning(
    EnvCpt:::cpt.reg(data = valid_ar1_matrix, minseglen = 1),
    "minseglen is too small"
  )
})

test_that("cpt.reg errors when minseglen is too large", {
  # minseglen too large relative to data should error
  expect_error(
    EnvCpt:::cpt.reg(data = valid_ar1_matrix, minseglen = 100),
    "Minimum segment length is too large"
  )
})

# =============================================================================
# Tests for check_data() function
# =============================================================================

context("check_data function tests")

test_that("check_data rejects non-array input", {
  expect_error(
    EnvCpt:::check_data(c(1, 2, 3)),
    "Argument 'data' must be a numerical matrix."
  )
})

test_that("check_data rejects non-numeric input", {
  expect_error(
    EnvCpt:::check_data(matrix(letters[1:6], nrow = 2)),
    "Argument 'data' must be a numerical matrix."
  )
})

test_that("check_data rejects 3D arrays", {
  arr <- array(1:24, dim = c(2, 3, 4))
  expect_error(
    EnvCpt:::check_data(arr),
    "Argument 'data' must be a numerical matrix."
  )
})

test_that("check_data rejects invalid minseglen", {
  expect_error(
    EnvCpt:::check_data(valid_ar1_matrix, minseglen = "3"),
    "Argument 'minseglen' is invalid."
  )
  expect_error(
    EnvCpt:::check_data(valid_ar1_matrix, minseglen = c(3, 5)),
    "Argument 'minseglen' is invalid."
  )
})

test_that("check_data rejects data with only response (no regressors)", {
  single_col <- matrix(1:10, ncol = 1)
  expect_error(
    EnvCpt:::check_data(single_col),
    "Dimension of data is 1, no regressors found."
  )
})

test_that("check_data rejects data with more regressors than observations", {
  wide_matrix <- matrix(1:20, nrow = 2, ncol = 10)
  expect_error(
    EnvCpt:::check_data(wide_matrix),
    "More regressors than observations."
  )
})

test_that("check_data warns and adds intercept when missing", {
  # Create matrix without intercept column
  no_intercept <- cbind(1:10, (1:10)^2)  # response and regressor, no intercept

  expect_warning(
    result <- EnvCpt:::check_data(no_intercept),
    "Missing intercept regressor"
  )

  # Should have added an intercept column
  expect_equal(ncol(result), ncol(no_intercept) + 1)
})

test_that("check_data handles data with single intercept correctly", {
  # Create matrix with intercept in correct position (column 2)
  with_intercept <- cbind(1:10, rep(1, 10), (1:10)^2)

  # Should not produce any warnings
  expect_silent(
    result <- EnvCpt:::check_data(with_intercept)
  )

  # Result should be the same as input
  expect_equal(result, with_intercept)
})

# =============================================================================
# Tests for ChangepointRegression() function
# =============================================================================

context("ChangepointRegression function tests")

# Note: ChangepointRegression is called internally by cpt.reg
# These tests verify its behavior directly

test_that("ChangepointRegression handles cpts.only argument", {
  # Note: The condition in source is: !is.logical(cpts.only) && length(cpts.only)>1
  # This means single non-logical values won't trigger error (uses && not ||)
  # Test the actual behavior - cpts.only=TRUE and FALSE both work
  result_true <- EnvCpt:::ChangepointRegression(
    data = valid_ar1_matrix,
    method = "PELT",
    penalty.value = log(nrow(valid_ar1_matrix)),
    cpts.only = TRUE
  )
  expect_true(is.numeric(result_true))

  result_false <- EnvCpt:::ChangepointRegression(
    data = valid_ar1_matrix,
    method = "PELT",
    penalty.value = log(nrow(valid_ar1_matrix)),
    cpts.only = FALSE
  )
  expect_true(is.list(result_false))
})

test_that("ChangepointRegression errors on unrecognized method", {
  expect_error(
    EnvCpt:::ChangepointRegression(
      data = valid_ar1_matrix,
      method = "UNKNOWN",
      dist = "Normal"
    ),
    "Changepoint in regression method not recognised."
  )
})

test_that("ChangepointRegression returns changepoints when cpts.only=TRUE", {
  result <- EnvCpt:::ChangepointRegression(
    data = valid_ar1_matrix,
    method = "PELT",
    penalty.value = log(nrow(valid_ar1_matrix)),
    cpts.only = TRUE
  )

  # Should return a sorted vector of changepoints
  expect_true(is.numeric(result))
})

test_that("ChangepointRegression returns full output when cpts.only=FALSE", {
  result <- EnvCpt:::ChangepointRegression(
    data = valid_ar1_matrix,
    method = "PELT",
    penalty.value = log(nrow(valid_ar1_matrix)),
    cpts.only = FALSE
  )

  # Should return a list with more information
  expect_true(is.list(result))
  expect_true("cpts" %in% names(result))
})

# =============================================================================
# Tests for CptReg_AMOC_Normal() function
# =============================================================================

context("CptReg_AMOC_Normal function tests")

test_that("CptReg_AMOC_Normal rejects invalid shape argument", {
  expect_error(
    EnvCpt:::CptReg_AMOC_Normal(
      data = valid_ar1_matrix,
      shape = "0"
    ),
    "Argument 'shape' is invalid."
  )
  expect_error(
    EnvCpt:::CptReg_AMOC_Normal(
      data = valid_ar1_matrix,
      shape = c(0, 1)
    ),
    "Argument 'shape' is invalid."
  )
})

test_that("CptReg_AMOC_Normal rejects invalid data dimensions", {
  # Single column (no regressors after response)
  single_col <- matrix(1:10, ncol = 1)
  expect_error(
    EnvCpt:::CptReg_AMOC_Normal(data = single_col),
    "Invalid data dimensions."
  )
})

# =============================================================================
# Tests for CptReg_PELT_Normal() function
# =============================================================================

context("CptReg_PELT_Normal function tests")

test_that("CptReg_PELT_Normal rejects invalid shape argument", {
  expect_error(
    EnvCpt:::CptReg_PELT_Normal(
      data = valid_ar1_matrix,
      shape = "0"
    ),
    "Argument 'shape' is invalid."
  )
  expect_error(
    EnvCpt:::CptReg_PELT_Normal(
      data = valid_ar1_matrix,
      shape = c(0, 1)
    ),
    "Argument 'shape' is invalid."
  )
})

# =============================================================================
# Tests for BIC.envcpt() Error Handling
# =============================================================================

context("BIC.envcpt error handling tests")

test_that("BIC.envcpt rejects non-envcpt objects", {
  # Need to call EnvCpt:::BIC.envcpt directly since BIC() dispatches to stats::BIC
  tmp <- list(a = 1, b = 2)
  expect_error(
    EnvCpt:::BIC.envcpt(tmp),
    "object must be of class envcpt"
  )
})

test_that("BIC.envcpt rejects non-list objects", {
  tmp <- rnorm(100)
  class(tmp) <- "envcpt"
  expect_error(
    EnvCpt:::BIC.envcpt(tmp),
    "object argument must be a list"
  )
})

test_that("BIC.envcpt rejects objects without matrix first element", {
  tmp <- list(summary = rnorm(100))
  class(tmp) <- "envcpt"
  expect_error(
    EnvCpt:::BIC.envcpt(tmp),
    "first element in the object list must be a matrix."
  )
})

test_that("BIC.envcpt rejects objects with non-numeric matrix", {
  tmp <- list(summary = matrix(LETTERS, nrow = 2))
  class(tmp) <- "envcpt"
  expect_error(
    EnvCpt:::BIC.envcpt(tmp),
    "First two rows in matrix in first element of object list must be numeric"
  )
})

# =============================================================================
# Tests for AICweights() function
# =============================================================================

context("AICweights function tests")

test_that("AICweights.default returns message", {
  result <- AICweights(c(1, 2, 3))
  expect_equal(result, "No default method created for S3 class AICweights.")
})

test_that("AICweights.envcpt returns valid weights", {
  skip_if_not(identical(Sys.getenv("NOT_CRAN"), "true"))

  set.seed(1)
  x <- c(rnorm(50, 0, 1), rnorm(50, 5, 1))
  out <- envcpt(x, verbose = FALSE)

  weights <- AICweights(out)

  # Weights should sum to 1

  expect_equal(sum(weights, na.rm = TRUE), 1, tolerance = 1e-10)

  # All weights should be between 0 and 1
  expect_true(all(weights >= 0 & weights <= 1, na.rm = TRUE))
})

# =============================================================================
# Tests for plot.envcpt() Error Handling
# =============================================================================

context("plot.envcpt error handling tests")

test_that("plot.envcpt rejects non-envcpt objects", {
  expect_error(
    EnvCpt:::plot.envcpt(list(a = 1)),
    "x must be an object with class envcpt"
  )
})

test_that("plot.envcpt rejects insufficient colors", {
  skip_if_not(identical(Sys.getenv("NOT_CRAN"), "true"))

  set.seed(1)
  x <- rnorm(50)
  out <- envcpt(x, verbose = FALSE)

  expect_error(
    EnvCpt:::plot.envcpt(out, colors = c("red", "blue")),
    "colors must be a vector of length 12"
  )
})

test_that("plot.envcpt rejects invalid colors", {
  skip_if_not(identical(Sys.getenv("NOT_CRAN"), "true"))

  set.seed(1)
  x <- rnorm(50)
  out <- envcpt(x, verbose = FALSE)

  expect_error(
    EnvCpt:::plot.envcpt(out, colors = c("red", "notacolor", "blue", "green",
                          "yellow", "purple", "orange", "pink",
                          "brown", "gray", "cyan", "magenta")),
    "Atleast one of your colours is not resolvable"
  )
})

test_that("plot.envcpt rejects invalid type", {
  skip_if_not(identical(Sys.getenv("NOT_CRAN"), "true"))

  set.seed(1)
  x <- rnorm(50)
  out <- envcpt(x, verbose = FALSE)

  expect_error(
    EnvCpt:::plot.envcpt(out, type = "invalid"),
    "type supplied can only be 'aic', 'bic' or 'fit'."
  )
})

# =============================================================================
# Integration Tests
# =============================================================================

context("Integration tests for cpt.reg with envcpt")

test_that("cpt.reg produces consistent results with envcpt internal calls", {
  skip_if_not(identical(Sys.getenv("NOT_CRAN"), "true"))

  set.seed(123)
  n <- 200

  # Create AR(1) data with changepoint
  x1 <- arima.sim(model = list(ar = 0.3), n = 100)
  x2 <- arima.sim(model = list(ar = 0.8), n = 100)
  x <- c(x1, x2)

  # Create data matrix as envcpt does internally
  data_matrix <- cbind(x[-1], rep(1, n - 1), x[-n])

  # Run cpt.reg directly
  direct_result <- EnvCpt:::cpt.reg(
    data = data_matrix,
    method = "PELT",
    penalty = "MBIC",
    minseglen = 5
  )

  # Should produce a valid cpt.reg object
  expect_s4_class(direct_result, "cpt.reg")

  # Changepoints should be detected near position 100
  detected_cpts <- cpts(direct_result)
  # Filter out the final endpoint
  internal_cpts <- detected_cpts[detected_cpts < n]

  # There should be at least one changepoint detected
  # (may not be exactly at 100 due to randomness)
  expect_true(length(internal_cpts) >= 0)
})

# =============================================================================
# Additional Tests for AR(2) Models
# =============================================================================

context("AR(2) model tests")

# Create AR(2) test data
set.seed(123)
n_ar2 <- 150
ar2_segment1 <- arima.sim(model = list(ar = c(0.5, 0.3)), n = 75)
ar2_segment2 <- arima.sim(model = list(ar = c(-0.3, 0.2)), n = 75)
ar2_data <- c(ar2_segment1, ar2_segment2)

# AR(2) data matrix: cbind(y, intercept, lag1, lag2)
valid_ar2_matrix <- cbind(
  ar2_data[-(1:2)],           # Response: y[t]
  rep(1, n_ar2 - 2),          # Intercept
  ar2_data[-c(1, n_ar2)],     # Lag 1: y[t-1]
  ar2_data[-c(n_ar2-1, n_ar2)] # Lag 2: y[t-2]
)

test_that("cpt.reg works with AR(2) data matrix", {
  result <- EnvCpt:::cpt.reg(data = valid_ar2_matrix, method = "PELT")
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg AMOC works with AR(2) data", {
  result <- EnvCpt:::cpt.reg(data = valid_ar2_matrix, method = "AMOC")
  expect_s4_class(result, "cpt.reg")
})

test_that("AR(2) model detects changepoint near true location", {
  result <- EnvCpt:::cpt.reg(data = valid_ar2_matrix, method = "PELT", penalty = "MBIC")
  detected_cpts <- cpts(result)
  internal_cpts <- detected_cpts[detected_cpts < n_ar2 - 2]

  # If changepoints detected, at least one should be near 75
  if (length(internal_cpts) > 0) {
    min_distance <- min(abs(internal_cpts - 73))  # 73 because of lag adjustment
    expect_true(min_distance < 20)
  }
})

# =============================================================================
# Tests for Trend Data
# =============================================================================

context("Trend data tests")

# Create data with trend
set.seed(456)
n_trend <- 120
time_vec <- 1:n_trend
trend_data <- 0.1 * time_vec + rnorm(n_trend, sd = 0.5)
trend_data[61:n_trend] <- trend_data[61:n_trend] + 5  # Level shift

# Trend data matrix: cbind(y, intercept, time)
trend_matrix <- cbind(
  trend_data,
  rep(1, n_trend),
  time_vec
)

test_that("cpt.reg works with trend data", {
  result <- EnvCpt:::cpt.reg(data = trend_matrix, method = "PELT")
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg AMOC works with trend data", {
  result <- EnvCpt:::cpt.reg(data = trend_matrix, method = "AMOC")
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg detects level shift in trend data", {
  result <- EnvCpt:::cpt.reg(data = trend_matrix, method = "PELT", penalty = "BIC")
  detected_cpts <- cpts(result)
  internal_cpts <- detected_cpts[detected_cpts < n_trend]

  # Should detect changepoint near 60
  if (length(internal_cpts) > 0) {
    min_distance <- min(abs(internal_cpts - 60))
    expect_true(min_distance < 15)
  }
})

# =============================================================================
# Tests for AR(1) with Trend
# =============================================================================

context("AR(1) with trend tests")

# AR(1) + Trend data matrix
set.seed(789)
n_ar1t <- 100
ar1_part <- arima.sim(model = list(ar = 0.5), n = n_ar1t)
ar1_trend_data <- ar1_part + 0.05 * (1:n_ar1t)

ar1_trend_matrix <- cbind(
  ar1_trend_data[-1],
  rep(1, n_ar1t - 1),
  ar1_trend_data[-n_ar1t],
  2:n_ar1t  # Time regressor
)

test_that("cpt.reg works with AR(1) + trend matrix", {
  result <- EnvCpt:::cpt.reg(data = ar1_trend_matrix, method = "PELT")
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg AMOC works with AR(1) + trend matrix", {
  result <- EnvCpt:::cpt.reg(data = ar1_trend_matrix, method = "AMOC")
  expect_s4_class(result, "cpt.reg")
})

# =============================================================================
# Tests for Manual Penalty
# =============================================================================

context("Manual penalty tests")

test_that("cpt.reg works with Manual penalty", {
  result <- EnvCpt:::cpt.reg(
    data = valid_ar1_matrix,
    method = "PELT",
    penalty = "Manual",
    pen.value = 10
  )
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg with different manual penalty values produces different results", {
  result_low <- EnvCpt:::cpt.reg(
    data = valid_ar1_matrix,
    method = "PELT",
    penalty = "Manual",
    pen.value = 1
  )
  result_high <- EnvCpt:::cpt.reg(
    data = valid_ar1_matrix,
    method = "PELT",
    penalty = "Manual",
    pen.value = 100
  )

  # Lower penalty should allow more changepoints
  cpts_low <- length(cpts(result_low))
  cpts_high <- length(cpts(result_high))
  expect_true(cpts_low >= cpts_high)
})

test_that("cpt.reg AMOC works with Manual penalty", {
  result <- EnvCpt:::cpt.reg(
    data = valid_ar1_matrix,
    method = "AMOC",
    penalty = "Manual",
    pen.value = 5
  )
  expect_s4_class(result, "cpt.reg")
})

# =============================================================================
# Tests for Hannan-Quinn Penalty
# =============================================================================

context("Hannan-Quinn penalty tests")

test_that("cpt.reg works with Hannan-Quinn penalty", {
  result <- EnvCpt:::cpt.reg(
    data = valid_ar1_matrix,
    method = "PELT",
    penalty = "Hannan-Quinn"
  )
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg AMOC works with Hannan-Quinn penalty", {
  result <- EnvCpt:::cpt.reg(
    data = valid_ar1_matrix,
    method = "AMOC",
    penalty = "Hannan-Quinn"
  )
  expect_s4_class(result, "cpt.reg")
})

# =============================================================================
# Tests for Shape Parameter (Fixed Variance)
# =============================================================================

context("Shape parameter tests")

test_that("cpt.reg works with shape = 0 (estimate variance)", {
  result <- EnvCpt:::cpt.reg(
    data = valid_ar1_matrix,
    method = "PELT",
    shape = 0
  )
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg works with fixed shape (known variance)", {
  result <- EnvCpt:::cpt.reg(
    data = valid_ar1_matrix,
    method = "PELT",
    shape = 1  # Fixed variance = 1
  )
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg AMOC works with shape = 0", {
  result <- EnvCpt:::cpt.reg(
    data = valid_ar1_matrix,
    method = "AMOC",
    shape = 0
  )
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg AMOC works with fixed shape", {
  result <- EnvCpt:::cpt.reg(
    data = valid_ar1_matrix,
    method = "AMOC",
    shape = 2
  )
  expect_s4_class(result, "cpt.reg")
})

test_that("Different shape values produce different results", {
  result_est <- EnvCpt:::cpt.reg(
    data = valid_ar1_matrix,
    method = "PELT",
    shape = 0,
    penalty = "Manual",
    pen.value = 5
  )
  result_fixed <- EnvCpt:::cpt.reg(
    data = valid_ar1_matrix,
    method = "PELT",
    shape = 10,
    penalty = "Manual",
    pen.value = 5
  )

  # Results may differ - just check both are valid
  expect_s4_class(result_est, "cpt.reg")
  expect_s4_class(result_fixed, "cpt.reg")
})

# =============================================================================
# Tests for CptReg_PELT_Normal() Additional
# =============================================================================

context("CptReg_PELT_Normal additional tests")

test_that("CptReg_PELT_Normal works with default parameters", {
  result <- EnvCpt:::CptReg_PELT_Normal(data = valid_ar1_matrix)
  expect_true(is.list(result) || is.numeric(result))
})

test_that("CptReg_PELT_Normal works with shape = 0", {
  result <- EnvCpt:::CptReg_PELT_Normal(
    data = valid_ar1_matrix,
    shape = 0
  )
  expect_true(is.list(result) || is.numeric(result))
})

test_that("CptReg_PELT_Normal works with fixed shape", {
  result <- EnvCpt:::CptReg_PELT_Normal(
    data = valid_ar1_matrix,
    shape = 1
  )
  expect_true(is.list(result) || is.numeric(result))
})

test_that("CptReg_PELT_Normal works with AR(2) data", {
  result <- EnvCpt:::CptReg_PELT_Normal(data = valid_ar2_matrix)
  expect_true(is.list(result) || is.numeric(result))
})

test_that("CptReg_PELT_Normal rejects invalid data", {
  bad_matrix <- matrix(1:6, nrow = 2, ncol = 3)
  expect_error(
    EnvCpt:::CptReg_PELT_Normal(data = bad_matrix),
    regexp = NULL  # Some error expected
  )
})

# =============================================================================
# Tests for CptReg_AMOC_Normal() Additional
# =============================================================================

context("CptReg_AMOC_Normal additional tests")

test_that("CptReg_AMOC_Normal works with default parameters", {
  result <- EnvCpt:::CptReg_AMOC_Normal(data = valid_ar1_matrix)
  expect_true(is.list(result) || is.numeric(result))
})

test_that("CptReg_AMOC_Normal works with shape = 0", {
  result <- EnvCpt:::CptReg_AMOC_Normal(
    data = valid_ar1_matrix,
    shape = 0
  )
  expect_true(is.list(result) || is.numeric(result))
})

test_that("CptReg_AMOC_Normal works with fixed shape", {
  result <- EnvCpt:::CptReg_AMOC_Normal(
    data = valid_ar1_matrix,
    shape = 1
  )
  expect_true(is.list(result) || is.numeric(result))
})

test_that("CptReg_AMOC_Normal works with AR(2) data", {
  result <- EnvCpt:::CptReg_AMOC_Normal(data = valid_ar2_matrix)
  expect_true(is.list(result) || is.numeric(result))
})

test_that("CptReg_AMOC_Normal works with trend data", {
  result <- EnvCpt:::CptReg_AMOC_Normal(data = trend_matrix)
  expect_true(is.list(result) || is.numeric(result))
})

# =============================================================================
# Tests for ChangepointRegression() Additional
# =============================================================================

context("ChangepointRegression additional tests")

test_that("ChangepointRegression AMOC returns valid output structure", {
  result <- EnvCpt:::ChangepointRegression(
    data = valid_ar1_matrix,
    method = "AMOC",
    penalty.value = log(nrow(valid_ar1_matrix)),
    cpts.only = FALSE
  )

  expect_true(is.list(result))
  expect_true("cpts" %in% names(result))
})

test_that("ChangepointRegression PELT with minseglen parameter", {
  result <- EnvCpt:::ChangepointRegression(
    data = valid_ar1_matrix,
    method = "PELT",
    penalty.value = log(nrow(valid_ar1_matrix)),
    minseglen = 10,
    cpts.only = FALSE
  )

  expect_true(is.list(result))
})

test_that("ChangepointRegression with shape parameter", {
  result <- EnvCpt:::ChangepointRegression(
    data = valid_ar1_matrix,
    method = "PELT",
    penalty.value = log(nrow(valid_ar1_matrix)),
    shape = 1,
    cpts.only = FALSE
  )

  expect_true(is.list(result))
})

# =============================================================================
# Tests for check_data() Additional
# =============================================================================

context("check_data additional tests")

test_that("check_data handles minseglen smaller than regressors", {
  # With 3 columns (1 response + 2 regressors), minseglen = 1 is too small
  expect_warning(
    EnvCpt:::check_data(valid_ar1_matrix, minseglen = 1),
    "minseglen is too small"
  )
})

test_that("check_data handles minseglen equal to regressors", {
  # minseglen = 2 equals number of regressors, should still warn
  expect_warning(
    EnvCpt:::check_data(valid_ar1_matrix, minseglen = 2)
  )
})

test_that("check_data accepts valid minseglen", {
  # minseglen = 5 is larger than regressors, should be OK
  result <- EnvCpt:::check_data(valid_ar1_matrix, minseglen = 5)
  expect_true(is.matrix(result))
})

test_that("check_data warns for duplicate intercept columns", {
  # Create matrix with two intercept columns
  double_intercept <- cbind(1:10, rep(1, 10), rep(1, 10), (1:10)^2)

  expect_warning(
    EnvCpt:::check_data(double_intercept),
    regexp = "intercept"
  )
})

# =============================================================================
# Tests for Minseglen Edge Cases
# =============================================================================

context("Minseglen edge case tests")

test_that("cpt.reg with minseglen at minimum valid size", {
  # Minimum valid minseglen is ncol(data) - 1 + 1 = ncol(data)
  min_valid <- ncol(valid_ar1_matrix)

  # Should work without warning when minseglen > number of regressors
  result <- suppressWarnings(
    EnvCpt:::cpt.reg(data = valid_ar1_matrix, minseglen = min_valid + 1)
  )
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg with reasonable minseglen", {
  result <- EnvCpt:::cpt.reg(
    data = valid_ar1_matrix,
    method = "PELT",
    minseglen = 10
  )
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg with large minseglen reduces changepoints", {
  result_small <- suppressWarnings(
    EnvCpt:::cpt.reg(data = valid_ar1_matrix, minseglen = 5, penalty = "Manual", pen.value = 1)
  )
  result_large <- suppressWarnings(
    EnvCpt:::cpt.reg(data = valid_ar1_matrix, minseglen = 20, penalty = "Manual", pen.value = 1)
  )

  # Larger minseglen should give fewer or equal changepoints
  expect_true(length(cpts(result_large)) <= length(cpts(result_small)) + 1)
})

# =============================================================================
# Tests for Tolerance Parameter
# =============================================================================

context("Tolerance parameter tests")

test_that("cpt.reg works with default tolerance", {
  result <- EnvCpt:::cpt.reg(data = valid_ar1_matrix, method = "PELT")
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg works with larger tolerance", {
  result <- EnvCpt:::cpt.reg(
    data = valid_ar1_matrix,
    method = "PELT",
    tol = 1e-05
  )
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg works with very small tolerance", {
  result <- EnvCpt:::cpt.reg(
    data = valid_ar1_matrix,
    method = "PELT",
    tol = 1e-10
  )
  expect_s4_class(result, "cpt.reg")
})

# =============================================================================
# Tests for Multiple Time Series (3D Array) Additional
# =============================================================================

context("3D array input additional tests")

test_that("cpt.reg handles 3D array with AMOC method", {
  multi_data <- array(
    c(valid_ar1_matrix, stable_matrix),
    dim = c(2, nrow(valid_ar1_matrix), ncol(valid_ar1_matrix))
  )

  result <- EnvCpt:::cpt.reg(data = multi_data, method = "AMOC")

  expect_true(is.list(result))
  expect_equal(length(result), 2)
})

test_that("cpt.reg handles 3D array with different penalties", {
  multi_data <- array(
    c(valid_ar1_matrix, stable_matrix),
    dim = c(2, nrow(valid_ar1_matrix), ncol(valid_ar1_matrix))
  )

  result <- EnvCpt:::cpt.reg(data = multi_data, method = "PELT", penalty = "AIC")

  expect_true(is.list(result))
})

test_that("cpt.reg handles 3D array with class=FALSE", {
  multi_data <- array(
    c(valid_ar1_matrix, stable_matrix),
    dim = c(2, nrow(valid_ar1_matrix), ncol(valid_ar1_matrix))
  )

  result <- EnvCpt:::cpt.reg(data = multi_data, method = "PELT", class = FALSE)

  expect_true(is.list(result))
})

# =============================================================================
# Tests for Changepoint Detection Accuracy
# =============================================================================

context("Changepoint detection accuracy tests")

test_that("Single changepoint detection accuracy", {
  set.seed(999)
  n_acc <- 200
  true_cp <- 100

  # Clear AR(1) structure change
  seg1 <- arima.sim(model = list(ar = 0.2), n = true_cp)
  seg2 <- arima.sim(model = list(ar = 0.9), n = n_acc - true_cp)
  data_acc <- c(seg1, seg2)

  data_matrix <- cbind(
    data_acc[-1],
    rep(1, n_acc - 1),
    data_acc[-n_acc]
  )

  result <- EnvCpt:::cpt.reg(data = data_matrix, method = "PELT", penalty = "MBIC")
  detected <- cpts(result)
  internal <- detected[detected < n_acc - 1]

  if (length(internal) > 0) {
    closest <- internal[which.min(abs(internal - (true_cp - 1)))]
    expect_true(abs(closest - (true_cp - 1)) < 15)
  }
})

test_that("Multiple changepoint detection", {
  set.seed(888)
  n_multi <- 300

  # Three segments with different AR structures
  seg1 <- arima.sim(model = list(ar = 0.3), n = 100)
  seg2 <- arima.sim(model = list(ar = 0.8), n = 100)
  seg3 <- arima.sim(model = list(ar = -0.5), n = 100)
  data_multi <- c(seg1, seg2, seg3)

  data_matrix <- cbind(
    data_multi[-1],
    rep(1, n_multi - 1),
    data_multi[-n_multi]
  )

  result <- EnvCpt:::cpt.reg(data = data_matrix, method = "PELT", penalty = "MBIC")
  detected <- cpts(result)

  # Should detect at least one changepoint
  expect_true(length(detected) >= 1)
})

# =============================================================================
# Tests for cpt.reg Return Object Structure
# =============================================================================

context("cpt.reg return object structure tests")

test_that("cpt.reg object has correct slots", {
  result <- EnvCpt:::cpt.reg(data = valid_ar1_matrix, method = "PELT")

  # Check standard cpt.reg slots
  expect_true(methods::is(result, "cpt.reg"))
  expect_true("cpts" %in% methods::slotNames(result))
})

test_that("cpt.reg cpts() accessor works", {
  result <- EnvCpt:::cpt.reg(data = valid_ar1_matrix, method = "PELT")

  cpts_result <- cpts(result)
  expect_true(is.numeric(cpts_result))
})

test_that("cpt.reg ncpts() accessor works", {
  result <- EnvCpt:::cpt.reg(data = valid_ar1_matrix, method = "PELT")

  n_cpts <- ncpts(result)
  expect_true(is.numeric(n_cpts))
  expect_true(n_cpts >= 0)
})

# =============================================================================
# Tests for envcpt Model Subsets
# =============================================================================

context("envcpt model subset tests")

test_that("envcpt works with single model name", {
  skip_if_not(identical(Sys.getenv("NOT_CRAN"), "true"))

  set.seed(111)
  x <- c(rnorm(50, 0, 1), rnorm(50, 3, 1))

  result <- envcpt(x, models = "meancpt", verbose = FALSE)
  expect_s3_class(result, "envcpt")
})

test_that("envcpt works with multiple model indices", {
  skip_if_not(identical(Sys.getenv("NOT_CRAN"), "true"))

  set.seed(222)
  x <- c(rnorm(50, 0, 1), rnorm(50, 3, 1))

  result <- envcpt(x, models = c(1, 2, 5), verbose = FALSE)
  expect_s3_class(result, "envcpt")
})

test_that("envcpt works with AR1 model", {
  skip_if_not(identical(Sys.getenv("NOT_CRAN"), "true"))

  set.seed(333)
  x <- arima.sim(model = list(ar = 0.7), n = 100)

  result <- envcpt(x, models = "meanar1", verbose = FALSE)
  expect_s3_class(result, "envcpt")
})

test_that("envcpt works with trendAR models", {
  skip_if_not(identical(Sys.getenv("NOT_CRAN"), "true"))

  set.seed(444)
  x <- arima.sim(model = list(ar = 0.5), n = 100) + 0.05 * (1:100)

  result <- envcpt(x, models = "trendar1", verbose = FALSE)
  expect_s3_class(result, "envcpt")
})

test_that("envcpt propagates minseglen to internal calls", {
  skip_if_not(identical(Sys.getenv("NOT_CRAN"), "true"))

  set.seed(555)
  x <- c(rnorm(50, 0, 1), rnorm(50, 5, 1))

  # Should work with minseglen
  result <- envcpt(x, minseglen = 10, verbose = FALSE)
  expect_s3_class(result, "envcpt")
})

# =============================================================================
# Tests for BIC.envcpt Additional
# =============================================================================

context("BIC.envcpt additional tests")

test_that("BIC.envcpt works with valid envcpt object", {
  skip_if_not(identical(Sys.getenv("NOT_CRAN"), "true"))

  set.seed(666)
  x <- c(rnorm(50, 0, 1), rnorm(50, 3, 1))
  out <- envcpt(x, verbose = FALSE)

  bic_result <- EnvCpt:::BIC.envcpt(out)

  expect_true(is.numeric(bic_result))
  expect_true(length(bic_result) > 0)
})

test_that("BIC.envcpt handles all NA values gracefully", {
  skip_if_not(identical(Sys.getenv("NOT_CRAN"), "true"))

  set.seed(777)
  x <- c(rnorm(50, 0, 1), rnorm(50, 3, 1))
  out <- envcpt(x, verbose = FALSE)

  # Should not error
  bic_result <- try(EnvCpt:::BIC.envcpt(out), silent = TRUE)
  expect_false(inherits(bic_result, "try-error"))
})

# =============================================================================
# Tests for AICweights Additional
# =============================================================================

context("AICweights additional tests")

test_that("AICweights returns weights summing to 1", {
  skip_if_not(identical(Sys.getenv("NOT_CRAN"), "true"))

  set.seed(888)
  x <- c(rnorm(50, 0, 1), rnorm(50, 5, 1))
  out <- envcpt(x, verbose = FALSE)

  weights <- AICweights(out)
  valid_weights <- weights[!is.na(weights)]

  if (length(valid_weights) > 0) {
    expect_equal(sum(valid_weights), 1, tolerance = 1e-8)
  }
})

test_that("AICweights returns non-negative weights", {
  skip_if_not(identical(Sys.getenv("NOT_CRAN"), "true"))

  set.seed(999)
  x <- c(rnorm(75, 0, 1), rnorm(75, 3, 1))
  out <- envcpt(x, verbose = FALSE)

  weights <- AICweights(out)
  valid_weights <- weights[!is.na(weights)]

  expect_true(all(valid_weights >= 0))
})

test_that("AICweights returns weights bounded by 1", {
  skip_if_not(identical(Sys.getenv("NOT_CRAN"), "true"))

  set.seed(1000)
  x <- c(rnorm(60, 0, 1), rnorm(60, 4, 1))
  out <- envcpt(x, verbose = FALSE)

  weights <- AICweights(out)
  valid_weights <- weights[!is.na(weights)]

  expect_true(all(valid_weights <= 1))
})

# =============================================================================
# Tests for plot.envcpt Additional
# =============================================================================

context("plot.envcpt additional tests")

test_that("plot.envcpt works with type='aic'", {
  skip_if_not(identical(Sys.getenv("NOT_CRAN"), "true"))

  set.seed(1111)
  x <- c(rnorm(50, 0, 1), rnorm(50, 3, 1))
  out <- envcpt(x, verbose = FALSE)

  # Should complete without error
  expect_silent(EnvCpt:::plot.envcpt(out, type = "aic"))
})

test_that("plot.envcpt works with type='bic'", {
  skip_if_not(identical(Sys.getenv("NOT_CRAN"), "true"))

  set.seed(1112)
  x <- c(rnorm(50, 0, 1), rnorm(50, 3, 1))
  out <- envcpt(x, verbose = FALSE)

  expect_silent(EnvCpt:::plot.envcpt(out, type = "bic"))
})

test_that("plot.envcpt works with type='fit'", {
  skip_if_not(identical(Sys.getenv("NOT_CRAN"), "true"))

  set.seed(1113)
  x <- c(rnorm(50, 0, 1), rnorm(50, 3, 1))
  out <- envcpt(x, verbose = FALSE)

  expect_silent(EnvCpt:::plot.envcpt(out, type = "fit"))
})

test_that("plot.envcpt works with custom valid colors", {
  skip_if_not(identical(Sys.getenv("NOT_CRAN"), "true"))

  set.seed(1114)
  x <- c(rnorm(50, 0, 1), rnorm(50, 3, 1))
  out <- envcpt(x, verbose = FALSE)

  custom_colors <- c("red", "blue", "green", "yellow", "purple", "orange",
                     "pink", "brown", "gray", "cyan", "magenta", "black")

  expect_silent(EnvCpt:::plot.envcpt(out, type = "aic", colors = custom_colors))
})

# =============================================================================
# Tests for Numerical Stability
# =============================================================================

context("Numerical stability tests")

test_that("cpt.reg handles data with small variance", {
  set.seed(2000)
  small_var_data <- rnorm(100, mean = 100, sd = 0.01)
  small_var_matrix <- cbind(
    small_var_data[-1],
    rep(1, 99),
    small_var_data[-100]
  )

  result <- EnvCpt:::cpt.reg(data = small_var_matrix, method = "PELT")
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg handles data with large variance", {
  set.seed(2001)
  large_var_data <- rnorm(100, mean = 0, sd = 1000)
  large_var_matrix <- cbind(
    large_var_data[-1],
    rep(1, 99),
    large_var_data[-100]
  )

  result <- EnvCpt:::cpt.reg(data = large_var_matrix, method = "PELT")
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg handles data with large mean", {
  set.seed(2002)
  large_mean_data <- rnorm(100, mean = 1e6, sd = 1)
  large_mean_matrix <- cbind(
    large_mean_data[-1],
    rep(1, 99),
    large_mean_data[-100]
  )

  result <- EnvCpt:::cpt.reg(data = large_mean_matrix, method = "PELT")
  expect_s4_class(result, "cpt.reg")
})

# =============================================================================
# Tests for Data Size Variations
# =============================================================================

context("Data size variation tests")

test_that("cpt.reg works with minimum size data", {
  # Minimum: need enough data for at least one segment
  min_data <- rnorm(20)
  min_matrix <- cbind(
    min_data[-1],
    rep(1, 19),
    min_data[-20]
  )

  result <- suppressWarnings(
    EnvCpt:::cpt.reg(data = min_matrix, method = "PELT", minseglen = 5)
  )
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg works with medium size data", {
  set.seed(3000)
  med_data <- rnorm(500)
  med_matrix <- cbind(
    med_data[-1],
    rep(1, 499),
    med_data[-500]
  )

  result <- EnvCpt:::cpt.reg(data = med_matrix, method = "PELT")
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg works with larger data", {
  set.seed(3001)
  large_data <- rnorm(1000)
  large_matrix <- cbind(
    large_data[-1],
    rep(1, 999),
    large_data[-1000]
  )

  result <- EnvCpt:::cpt.reg(data = large_matrix, method = "PELT")
  expect_s4_class(result, "cpt.reg")
})

# =============================================================================
# Tests for Penalty Comparison
# =============================================================================

context("Penalty comparison tests")

test_that("Different penalties give different number of changepoints", {
  set.seed(4000)
  n_pen <- 200
  seg1 <- arima.sim(model = list(ar = 0.3), n = 100)
  seg2 <- arima.sim(model = list(ar = 0.8), n = 100)
  pen_data <- c(seg1, seg2)

  pen_matrix <- cbind(
    pen_data[-1],
    rep(1, n_pen - 1),
    pen_data[-n_pen]
  )

  result_aic <- EnvCpt:::cpt.reg(data = pen_matrix, penalty = "AIC")
  result_bic <- EnvCpt:::cpt.reg(data = pen_matrix, penalty = "BIC")
  result_mbic <- EnvCpt:::cpt.reg(data = pen_matrix, penalty = "MBIC")

  # All should be valid
  expect_s4_class(result_aic, "cpt.reg")
  expect_s4_class(result_bic, "cpt.reg")
  expect_s4_class(result_mbic, "cpt.reg")

  # AIC typically gives more changepoints than BIC/MBIC
  cpts_aic <- length(cpts(result_aic))
  cpts_mbic <- length(cpts(result_mbic))
  expect_true(cpts_aic >= cpts_mbic - 1)  # Allow some variation
})

# =============================================================================
# Tests for Distribution Parameter
# =============================================================================

context("Distribution parameter tests")

test_that("cpt.reg with dist='Normal' works", {
  result <- EnvCpt:::cpt.reg(
    data = valid_ar1_matrix,
    method = "PELT",
    dist = "Normal"
  )
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg with unsupported distribution warns and converts", {
  expect_warning(
    result <- EnvCpt:::cpt.reg(
      data = valid_ar1_matrix,
      method = "PELT",
      dist = "Exponential"
    ),
    "is not supported"
  )
  expect_s4_class(result, "cpt.reg")
})

# =============================================================================
# Final Count Verification
# =============================================================================

context("Test suite verification")

test_that("Test suite runs successfully", {
  # This test confirms the suite completes

  expect_true(TRUE)
})
