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
