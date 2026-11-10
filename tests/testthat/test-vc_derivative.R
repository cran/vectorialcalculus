test_that("polynomial and trigonometric reference derivatives agree", {
  for (method in c("forward", "backward", "central")) {
    expect_equal(vc_derivative(function(x) x^3, 2, method)$estimate, 12, tolerance = 1e-5)
    expect_equal(vc_derivative(sin, 0.7, method)$estimate, cos(0.7), tolerance = 1e-5)
    expect_equal(vc_derivative(function(x) 5, 0, method)$estimate, 0)
    expect_equal(vc_derivative(function(x) -3*x, -2, method)$estimate, -3, tolerance = 1e-6)
  }
  expect_equal((sin |> vc_derivative(x = 0))$estimate, 1, tolerance = 1e-8)
})

test_that("explicit stencils and sensitivity tolerances behave as documented", {
  a <- vc_derivative(function(x) x^2, 1, "forward", h = 0.1, warn = FALSE)
  b <- vc_derivative(function(x) x^2, 1, "backward", h = 0.1, warn = FALSE)
  c <- vc_derivative(function(x) x^2, 1, "central", h = 0.1, warn = FALSE)
  expect_equal(a$diagnostics$estimate, 2 + c(0.1, 0.05, 0.025))
  expect_equal(b$diagnostics$estimate, 2 - c(0.1, 0.05, 0.025))
  expect_equal(c$estimate, 2)
  expect_false(a$stable)
  expect_true(vc_derivative(function(x) x^2, 1, "forward", h = 0.1, atol = 0.1)$stable)
  expect_warning(vc_derivative(exp, 1, h = 1), "sensitivity")
  expect_false(vc_derivative(exp, 1, h = 1, warn = FALSE)$stable)
  expect_true(is.na(a$side_agreement))
})

test_that("cusps and discontinuities are not certified by central cancellation", {
  a <- vc_derivative(abs, 0, warn = FALSE)
  expect_equal(a$estimate, 0)
  expect_true(a$stable)
  expect_false(a$side_agreement)
  expect_false(a$reliable)
  expect_warning(vc_derivative(abs, 0), "sensitivity")
  expect_equal(vc_derivative(abs, 0, "forward")$estimate, 1)
  expect_equal(vc_derivative(abs, 0, "backward")$estimate, -1)
  jump <- function(x) if (x == 0) 1 else 0
  expect_false(vc_derivative(jump, 0, warn = FALSE)$reliable)
  expect_match(a$note, "do not establish differentiability")
})

test_that("domain boundaries respect the selected side", {
  f <- function(x) if (x < 0) NA_real_ else x^2
  expect_true(is.finite(vc_derivative(f, 0, "forward")$estimate))
  expect_error(vc_derivative(f, 0), "finite real scalar")
  g <- function(x) if (x > 0) NA_real_ else x^2
  expect_true(is.finite(vc_derivative(g, 0, "backward")$estimate))
})

test_that("invalid inputs and outputs are rejected", {
  expect_error(vc_derivative(1, 0), "function")
  for (x in list(NA_real_, Inf, NaN, 1i, numeric(), c(1,2), "1", TRUE))
    expect_error(vc_derivative(sin, x), "finite real scalar")
  for (h in list(0, -1, NA_real_, Inf, 1i, c(1,2), "1"))
    expect_error(vc_derivative(sin, 0, h = h), "positive finite")
  expect_error(vc_derivative(sin, 0, method = "other"))
  for (tol in list(-1, Inf, NA_real_, 1i, c(1,2))) {
    expect_error(vc_derivative(sin, 0, atol = tol), "atol")
    expect_error(vc_derivative(sin, 0, rtol = tol), "rtol")
  }
  for (w in list(NA, 1, c(TRUE,FALSE)))
    expect_error(vc_derivative(sin, 0, warn = w), "warn")
  for (value in list(NA_real_, Inf, NaN, 1i, c(1,2), "1", TRUE)) {
    f <- function(x) value
    expect_error(vc_derivative(f, 0), "finite real scalar")
  }
})

test_that("unrepresentable steps and overflowing quotients fail explicitly", {
  expect_error(vc_derivative(identity, 1, h = 1e-30), "collapsed")
  expect_error(vc_derivative(identity, .Machine$double.xmax, h = .Machine$double.xmax), "overflowed")
  expect_error(vc_derivative(identity, 0, h = 5e-324), "underflow")
  f <- function(x) if (x > 0) 1e308 else -1e308
  expect_error(vc_derivative(f, 0), "Nonfinite")
})

test_that("results are plain data without graphical dependencies", {
  graphics_before <- intersect(c("plotly", "ggplot2"), loadedNamespaces())
  a <- vc_derivative(sin, 1)
  expect_identical(class(a), "list")
  expect_s3_class(a$diagnostics, "data.frame")
  expect_equal(nrow(a$diagnostics), 3)
  expect_length(a$changes, 2)
  expect_true(a$reliable)
  expect_identical(intersect(c("plotly", "ggplot2"), loadedNamespaces()), graphics_before)
})

test_that("rounded displacements use their actual spacing", {
  x <- 1
  h <- 7 * .Machine$double.eps
  for (method in c("forward", "backward", "central")) {
    out <- vc_derivative(identity, x, method, h = h)
    expect_equal(out$estimate, 1, tolerance = 0)
    expect_true(out$stable)
  }
})

test_that("singular endpoint slopes fail sensitivity diagnostics", {
  out <- vc_derivative(sqrt, 0, "forward", h = 1e-4, warn = FALSE)
  expect_true(is.finite(out$estimate))
  expect_false(out$stable)
  expect_false(out$reliable)
  expect_warning(vc_derivative(sqrt, 0, "forward", h = 1e-4), "sensitivity")
})

test_that("passing diagnostics do not guarantee accuracy under cancellation", {
  # The exact real derivative is one, but these increments are lost in doubles.
  out <- vc_derivative(function(x) 1e20 + x, 0)
  expect_equal(out$estimate, 0)
  expect_true(out$reliable)
  expect_match(out$note, "do not establish differentiability or bound the error")
})
