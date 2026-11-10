test_that("Taylor polynomials reproduce linear and quadratic references", {
  linear <- function(x) 3*x - 2
  a <- vc_taylor(linear, 0, order = 1, x = c(-2,0,2))
  expect_equal(a$coefficients, c(constant=-2, linear=3), tolerance=1e-7)
  expect_equal(a$comparison$error, rep(0,3), tolerance=1e-7)
  quad <- function(x) x^2
  b <- vc_taylor(quad, 0, order = 2, x = c(-2,0,2), warn=FALSE)
  expect_equal(b$coefficients, c(constant=0,linear=0,quadratic=1), tolerance=1e-12)
  expect_equal(b$polynomial(c(-2,0,2)), quad(c(-2,0,2)), tolerance=1e-12)
  expect_false(b$reliable)
  expect_true(b$second_derivative$reliable)
  expect_false(b$derivative$side_agreement)
  expect_equal(b$second_derivative$estimate, 2)
  expect_null(a$second_derivative)
  expect_null(vc_taylor(sin, 0, order=1)$comparison)
  # Shifted center and nonzero linear coefficient exercise local powers.
  b <- vc_taylor(quad, 2, x = c(1.8,2,2.2))
  expect_equal(b$coefficients, c(constant=4, linear=4, quadratic=1), tolerance=1e-5)
  expect_equal(b$comparison$error, rep(0,3), tolerance=1e-6)
})

test_that("trigonometric and exponential coefficients match analytic derivatives", {
  for (center in c(0, .7)) {
    for (f in list(sin, cos, exp)) {
      a <- vc_taylor(f, center, order=2, warn=FALSE)
      first <- if (identical(f,sin)) cos(center) else if (identical(f,cos)) -sin(center) else exp(center)
      second <- if (identical(f,exp)) exp(center) else -f(center)
      expect_equal(unname(a$coefficients[2]), first, tolerance=1e-6)
      expect_equal(unname(a$coefficients[3]), second/2, tolerance=1e-5)
      expect_identical(a$derivative, vc_derivative(f, center, warn=FALSE))
      expect_identical(a$reliable, !(identical(f,cos) && center == 0))
    }
  }
  expect_equal((sin |> vc_taylor(x0=0, order=1))$polynomial(.2), .2, tolerance=1e-8)
})

test_that("comparison reports distance, signed error and decreasing local error", {
  a <- vc_taylor(exp, 0, x=c(.1,.2,.4))
  expect_equal(a$comparison$error, exp(a$comparison$x)-a$polynomial(a$comparison$x))
  expect_equal(a$comparison$abs_error, abs(a$comparison$error))
  expect_equal(a$comparison$distance, c(.1,.2,.4))
  expect_true(all(diff(a$comparison$abs_error)>0))
  first <- vc_taylor(exp, 0, order=1, x=.2)
  expect_lt(a$comparison$abs_error[2], first$comparison$abs_error)
  expect_true(all(a$comparison$status == "ok"))
})

test_that("second differences use direct stencils and expose sensitivity", {
  f <- function(x) x^4
  a <- vc_taylor(f, 0, h2=1, warn=FALSE)
  expect_equal(a$second_derivative$diagnostics$estimate, 2*c(1,.5,.25)^2)
  expect_false(a$second_derivative$stable)
  expect_false(a$reliable)
  expect_warning(vc_taylor(f, 0, h2=1), "Second-derivative")
  expect_true(vc_taylor(f, 0, h2=1, atol=2)$second_derivative$stable)
  # Tiny steps cause cancellation for exp: a numerically stable zero does
  # not justify a second-derivative coefficient; the roundoff flag rejects it.
  a <- vc_taylor(exp, 0, h2=1e-10, warn=FALSE)
  expect_false(a$second_derivative$resolved)
  expect_false(a$reliable)
  expect_warning(vc_taylor(exp, 0, h2=1e-10), "Second-derivative")
  a <- vc_taylor(function(x) 1, 0, warn=FALSE)
  expect_equal(a$second_derivative$estimate, 0)
  expect_true(a$second_derivative$stable)
  expect_false(a$second_derivative$resolved)
  expect_match(a$note, "do not prove smoothness")
  expect_false(vc_taylor(abs, 0, warn=FALSE)$reliable)
})

test_that("domain restrictions and unavailable comparisons are explicit", {
  f <- function(x) if (x < 0) NA_real_ else x^2
  a <- vc_taylor(f, 0, order=1, method="forward", x=c(-1,0,1))
  expect_true(is.na(a$comparison$value[1]))
  expect_true(is.na(a$comparison$error[1]))
  expect_equal(a$comparison$status[1], "function_not_finite")
  expect_equal(a$comparison$status[2:3], c("ok","ok"))
  expect_error(vc_taylor(f, 0, order=2, method="forward"), "finite real scalar")
  g <- function(x) if (x > 1) stop("outside") else exp(x)
  a <- vc_taylor(g, 0, x=c(0,2))
  expect_true(is.na(a$comparison$error[2]))
  expect_equal(a$comparison$status[2], "function_not_finite")
})

test_that("Taylor validates inputs, stencils and overflowing evaluations", {
  for (order in list(0,3,1.5,NA_real_,1i,c(1,2),"2"))
    expect_error(vc_taylor(sin, 0, order=order), "order")
  for (h2 in list(0,-1,Inf,NA_real_,1i,c(1,2)))
    expect_error(vc_taylor(sin,0,h2=h2), "h2")
  for (x in list(numeric(),c(0,Inf),NA_real_,1i,"1",matrix(1,1,1)))
    expect_error(vc_taylor(sin,0,x=x), "finite real vector")
  expect_error(vc_taylor(sin,0,warn=NA), "warn")
  expect_error(vc_taylor(sin,1,h2=1e-30), "collapsed")
  expect_error(vc_taylor(sin,0,h2=5e-324), "underflow")
  expect_error(vc_taylor(identity,.Machine$double.xmax/2,h2=.Machine$double.xmax), "overflowed")
  a <- vc_taylor(function(x) x^2, 0, warn=FALSE)
  expect_error(a$polynomial(Inf), "finite real vector")
  expect_error(a$polynomial(1e308), "overflowed")
  a <- vc_taylor(function(x) x^2, 0, x=1e308, warn=FALSE)
  expect_equal(a$comparison$status, "polynomial_not_finite")
  expect_true(is.na(a$comparison$error))
})

test_that("Taylor calculation has no graphical dependency or class requirement", {
  before <- intersect(c("plotly","ggplot2"), loadedNamespaces())
  a <- vc_taylor(exp, 0)
  expect_identical(class(a), "list")
  expect_s3_class(vc_taylor(exp,0,x=0)$comparison, "data.frame")
  expect_identical(intersect(c("plotly","ggplot2"), loadedNamespaces()), before)
})


test_that("finite values with overflowing sampled error remain explicitly unavailable", {
  f <- function(x) if (x == -1) 1e308 else 1e308*x
  a <- vc_taylor(f,0,order=1,x=-1)
  expect_true(is.finite(a$comparison$value))
  expect_true(is.finite(a$comparison$approximation))
  expect_true(is.na(a$comparison$error))
  expect_identical(a$comparison$status,"error_overflow")
})
