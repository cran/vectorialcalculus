test_that("secants reproduce the line and are invariant under endpoint reversal", {
  f <- function(x) x^2
  a <- vc_secant(f, 1, 2)
  b <- vc_secant(f, 2, 1)
  expect_equal(a$slope, 3)
  expect_equal(a$y0 + a$slope * a$dx, a$y1)
  expect_equal(a$slope, b$slope)
  expect_equal(a$dx, -b$dx)
  expect_equal((sin |> vc_secant(x0 = 0, x1 = pi/2))$slope, 2/pi)
  expect_equal(vc_secant(function(x) 7, -1, 2)$slope, 0)
})

test_that("tangents reuse derivative diagnostics and preserve warnings", {
  f <- function(x) x^3
  a <- vc_tangent(f, 2)
  expect_equal(a$slope, 12, tolerance = 1e-6)
  expect_identical(a$derivative, vc_derivative(f, 2))
  expect_equal(a$y0, 8)
  expect_equal((sin |> vc_tangent(x0 = 0))$slope, 1, tolerance = 1e-8)
  expect_warning(vc_tangent(abs, 0), "sensitivity")
  expect_false(vc_tangent(abs, 0, warn = FALSE)$derivative$reliable)
  expect_equal(vc_tangent(abs, 0, "backward")$slope, -1)
  expect_error(vc_tangent(sin, 0, h = 0), "positive finite")
})

test_that("differentials distinguish linear and actual changes", {
  f <- function(x) x^2
  for (dx in c(-0.2, 0, 0.2)) {
    a <- vc_differential(f, 1, dx)
    expect_equal(a$dy, 2*dx, tolerance = 1e-8)
    expect_equal(a$delta_y, 2*dx + dx^2, tolerance = 1e-12)
    expect_equal(a$remainder, dx^2, tolerance = 1e-8)
    expect_identical(a$tangent, vc_tangent(f, 1))
  }
  a <- vc_differential(f, 1, 0.2)
  b <- vc_differential(f, 1, 0.1)
  expect_equal(b$remainder/a$remainder, 0.25, tolerance = 1e-8)
  expect_equal((sin |> vc_differential(x0 = 0, dx = 0.1))$dy, 0.1, tolerance = 1e-8)
  expect_warning(vc_differential(abs, 0, 0.1), "sensitivity")
  expect_false(vc_differential(abs, 0, 0.1, warn = FALSE)$tangent$derivative$reliable)
  expect_equal(vc_differential(identity, 1, 7*.Machine$double.eps)$actual_dx,
               (1 + 7*.Machine$double.eps) - 1)
})

test_that("local models reject invalid endpoints, outputs and overflow", {
  expect_error(vc_secant(sin, 1, 1), "distinct")
  for (bad in list(NA_real_, Inf, 1i, c(1,2), "1")) {
    expect_error(vc_secant(sin, bad, 2), "finite real scalar")
    expect_error(vc_secant(sin, 0, bad), "finite real scalar")
    expect_error(vc_differential(sin, 0, bad), "dx")
  }
  expect_error(vc_secant(1, 0, 1), "function")
  expect_error(vc_secant(function(x) c(x,x), 0, 1), "finite real scalar")
  expect_error(vc_secant(identity, -.Machine$double.xmax, .Machine$double.xmax), "displacement")
  expect_error(vc_secant(function(x) if(x == 0) -1e308 else 1e308, 0, 1), "Ordinate")
  expect_error(vc_secant(function(x) if(x == 0) 0 else 1e308, 0, 1e-308), "slope")
  expect_error(vc_differential(identity, 1, 1e-30), "collapsed")
  expect_error(vc_differential(identity, .Machine$double.xmax/2, .Machine$double.xmax), "Endpoint")
  expect_error(vc_differential(function(x) 2*x, 0, 1e308), "finite real scalar")
  f <- function(x) if (x < 0) NA_real_ else x^2
  expect_error(vc_differential(f, 0, -1, method = "forward"), "finite real scalar")
})

test_that("calculations do not load graphical namespaces", {
  before <- intersect(c("plotly", "ggplot2"), loadedNamespaces())
  expect_type(vc_secant(sin, 0, 1), "list")
  expect_type(vc_tangent(sin, 0), "list")
  expect_type(vc_differential(sin, 0, 0.1), "list")
  expect_identical(intersect(c("plotly", "ggplot2"), loadedNamespaces()), before)
})
