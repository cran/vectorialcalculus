test_that("interactive plots preserve mathematical anchor and increment geometry", {
  skip_if_not_installed("plotly")
  f <- function(x) x^2
  a <- vc_secant(f, 1, 2)
  p <- plotly::plotly_build(vc_plot_secant(a, f, n = 51))
  expect_s3_class(p, "plotly")
  expect_length(p$x$data, 3)
  expect_equal(as.numeric(p$x$data[[3]]$x), c(1, 2))
  expect_equal(as.numeric(p$x$data[[3]]$y), c(1, 4))
  line <- p$x$data[[2]]
  expect_equal(as.numeric(line$y), a$y0 + a$slope*(as.numeric(line$x)-a$x0))
  expect_identical(p$x$data[[1]]$connectgaps, FALSE)
  p <- plotly::plotly_build(vc_tangent(f, 1) |> vc_plot_tangent(f = f))
  expect_match(p$x$data[[2]]$name, "Approximate")
  expect_match(p$x$layout$title$text, "not proof")
  for (dx in c(-0.2, 0, 0.2)) {
    a <- vc_differential(f, 1, dx)
    p <- plotly::plotly_build(vc_plot_differential(a, f))
    expect_length(p$x$data, 6)
    expect_equal(as.numeric(p$x$data[[4]]$x), c(a$x0, a$x1))
    expect_equal(as.numeric(p$x$data[[5]]$y), c(a$y0, a$y0+a$dy))
    expect_equal(as.numeric(p$x$data[[6]]$y), c(a$y0, a$y1))
  }
})

test_that("plots label failed diagnostics and one-sided choices explicitly", {
  skip_if_not_installed("plotly")
  a <- vc_tangent(abs, 0, warn = FALSE)
  p <- plotly::plotly_build(vc_plot_tangent(a, abs))
  expect_match(p$x$layout$title$text, "FAILED")
  p <- plotly::plotly_build(vc_plot_tangent(vc_tangent(abs, 0, "forward"), abs))
  expect_match(p$x$layout$title$text, "forward")
})

test_that("curve sampling leaves invalid values as gaps and validates anchors", {
  skip_if_not_installed("plotly")
  f <- function(x) if (abs(x) < 0.2) NA_real_ else x^2
  a <- vc_secant(f, -1, 1)
  p <- plotly::plotly_build(vc_plot_secant(a, f, xlim = c(-1,1), n = 21))
  expect_true(anyNA(p$x$data[[1]]$y))
  expect_false(p$x$data[[1]]$connectgaps)
  expect_error(vc_plot_secant(a, sin), "does not match")
  g <- function(x) if (x == -1 || x == 1) 1 else NA_real_
  expect_error(vc_plot_secant(vc_secant(g, -1, 1), g, xlim = c(-2,2), n = 2), "No finite")
})

test_that("plot customization and invalid arguments are handled", {
  skip_if_not_installed("plotly")
  a <- vc_secant(sin, 0, 1)
  expect_error(vc_plot_tangent(a, sin), "must come from")
  for (bad in list(c(1,0), c(0,Inf), c(0,0.5), 1i, NA_real_))
    expect_error(vc_plot_secant(a, sin, xlim = bad), "xlim")
  for (bad in list(1, 2.5, Inf, NA_real_, "4"))
    expect_error(vc_plot_secant(a, sin, n = bad), "integer")
  expect_error(vc_plot_secant(a, sin, colors = "red"), "colors")
  expect_error(vc_plot_secant(a, sin, title = NA_character_), "labels")
  expect_error(vc_plot_secant(a, sin, xlab = 1), "labels")
  colors <- c(curve="red", line="blue", point="black", differential="green", change="purple")
  p <- plotly::plotly_build(vc_plot_secant(a, sin, colors = colors, title = "Example", xlab = "t"))
  expect_identical(p$x$data[[1]]$line$color, "red")
  expect_match(p$x$layout$title$text, "Example")
  expect_identical(p$x$layout$xaxis$title, "t")
  broken <- a; broken$slope <- Inf
  expect_error(vc_plot_secant(broken, sin), "nonfinite")
})

