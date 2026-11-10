test_that("Gauss compares outward flux and divergence on the unit cube", {
  r <- vc_gauss(identity,c(0,1),c(0,1),c(0,1),h=.001,
    n_outer=2,n_middle=2,n_inner=2)
  expect_equal(r$outward_flux,3)
  expect_equal(r$volume_divergence,3)
  expect_equal(unname(r$face_flux),c(0,1,0,1,0,1))
  expect_true(r$within_tolerance)
  expect_true(r$derivative_sensitivity_ok)
  expect_equal(r$derivative_evaluations,125)
  expect_identical(r$derivative_source,"numerical")
  expect_true(all(vapply(r$faces,function(z)z$tolerance_met,logical(1))))
})

test_that("all six outward signs and coordinate widths work on a shifted box", {
  f <- function(p)p^2
  ranges <- list(c(1,2),c(-2,1),c(3,5))
  r <- vc_gauss(f,ranges[[1]],ranges[[2]],ranges[[3]],function(p)2*sum(p),
    n_outer=2,n_middle=4,n_inner=6)
  expect_equal(unname(r$face_flux),c(-6,24,-8,2,-27,75))
  expect_equal(r$volume_divergence,60)
  expect_equal(r$outward_flux,60)
  expect_equal(unname(r$normals),rbind(c(-1,0,0),c(1,0,0),c(0,-1,0),
    c(0,1,0),c(0,0,-1),c(0,0,1)))
  expect_true(is.na(r$derivative_sensitivity_ok))
  expect_identical(r$derivative_evaluations,0L)
  wrong <- vc_gauss(f,ranges[[1]],ranges[[2]],ranges[[3]],function(p)0,
    n_outer=2,n_middle=2,n_inner=2)
  expect_false(wrong$within_tolerance)
  expect_equal(wrong$difference,60)
})

test_that("zero total flux can arise from nonzero signed face contributions", {
  r <- vc_gauss(function(p)c(2,3,4),c(0,1),c(0,1),c(0,1),
    method="trapezoid",n_outer=1,n_middle=1,n_inner=1)
  expect_equal(unname(r$face_flux),c(-2,2,-3,3,-4,4))
  expect_equal(r$outward_flux,0)
  expect_equal(r$volume_divergence,0)
  expect_error(vc_gauss(identity,c(1,0),c(0,1),c(0,1)),"strictly increasing")
  expect_error(vc_gauss(identity,c(0,1),c(0,0),c(0,1)),"strictly increasing")
  expect_error(vc_gauss(function(p)c(1,2),c(0,1),c(0,1),c(0,1)),"three")
})

test_that("Gauss plots store outward normal vectors and report signed flux", {
  skip_if_not_installed("plotly")
  r <- vc_gauss(identity,c(0,1),c(0,1),c(0,1),function(p)3,
    n_outer=2,n_middle=2,n_inner=2)
  p <- plotly::plotly_build(vc_plot_gauss(r))
  arrows <- p$x$data[[7]]
  expect_identical(arrows$type,"cone")
  expect_equal(as.numeric(arrows$x),c(0,1,.5,.5,.5,.5))
  expect_equal(sign(as.numeric(arrows$u)),c(-1,1,0,0,0,0))
  expect_match(arrows$text[2],"x_max / flux = 1")
  expect_match(p$x$layout$annotations[[1]]$text,"Volume divergence = 3")
})
