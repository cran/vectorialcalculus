test_that("Green compares oriented circulation with the planar rotation", {
  f <- function(p)c(-p[2],p[1])
  lo <- function(x)0
  hi <- function(x)x
  ccw <- vc_green(f,lo,hi,c(0,1),n_boundary=4,n_outer=2,n_inner=2)
  expect_equal(ccw$circulation,1)
  expect_equal(ccw$interior_integral,1)
  expect_true(ccw$within_tolerance)
  expect_true(ccw$derivative_sensitivity$interior)
  expect_gt(ccw$derivative_sensitivity$interior_evaluations,0)
  expect_equal(ccw$boundary$left$estimate,0)
  expect_equal(nrow(ccw$boundary$left$samples),0)
  expect_equal(ccw$boundary$upper$samples$x[c(1,9)],c(1,0))
  cw <- vc_green(f,lo,hi,c(0,1),rotation=function(p)2,orientation="clockwise",
    lower_derivative=function(x)0,upper_derivative=function(x)1,
    n_boundary=4,n_outer=2,n_inner=2)
  expect_equal(cw$circulation,-1)
  expect_equal(cw$interior_integral,-1)
  expect_identical(names(cw$boundary),c("left","upper","right","lower"))
  expect_equal(cw$boundary$upper$samples$x[c(1,9)],c(0,1))
  expect_true(is.na(cw$derivative_sensitivity$interior))
})

test_that("curved boundaries and discrepant supplied rotation are explicit", {
  f <- function(p)c(0,p[1]^2)
  lo <- function(x)0
  hi <- function(x)x^2
  result <- vc_green(f,lo,hi,c(0,1),lower_derivative=function(x)0,
    upper_derivative=function(x)2*x,rotation=function(p)2*p[1],
    boundary_method="adaptive",n_outer=2,n_inner=2)
  expect_equal(result$circulation,.5)
  expect_equal(result$interior_integral,.5)
  expect_true(result$within_tolerance)
  mismatch <- vc_green(f,lo,hi,c(0,1),rotation=function(p)0,
    boundary_method="adaptive",n_outer=2,n_inner=2)
  expect_false(mismatch$within_tolerance)
  expect_equal(mismatch$difference,.5)
  expect_error(vc_green(f,hi,lo,c(0,1)),"exceeds")
  expect_error(vc_green(f,lo,hi,c(1,0)),"strictly increasing")
  expect_error(vc_green(function(p)1,lo,hi,c(0,1),n_boundary=2),"two components")
})

test_that("the Green plot uses stored oriented segments", {
  skip_if_not_installed("plotly")
  lo <- function(x)0
  hi <- function(x)1
  f <- function(p)c(-p[2],p[1])
  result <- vc_green(f,lo,hi,c(0,1),rotation=function(p)2,
    n_boundary=4,n_outer=2,n_inner=2)
  region <- vc_region(lo,hi,c(0,1))
  a <- plotly::plotly_build(vc_plot_green(result,region))$x$layout$annotations
  expect_length(a,5)
  expect_gt(a[[1]]$x,a[[1]]$ax)
  expect_gt(a[[2]]$y,a[[2]]$ay)
  expect_lt(a[[3]]$x,a[[3]]$ax)
  expect_lt(a[[4]]$y,a[[4]]$ay)
  expect_match(a[[5]]$text,"Circulation = 2")
  expect_error(vc_plot_green(result,vc_region(lo,hi,c(0,2))),"must match")
})
