test_that("scalar line integrals retain signed density under reversal", {
  curve <- function(t)c(t,2*t)
  a <- vc_line_integral(function(p)-p[1],curve,c(0,2))
  b <- vc_line_integral(function(p)-p[1],curve,c(2,0),
    velocity=function(t)c(1,2),method="simpson",n=4)
  expect_equal(a$estimate,-2*sqrt(5))
  expect_equal(b$estimate,a$estimate)
  expect_equal(b$samples$t[c(1,nrow(b$samples))],c(2,0))
  expect_equal(b$samples$density,-b$samples$x)
  expect_identical(b$derivative_source,"supplied")
  expect_true(a$derivative_sensitivity_ok)
  expect_equal(length(a$derivatives),nrow(a$samples))
})

test_that("work changes sign and supports explicit parameter dependence", {
  curve <- function(t)c(cos(t),sin(t))
  v <- function(t)c(-sin(t),cos(t))
  rotation <- function(p)c(-p[2],p[1])
  a <- vc_line_integral(rotation,curve,c(0,2*pi),"work",velocity=v)
  b <- vc_line_integral(rotation,curve,c(2*pi,0),"work",velocity=v)
  expect_equal(a$estimate,2*pi)
  expect_equal(b$estimate,-a$estimate)
  expect_equal(b$integral$estimate,a$estimate)
  work <- vc_line_integral(function(p,t)c(t,2*t,0),function(t)c(t,2*t,3*t),
    c(0,2),"work",breaks=1,time_dependent=TRUE,method="simpson",n=4)
  expect_equal(work$estimate,10)
  expect_length(work$integral$results,2)
  expect_equal(work$dimension,3L)
  expect_equal(work$samples$t[c(1,nrow(work$samples))],c(0,2))
})

test_that("zero intervals, validation and stored graphics are explicit", {
  never <- function(...)stop("not evaluated")
  zero <- vc_line_integral(never,never,c(0,0))
  expect_equal(zero$estimate,0)
  expect_equal(nrow(zero$samples),0)
  expect_error(vc_line_integral(identity,function(t)t,c(0,1)),"two or three")
  expect_error(vc_line_integral(identity,function(t)c(t,t),c(0,1)),"Scalar field")
  expect_error(vc_line_integral(function(p)1,function(t)c(t,t),c(0,1),"work"),"Field dimension")
  expect_error(vc_line_integral(function(p)1,function(t)c(t,t),c(0,1),
    velocity=function(t)c(1,1,1)),"Velocity dimension")
  skip_if_not_installed("plotly")
  r <- vc_line_integral(function(p)p[1],function(t)c(t,t^2),c(-1,1),
    method="simpson",n=4)
  p <- plotly::plotly_build(vc_plot_line_integral(r))
  expect_identical(vapply(p$x$data,function(d)d$type,character(1)),
    c("scatter3d","surface","scatter3d"))
  expect_equal(as.numeric(p$x$data[[2]]$z[2,]),r$samples$density)
  expect_error(vc_plot_line_integral(zero),"No quadrature samples")
})
