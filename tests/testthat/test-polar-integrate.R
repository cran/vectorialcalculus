test_that("polar Jacobian, Cartesian field and signed orientation are correct", {
  z <- function(t)0
  one <- function(t)1
  disk <- vc_polar_integrate(function(v)1,c(0,2*pi),z,one)
  expect_equal(disk$estimate,pi)
  expect_true(disk$tolerance_met)
  expect_equal(vc_polar_integrate(function(v)sum(v^2),c(0,2*pi),z,one)$estimate,pi/2)
  expect_equal(vc_polar_integrate(function(v)v[1],c(0,pi/2),z,one)$estimate,1/3,tolerance=1e-6)
  expect_equal(vc_polar_integrate(function(v)-1,c(0,pi/2),one,function(t)2)$estimate,-3*pi/4)
  expect_equal(vc_polar_integrate(function(v)1,c(2*pi,0),z,one)$estimate,-pi)
  expect_equal(vc_polar_integrate(function(v)1,c(0,pi),z,function(t)sin(t))$estimate,pi/4,tolerance=1e-6)
  nested <- vc_double_integrate(function(v)v[2],c(0,2*pi),z,one)
  expect_equal(disk$outer,nested$outer)
  expect_equal(disk$inner,nested$inner)
  expect_equal(disk$sections,nested$sections)
  expect_equal(disk$error_estimate,nested$error_estimate)
  expect_equal(disk$type,"polar_integral")
})
test_that("polar failures and empty regions retain explicit semantics", {
  bad <- function(...)stop("unexpected evaluation")
  z <- function(t)0
  one <- function(t)1
  expect_equal(vc_polar_integrate(bad,c(0,0),bad,bad)$estimate,0)
  expect_equal(vc_polar_integrate(bad,c(0,1),one,one)$estimate,0)
  expect_error(vc_polar_integrate(sum,c(0,7),z,one),"span")
  expect_error(vc_polar_integrate(sum,c(0,1),function(t)-1,one),"nonnegative")
  expect_error(vc_polar_integrate(sum,c(0,1),one,z),"lower")
  expect_error(vc_polar_integrate(sum,c(0,Inf),z,one),"finite")
  expect_error(vc_polar_integrate(sum,c(0,1),z,function(t)NA_real_))
  expect_error(vc_polar_integrate(function(v)1/sqrt(sum(v^2)),c(0,1),z,one))
  expect_error(vc_polar_integrate(function(v)c(1,2),c(0,1),z,one))
  expect_error(vc_polar_integrate(function(v)1e308,c(0,1),one,function(t)3))
  expect_error(vc_polar_integrate(sum,c(0,0),z,one,n_inner=3))
  expect_error(vc_polar_integrate(sum,c(0,0),z,one,abs_tol=0,rel_tol=0))
  expect_error(vc_polar_integrate(NULL,c(0,0),z,one))
  expect_error(vc_polar_integrate(sum,c(0,0),z,one,method="adaptive"))
})
test_that("trapezoid refinement improves a nonpolynomial radial field", {
  calc <- function(n)vc_polar_integrate(function(v)exp(sqrt(sum(v^2))),c(0,1),
    function(t)0,function(t)1,method="trapezoid",n_outer=2,n_inner=n,
    abs_tol=1e-12,rel_tol=1e-12)
  a <- calc(2); b <- calc(16)
  expect_lt(abs(b$estimate-1),abs(a$estimate-1))
  expect_false(a$tolerance_met)
  expect_gt(a$inner_error,0)
  expect_equal(a$error_estimate,a$outer_error+a$inner_error)
})
test_that("polar geometry transforms stored radii and angles without evaluation", {
  skip_if_not_installed("plotly")
  for(interval in list(c(0,2*pi),c(pi/2,0))) {
    r <- vc_polar_integrate(function(v)1,interval,function(t)1,function(t)2)
    p <- plotly::plotly_build(vc_plot_polar_integral(r))
    d <- r$sections[order(r$sections$outer),]
    t <- c(d$outer,rev(d$outer),d$outer[1])
    radius <- c(d$upper,rev(d$lower),d$upper[1])
    expect_equal(as.numeric(p$x$data[[1]]$x),radius*cos(t))
    expect_equal(as.numeric(p$x$data[[1]]$y),radius*sin(t))
    expect_equal(length(p$x$data),2)
    expect_equal(p$x$data[[1]]$fill,"toself")
    expect_equal(p$x$layout$yaxis$scaleanchor,"x")
    # Polygon shoelace area excludes the annular hole (chord approximation).
    x <- as.numeric(p$x$data[[1]]$x); y <- as.numeric(p$x$data[[1]]$y)
    area <- abs(sum(x[-length(x)]*y[-1]-x[-1]*y[-length(y)])/2)
    expect_equal(area,1.5*abs(diff(interval)),tolerance=.01)
  }
  expect_error(vc_plot_polar_integral(list()),"Invalid")
  expect_error(vc_plot_polar_integral(vc_polar_integrate(sum,c(0,0),identity,identity)),"No stored")
})
