test_that("Cartesian orders, signed fields and outer orientation agree analytically", {
  f <- function(v)v[1]+2*v[2]
  a <- vc_double_integrate(f,c(0,1),function(x)0,function(x)x,n_outer=2,n_inner=2)
  b <- vc_double_integrate(f,c(0,1),function(y)y,function(y)1,order="xy",n_outer=2,n_inner=2)
  expect_equal(a$estimate,2/3)
  expect_equal(b$estimate,2/3)
  expect_true(a$tolerance_met)
  expect_equal(vc_double_integrate(f,c(1,0),function(x)0,function(x)x)$estimate,-2/3)
  expect_equal(vc_double_integrate(function(v)-1,c(0,2),function(x)-1,function(x)2)$estimate,-6)
  expect_equal(length(a$inner),5)
  expect_equal(nrow(a$sections),5)
  expect_equal(a$outer$estimate,a$estimate)
  expect_equal(a$inner[[1]]$type,"integral")
  expect_equal(a$inner[[1]]$status,"zero_interval")
  expect_equal(a$error_estimate,a$outer_error+a$inner_error)
})

test_that("nested diagnostics preserve inner uncertainty and engine equivalence", {
  f <- function(v)exp(v[1]+v[2])
  r <- vc_double_integrate(f,c(0,1),function(x)0,function(x)1,
    method="trapezoid",n_outer=2,n_inner=2,abs_tol=1e-12,rel_tol=1e-12)
  expect_false(r$tolerance_met)
  expect_true(r$inner_error>0)
  expect_true(r$outer_error>0)
  for(i in seq_along(r$inner)) {
    t <- r$sections$outer[i]
    expect_equal(r$inner[[i]],vc_integrate(function(y)f(c(t,y)),c(0,1),
      "trapezoid",2,1e-12,1e-12))
  }
  d <- r$outer$table
  ids <- match(d$x,r$sections$outer)
  expect_equal(r$inner_error,sum(abs(d$weight)*r$sections$error_estimate[ids]))
  fine <- vc_double_integrate(f,c(0,1),function(x)0,function(x)1,
    method="trapezoid",n_outer=16,n_inner=16)
  expect_lt(abs(fine$estimate-(exp(1)-1)^2),abs(r$estimate-(exp(1)-1)^2))
  expect_equal(vc_double_integrate(function(v)sin(v[1])*cos(v[2]),c(0,pi),
    function(x)0,function(x)pi/2)$estimate,2,tolerance=2e-6)
})

test_that("empty regions and invalid inputs are explicit", {
  bad <- function(...)stop("must not evaluate")
  z <- vc_double_integrate(bad,c(1,1),bad,bad)
  expect_equal(z$estimate,0)
  expect_length(z$inner,0)
  expect_equal(z$status,"zero_interval")
  expect_equal(vc_double_integrate(bad,c(0,1),function(x)2,function(x)2)$estimate,0)
  expect_error(vc_double_integrate(sum,c(0,1),function(x)1,function(x)0),"lower")
  expect_error(vc_double_integrate(sum,c(0,Inf),identity,identity),"finite")
  expect_error(vc_double_integrate(sum,c(0,1),function(x)NA_real_,identity))
  expect_error(vc_double_integrate(function(v)Inf,c(0,1),function(x)0,function(x)1))
  expect_error(vc_double_integrate(function(v)c(1,2),c(0,1),function(x)0,function(x)1))
  for(n in c(0,1,3,2.5,NA,Inf))
    expect_error(vc_double_integrate(sum,c(0,0),identity,identity,n_inner=n))
  expect_error(vc_double_integrate(sum,c(0,0),identity,identity,abs_tol=0,rel_tol=0))
  expect_error(vc_double_integrate(sum,c(0,0),identity,identity,order="bad"))
  expect_error(vc_double_integrate(sum,c(0,0),identity,identity,method="adaptive"))
  expect_error(vc_double_integrate(sum,c(1e16,1e16+2),function(x)0,function(x)1),"grid|Partition|represent")
})

test_that("small refinement indicators do not certify the integral", {
  r <- vc_double_integrate(function(v)cos(8*pi*v[1])^2,c(0,1),
    function(x)0,function(x)1,n_outer=2,n_inner=2)
  expect_equal(r$estimate,1)
  expect_lt(r$error_estimate,1e-12)
  expect_true(r$tolerance_met)
  expect_gt(abs(r$estimate-.5),.49)
  # Independently checked curved region area and asymmetric rectangle.
  expect_equal(vc_double_integrate(function(v)1,c(0,1),function(x)0,
    function(x)x^2)$estimate,1/3)
  expect_equal(vc_double_integrate(function(v)v[1]+2*v[2],c(-1,2),
    function(y)0,function(y)3,order="xy")$estimate,22.5)
})

test_that("domain views use stored physical coordinates for either order", {
  skip_if_not_installed("plotly")
  for(o in c("yx","xy")) {
    r <- vc_double_integrate(sum,c(0,2),function(t)-1,function(t)t+1,order=o,n_outer=2,n_inner=2)
    p <- plotly::plotly_build(vc_plot_double_integral(r))
    d <- r$sections[order(r$sections$outer),]
    outer <- c(d$outer,rev(d$outer),d$outer[1])
    inner <- c(d$lower,rev(d$upper),d$lower[1])
    expect_equal(as.numeric(p$x$data[[1]]$x),if(o=="yx")outer else inner)
    expect_equal(as.numeric(p$x$data[[1]]$y),if(o=="yx")inner else outer)
    expect_equal(length(p$x$data),2)
    expect_equal(p$x$data[[1]]$fill,"toself")
  }
  expect_error(vc_plot_double_integral(list()),"Invalid")
  z <- vc_double_integrate(sum,c(0,0),identity,identity)
  expect_error(vc_plot_double_integral(z),"No stored")
})

