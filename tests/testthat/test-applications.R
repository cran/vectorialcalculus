test_that("geometric area differs from signed integration", {
  for(method in c("adaptive","trapezoid","simpson")) {
    a <- identity |> vc_area(c(-1,1),breaks=0,method=method)
    expect_equal(a$estimate,1)
    expect_equal(a$signed_estimate,0,tolerance=1e-12)
    expect_true(a$area$tolerance_met)
    expect_equal(vc_area(identity,c(1,-1),breaks=0,method=method)$estimate,1)
  }
  expect_equal(vc_area(function(x)x^2,c(0,2))$estimate,8/3)
  expect_equal(vc_area(sin,c(0,pi))$estimate,2)
  expect_equal(vc_area(sin,c(0,2*pi),breaks=pi)$estimate,4)
  expect_equal(vc_area(identity,c(0,1),g=function(x)x^2)$estimate,1/6)
  expect_equal(vc_area(function(x)x^2,c(0,1),g=identity)$estimate,1/6)
  expect_equal(vc_area(identity,c(-1,1),g=identity)$estimate,0)
  expect_equal(vc_area(function(x)stop("no evaluation"),c(1,1))$estimate,0)
  a <- vc_area(function(x)sin(100*x),c(0,1),subdivisions=1)
  expect_false(a$area$tolerance_met)
})

test_that("washers and shells give known volumes without overlap", {
  for(q in c("adaptive","trapezoid","simpson")) {
    a <- identity |> vc_revolution(c(0,1),quadrature=q,n=256)
    expect_equal(a$estimate,pi/3,tolerance=1e-5)
    expect_identical(a$method,"washers")
    expect_equal(vc_revolution(function(x)2,c(0,1),g=function(x)1,quadrature=q)$estimate,3*pi)
    # Straddling the axis produces a disk, not a difference of squared heights.
    expect_equal(vc_revolution(function(x)2,c(0,1),g=function(x)-1,quadrature=q)$estimate,4*pi)
    expect_equal(vc_revolution(function(x)3,c(0,1),g=function(x)1,axis_at=2,quadrature=q)$estimate,pi)
  }
  sphere <- function(x)sqrt(1-x^2)
  expect_equal(vc_revolution(sphere,c(-1,1))$estimate,4*pi/3)
  expect_equal(vc_revolution(sphere,c(-1,1),variable="y",axis="vertical")$estimate,4*pi/3)
  expect_equal(vc_revolution(function(x)3,c(0,1),axis_at=2)$estimate,4*pi)
  expect_equal(vc_revolution(function(x)3,c(0,1),g=function(x)2,axis_at=2)$estimate,pi)
  for(variable in c("x","y")) {
    axis <- if(variable=="x")"vertical" else "horizontal"
    a <- vc_revolution(function(t)1-t,c(0,1),variable=variable,axis=axis)
    expect_identical(a$method,"shells")
    expect_equal(a$estimate,pi/3)
    expect_equal(vc_revolution(function(t)1,c(-1,0),variable=variable,axis=axis)$estimate,pi)
  }
  expect_equal(vc_revolution(sphere,c(0,1),g=function(x)-sphere(x),axis="vertical")$estimate,4*pi/3)
  expect_equal(vc_revolution(identity,c(1,0))$estimate,pi/3)
  expect_equal(vc_revolution(identity,c(0,1),g=identity)$estimate,0)
  expect_equal(vc_revolution(function(t)stop("no evaluation"),c(1,1))$estimate,0)
  expect_error(vc_revolution(identity,c(-1,1),axis="vertical"),"double-count")
})

test_that("graph lengths and lateral/capped surfaces recover analytic cases", {
  for(method in c("adaptive","trapezoid","simpson")) {
    a <- identity |> vc_arc_length(c(0,1),method=method)
    expect_equal(a$estimate,sqrt(2),tolerance=1e-9)
    expect_true(a$derivative_sensitivity_ok)
    expect_equal(vc_arc_length(function(t)3,c(-1,1),method=method)$estimate,2)
    a <- identity |> vc_surface_area(c(0,1),df=function(t)1,method=method)
    expect_equal(a$estimate,pi*sqrt(2))
    expect_equal(a$cap_area,0)
    expect_true(is.na(a$derivative_sensitivity_ok))
  }
  exact <- (2*sqrt(5)+asinh(2))/4
  expect_equal(vc_arc_length(function(t)t^2,c(0,1),df=function(t)2*t)$estimate,exact)
  expect_equal(vc_arc_length(function(t)t^2,c(0,1))$estimate,exact,tolerance=1e-7)
  expect_equal(vc_arc_length(abs,c(-1,1),df=function(t)sign(t),breaks=0)$estimate,2*sqrt(2))
  expect_equal(vc_surface_area(function(t)2,c(-1,1),df=function(t)0)$estimate,8*pi)
  sphere <- function(t)sqrt(1-t^2)
  expect_equal(vc_surface_area(sphere,c(-1,1),df=function(t)-t/sqrt(1-t^2))$estimate,4*pi)
  expect_equal(vc_surface_area(function(t)3,c(0,2),axis_at=2,df=function(t)0)$estimate,4*pi)
  a <- vc_surface_area(identity,c(0,1),df=function(t)1,end_caps=TRUE)
  expect_equal(a$lateral_area,pi*sqrt(2))
  expect_equal(a$cap_area,pi)
  expect_equal(a$estimate,pi*(sqrt(2)+1))
  expect_equal(vc_arc_length(identity,c(1,0))$estimate,sqrt(2),tolerance=1e-9)
  expect_equal(vc_surface_area(function(t)stop("no evaluation"),c(1,1),end_caps=TRUE)$estimate,0)
  # A central quotient at a corner does not become a reliable graph derivative.
  a <- vc_arc_length(abs,c(-1,1),method="simpson",n=2)
  expect_false(a$derivative_sensitivity_ok)
})

test_that("application inputs and densities are validated", {
  for(b in list(c(0,0),c(1,0),c(-1,0),c(0,2),Inf,1i))
    expect_error(vc_area(identity,c(-1,1),breaks=b))
  expect_error(vc_area(identity,c(0,Inf)))
  expect_error(vc_area(identity,c(0,1),g=1))
  expect_error(vc_area(function(t)NA_real_,c(0,1)))
  expect_error(vc_revolution(identity,c(0,1),axis_at=Inf))
  expect_error(vc_revolution(function(t)1e308,c(0,1)))
  expect_error(vc_arc_length(identity,c(0,1),df=function(t)Inf))
  expect_error(vc_arc_length(identity,c(0,1),h=0))
  expect_error(vc_surface_area(identity,c(0,1),end_caps=NA))
  expect_error(vc_surface_area(identity,c(0,1),axis_at=Inf))
  expect_error(vc_revolution(identity,c(0,1),quadrature="simpson",n=3))
  expect_error(vc_arc_length(identity,c(0,1),abs_tol=0,rel_tol=0))
})

test_that("application refinement approaches known geometric measures", {
  f <- function(t)t^2
  a <- vc_area(f,c(0,1),method="trapezoid",n=4)
  b <- vc_area(f,c(0,1),method="trapezoid",n=8)
  expect_lt(abs(b$estimate-1/3),abs(a$estimate-1/3))
  a <- vc_revolution(identity,c(0,1),quadrature="trapezoid",n=4)
  b <- vc_revolution(identity,c(0,1),quadrature="trapezoid",n=8)
  expect_equal(abs(a$estimate-pi/3)/abs(b$estimate-pi/3),4)
  exact <- (2*sqrt(5)+asinh(2))/4
  a <- vc_arc_length(f,c(0,1),df=function(t)2*t,method="trapezoid",n=4)
  b <- vc_arc_length(f,c(0,1),df=function(t)2*t,method="trapezoid",n=8)
  expect_lt(abs(b$estimate-exact),abs(a$estimate-exact))
  expect_lt(b$integral$error_estimate,a$integral$error_estimate)
})

test_that("geometric area widgets display both envelopes and original graphs", {
  skip_if_not_installed("plotly")
  a <- vc_area(identity,c(-1,1),breaks=0)
  p <- plotly::plotly_build(vc_plot_area(a,identity,n=10))
  expect_length(p$x$data,4)
  expect_equal(as.numeric(p$x$data[[1]]$y),pmin(as.numeric(p$x$data[[1]]$x),0))
  expect_equal(as.numeric(p$x$data[[2]]$y),pmax(as.numeric(p$x$data[[2]]$x),0))
  expect_identical(p$x$data[[2]]$fill,"tonexty")
  expect_true(0 %in% p$x$data[[1]]$x)
  expect_error(vc_plot_area(a,identity,alpha=2))
})

test_that("revolution meshes, axes and section frames match their formulas", {
  skip_if_not_installed("plotly")
  for(variable in c("x","y")) {
    axis <- if(variable=="x")"horizontal" else "vertical"
    a <- vc_revolution(identity,c(0,1),variable=variable,axis=axis,axis_at=0)
    p <- plotly::plotly_build(vc_plot_revolution(a,identity,n=7,nt=13,sections=c(.25,.5,.75)))
    outer <- Filter(function(z)identical(z$name,"Outer radial boundary"),p$x$data)[[1]]
    axial <- if(axis=="horizontal")outer$x else outer$y
    radial <- if(axis=="horizontal")outer$y else outer$x
    expect_equal(radial^2+outer$z^2,axial^2,tolerance=1e-12,ignore_attr=TRUE)
    expect_length(p$x$frames,3)
    for(i in seq_along(p$x$frames)) {
      z <- p$x$frames[[i]]$data[[1]]
      u <- if(axis=="horizontal")z$x else z$y
      expect_true(all(abs(u[is.finite(u)]-c(.25,.5,.75)[i])<1e-12))
    }
  }
  a <- vc_revolution(function(t)2,c(0,1),g=function(t)1)
  p <- plotly::plotly_build(vc_plot_revolution(a,function(t)2,g=function(t)1,n=5,nt=9))
  inner <- Filter(function(z)identical(z$name,"Inner radial boundary"),p$x$data)[[1]]
  expect_equal(inner$y^2+inner$z^2,matrix(1,9,5),ignore_attr=TRUE)
  a <- vc_revolution(function(t)1,c(0,1),axis="vertical",axis_at=-1)
  p <- plotly::plotly_build(vc_plot_revolution(a,function(t)1,n=5,nt=9,sections=.5))
  upper <- Filter(function(z)identical(z$name,"Upper axial boundary"),p$x$data)[[1]]
  expect_equal(upper$y,matrix(1,9,5),ignore_attr=TRUE)
  expect_equal((upper$x+1)^2+upper$z^2,matrix(rep(seq(1,2,length.out=5)^2,each=9),9,5),ignore_attr=TRUE)
  expect_length(p$x$frames,0)
  p <- plotly::plotly_build(vc_plot_revolution(a,function(t)1,n=5,nt=9,sections=c(.25,.5),animate=FALSE))
  expect_length(p$x$frames,0)
  section <- Filter(function(z)identical(z$name,"Selected section"),p$x$data)[[1]]
  radii2 <- (section$x+1)^2+section$z^2
  expect_true(all(abs(radii2[is.finite(radii2)]-1.5^2)<1e-12))
  expect_error(vc_plot_revolution(a,function(t)1,sections=2))
  expect_error(vc_plot_revolution(a,function(t)1,nt=3))
  expect_error(vc_plot_revolution(a,function(t)1,show_caps=NA))
})
