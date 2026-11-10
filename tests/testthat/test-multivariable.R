test_that("partials and gradients reuse the scalar engine in arbitrary dimension", {
  f <- function(v)v[1]^2+2*v[2]^2
  p <- c(x=1,y=2)
  for(method in c("central","forward","backward")) {
    a <- f |> vc_partial(p,2,method=method)
    expect_equal(a$estimate,8,tolerance=1e-6)
    expect_identical(a$derivative$method,method)
    g <- vc_gradient(f,p,method=method)
    expect_equal(unname(g$estimate),c(2,8),tolerance=1e-6)
    expect_identical(names(g$estimate),names(p))
    expect_true(g$reliable)
    expect_length(g$results,2)
  }
  f4 <- function(v)sum(seq_along(v)*v^2)
  g <- vc_gradient(f4,c(1,2,3,4))
  expect_equal(unname(g$estimate),c(2,8,18,32),tolerance=1e-6)
  g <- vc_gradient(function(v)exp(v[1])+v[2]^2,c(0,1),method=c("forward","central"),h=c(1e-6,1e-4))
  expect_equal(unname(g$estimate),c(1,2),tolerance=1e-5)
  expect_equal(g$table$h,c(1e-6,1e-4))
  g <- vc_gradient(function(v)v["x"]^2+v["y"]^2,c(x=1,y=2))
  expect_equal(unname(g$estimate),c(2,4),tolerance=1e-6)
  expect_equal(unname(vc_gradient(function(v)v^3,2)$estimate),12,tolerance=1e-6)
})

test_that("directional restriction estimates and scaled normalization are correct", {
  f <- function(v)v[1]^2+2*v[2]^2
  g <- vc_gradient(f,c(1,2))
  for(method in c("central","forward","backward")) {
    d <- f |> vc_directional(c(1,2),c(3,4),method=method,gradient=g)
    expect_equal(d$estimate,7.6,tolerance=1e-6)
    expect_equal(d$direction,c(.6,.8))
    expect_equal(d$gradient_dot,7.6,tolerance=1e-6)
    expect_identical(d$gradient_comparison,"compatible")
  }
  expect_equal(vc_directional(f,c(1,2),c(-3,-4))$estimate,-7.6,tolerance=1e-6)
  for(direction in list(c(1e308,1e308),c(1e-300,1e-300))) {
    d <- vc_directional(f,c(1,2),direction)
    expect_equal(d$direction,rep(1/sqrt(2),2))
    expect_equal(d$estimate,10/sqrt(2),tolerance=1e-6)
  }
  expect_identical(vc_directional(f,c(1,2),c(1,0))$gradient_comparison,"not_requested")
  f3 <- function(v)v[1]^2+2*v[2]^2+3*v[3]^2
  expect_equal(vc_directional(f3,c(1,2,-1),c(1,0,1))$estimate,-4/sqrt(2),tolerance=1e-6)
})

test_that("coordinate checks do not force a gradient formula at nonsmooth points", {
  f <- function(v)if(all(v==0))0 else v[1]^3/sum(v^2)
  g <- vc_gradient(f,c(0,0))
  expect_true(g$reliable)
  expect_equal(unname(g$estimate),c(1,0))
  d <- vc_directional(f,c(0,0),c(1,1),gradient=g)
  expect_true(d$reliable)
  expect_equal(d$estimate,1/(2*sqrt(2)))
  expect_equal(d$gradient_dot,1/sqrt(2))
  expect_identical(d$gradient_comparison,"mismatch")
  cusp <- function(v)if(all(v==0))0 else prod(v)/sqrt(sum(v^2))
  g <- vc_gradient(cusp,c(0,0))
  expect_true(g$reliable)
  d <- vc_directional(cusp,c(0,0),c(1,1),gradient=g,warn=FALSE)
  expect_false(d$reliable)
  expect_identical(d$gradient_comparison,"inconclusive")
  expect_match(g$note,"do not establish")
})

test_that("invalid vectors, directions and unrepresentable stencils are rejected", {
  f <- function(v)sum(v^2)
  for(point in list(numeric(),c(1,Inf),NA_real_,matrix(1),1i,c(x=1,x=2)))
    expect_error(vc_gradient(f,point))
  expect_error(vc_partial(f,c(1,2),3))
  expect_error(vc_partial(f,c(1,2),1.5))
  expect_error(vc_gradient(f,c(1,2),h=c(1,2,3)))
  expect_error(vc_gradient(f,c(1,2),method=c("central","bad")))
  expect_error(vc_gradient(function(v)NA_real_,c(1,2)))
  expect_error(vc_gradient(function(v)c(1,2),c(1,2)))
  expect_error(vc_directional(f,c(1,2),c(0,0)))
  expect_error(vc_directional(f,c(1,2),1))
  expect_error(vc_directional(f,c(1,2),c(1e308,1e-308)),"underflow")
  expect_error(vc_directional(f,c(1,2),c(1,1),h=1e-320),"collapsed")
  expect_error(vc_directional(f,c(1,2),c(1,1e-300)),"collapsed")
  expect_error(vc_directional(f,c(1,2),c(1,1),gradient=vc_gradient(f,c(2,2))),"same point")
  expect_error(suppressWarnings(vc_partial(function(v)sqrt(v[1])+v[2],c(0,1),1)))
})

test_that("gradient contours, vectors and directional frames preserve coordinates", {
  skip_if_not_installed("plotly")
  f <- function(v)v[1]^2+2*v[2]^2
  point <- c(1,2)
  g <- vc_gradient(f,point)
  dirs <- lapply(list(c(1,0),c(0,1),c(-1,0)),function(u)vc_directional(f,point,u,gradient=g))
  p <- plotly::plotly_build(vc_plot_gradient(g,f,directions=dirs,n=7))
  contour <- p$x$data[[1]]
  expect_equal(as.numeric(contour$z[2,3]),as.numeric(contour$x[3]^2+2*contour$y[2]^2))
  gradient <- Filter(function(z)identical(z$name,"Numerical gradient"),p$x$data)[[1]]
  expect_equal(as.numeric(gradient$x),c(1,3),tolerance=1e-6)
  expect_equal(as.numeric(gradient$y),c(2,10),tolerance=1e-6)
  expect_length(p$x$frames,3)
  for(i in seq_along(dirs)) {
    a <- p$x$frames[[i]]$data[[1]]
    expect_equal(c(a$x[2]-a$x[1],a$y[2]-a$y[1]),dirs[[i]]$direction)
  }
  expect_length(plotly::plotly_build(vc_plot_gradient(g,f,directions=dirs,animate=FALSE,n=5))$x$frames,0)
  g3 <- vc_gradient(function(v)sum(v^2),c(1,2,3))
  p3 <- plotly::plotly_build(vc_plot_gradient(g3,function(v)sum(v^2)))
  expect_identical(p3$x$data[[1]]$type,"scatter3d")
  expect_equal(as.numeric(p3$x$data[[1]]$z),c(3,9),tolerance=1e-6)
  f3 <- function(v)sum(v^2)
  d3 <- lapply(list(c(1,0,0),c(0,0,1)),function(u)vc_directional(f3,c(1,2,3),u,gradient=g3))
  p3 <- plotly::plotly_build(vc_plot_gradient(g3,f3,directions=d3))
  expect_length(p3$x$frames,2)
  for(i in seq_along(d3)) {
    a <- p3$x$frames[[i]]$data[[1]]
    expect_equal(c(a$x[2]-a$x[1],a$y[2]-a$y[1],a$z[2]-a$z[1]),d3[[i]]$direction)
  }
  expect_error(vc_plot_gradient(g,function(v)f(v)+1),"does not match")
  expect_error(vc_plot_gradient(g,f,directions=list(vc_directional(f,c(2,2),c(1,0)))),"same point")
})

test_that("directional surface paths and tangent heights match stored results", {
  skip_if_not_installed("plotly")
  f <- function(v)v[1]^2+2*v[2]^2
  d <- vc_directional(f,c(1,2),c(3,4))
  p <- plotly::plotly_build(vc_plot_directional(d,f,n=9,t_range=c(-.5,.5)))
  surface <- p$x$data[[1]]
  expect_equal(as.numeric(surface$z[2,3]),as.numeric(surface$x[3]^2+2*surface$y[2]^2))
  curve <- Filter(function(z)identical(z$name,"Actual line restriction (sampled)"),p$x$data)[[1]]
  tangent <- Filter(function(z)identical(z$name,"Numerical directional tangent"),p$x$data)[[1]]
  expect_equal(as.numeric(curve$z),as.numeric(curve$x)^2+2*as.numeric(curve$y)^2)
  t <- (as.numeric(tangent$x)-1)/.6
  expect_equal(as.numeric(tangent$z),9+d$estimate*t,tolerance=1e-12)
  expect_error(vc_plot_directional(d,f,t_range=c(1,2)),"contain zero")
  expect_error(vc_plot_directional(d,f,alpha=2))
})
