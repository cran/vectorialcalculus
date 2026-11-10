test_that("graph integrals distinguish area, signed density and oriented flux", {
  graph <- function(p)2*p[1]+3*p[2]
  gradient <- function(p)c(2,3)
  area <- vc_surface_graph_integral(graph,c(1,0),c(0,2),n_outer=2,n_inner=2)
  expect_equal(area$estimate,2*sqrt(14))
  expect_true(area$derivative_sensitivity_ok)
  expect_true(area$integral$tolerance_met)
  expect_equal(area$samples$z,2*area$samples$x+3*area$samples$y)
  scalar <- vc_surface_graph_integral(graph,c(0,1),c(0,2),function(p)-p[1],
    gradient=gradient,orientation="down",n_outer=2,n_inner=2)
  expect_equal(scalar$estimate,-sqrt(14))
  expect_identical(scalar$derivative_source,"supplied")
  expect_true(is.na(scalar$derivative_sensitivity_ok))
  for (direction in c("up","down")) {
    flux <- vc_surface_graph_integral(graph,c(0,1),c(0,2),function(p)c(1,2,3),
      "flux",gradient,direction,n_outer=2,n_inner=2)
    expect_equal(flux$estimate,if(direction=="up")-10 else 10)
    expect_identical(flux$orientation,direction)
  }
})

test_that("parametric geometry and orientation are independent of bound order", {
  cylinder <- function(p)c(2*cos(p[1]),2*sin(p[1]),p[2])
  j <- function(p)cbind(c(-2*sin(p[1]),2*cos(p[1]),0),c(0,0,1))
  f <- function(p)c(p[1],p[2],0)
  area <- vc_surface_integral(cylinder,c(0,2*pi),c(0,3),jacobian=j,n_outer=2,n_inner=2)
  expect_equal(area$estimate,12*pi)
  flux <- vc_surface_integral(cylinder,c(2*pi,0),c(3,0),f,"flux",j,
    n_outer=2,n_inner=2)
  expect_equal(flux$estimate,24*pi)
  swapped <- vc_surface_integral(function(p)cylinder(rev(p)),c(0,3),c(0,2*pi),f,"flux",
    function(p)j(rev(p))[,2:1],n_outer=2,n_inner=2)
  expect_equal(swapped$estimate,-flux$estimate)
  degenerate <- vc_surface_integral(function(p)c(sum(p),0,0),c(0,1),c(0,1),
    n_outer=2,n_inner=2)
  expect_equal(degenerate$estimate,0)
  expect_true(all(degenerate$samples$area_factor==0))
  never <- function(...)stop("not evaluated")
  zero <- vc_surface_integral(never,c(0,1),c(1,1),n_outer=2,n_inner=2)
  expect_equal(zero$estimate,0)
  expect_equal(nrow(zero$samples),0)
  expect_error(vc_surface_integral(cylinder,c(0,1),c(0,1),type="flux"),"requires")
  expect_error(vc_surface_integral(cylinder,c(0,1),c(0,1),jacobian=function(p)diag(2),
    n_outer=2,n_inner=2),"3 by 2")
})

test_that("stored patches use rows in the second parameter and preserve traversal", {
  surface <- function(p)c(p[1],p[2],prod(p))
  patch <- vc_surface_sample(surface,c(2,0),c(-1,1),n_u=3,n_v=5)
  expect_equal(dim(patch$grid$x),c(5,3))
  expect_equal(patch$grid$x[1,],c(2,1,0))
  expect_equal(patch$grid$y[,1],seq(-1,1,length.out=5))
  expect_equal(patch$grid$z,outer(patch$v,patch$u))
  tube <- vc_cylindrical_surface(function(t)c(t,t^2,t/2),c(0,2),c(2,-1),
    c(0,0,4),n_t=3,n_s=4)
  expect_equal(tube$direction,c(0,0,1))
  expect_equal(tube$grid$z[,2],c(2,1,0,-1)+.5)
  expect_equal(tube$grid$x[2,],tube$curve$samples$x)
  expect_error(vc_cylindrical_surface(identity,c(0,1),c(0,1),c(0,0,0)),"nonzero")
  skip_if_not_installed("plotly")
  p <- plotly::plotly_build(vc_plot_surface(tube))
  expect_equal(p$x$data[[1]]$z,tube$grid$z,ignore_attr=TRUE)
  expect_identical(p$x$data[[2]]$name,"Generating curve")
  point <- vc_surface_sample(surface,c(0,0),c(0,0))
  expect_identical(plotly::plotly_build(vc_plot_surface(point))$x$data[[1]]$mode,"markers")
})
