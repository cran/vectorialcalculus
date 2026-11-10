test_that("spherical Jacobian and coordinate convention match analytic moments", {
 calc <- function(f,interval=c(0,2*pi),polar=function(t)pi,lo=function(v)0,hi=function(v)1,n=8)
   vc_spherical_integrate(f,interval,function(t)0,polar,lo,hi,n_outer=8,n_middle=32,n_inner=n)
 a <- calc(function(v)1)
 expect_equal(a$estimate,4*pi/3,tolerance=1e-6)
 expect_equal(calc(function(v)1,polar=function(t)pi/2)$estimate,2*pi/3,tolerance=1e-6)
 expect_equal(calc(function(v)1,lo=function(v)1,hi=function(v)2)$estimate,28*pi/3,tolerance=1e-6)
 expect_equal(calc(function(v)sum(v^2),n=32)$estimate,4*pi/5,tolerance=3e-6)
 expect_equal(calc(function(v)v[3],polar=function(t)pi/2)$estimate,pi/4,tolerance=1e-6)
 expect_equal(calc(function(v)v[1],c(0,pi/2),function(t)pi/2)$estimate,pi/16,tolerance=3e-6)
 expect_equal(calc(function(v)v[2],c(0,pi/2),function(t)pi/2)$estimate,pi/16,tolerance=3e-6)
 expect_equal(calc(function(v)-1,c(2*pi,0))$estimate,4*pi/3,tolerance=1e-6)
 # Raw nested engine equivalence includes full diagnostic objects.
 r <- vc_spherical_integrate(function(v)1,c(0,1),function(t).2,function(t)1,
   function(v)1,function(v)2,n_outer=2,n_middle=2,n_inner=2)
 raw <- vc_triple_integrate(function(v)v[3]^2*sin(v[2]),c(0,1),function(t).2,function(t)1,
   function(v)1,function(v)2,n_outer=2,n_middle=2,n_inner=2)
 expect_equal(r$outer,raw$outer)
 expect_equal(r$slices,raw$slices)
 expect_equal(r$error_estimate,raw$error_estimate)
 expect_equal(r$type,"spherical_integral")
})
test_that("spherical refinement preserves heuristic uncertainty", {
 calc <- function(n)vc_spherical_integrate(function(v)1,c(0,1),function(t)0,function(t)pi,
   function(v)0,function(v)1,method="trapezoid",n_outer=2,n_middle=n,n_inner=2,abs_tol=1e-12,rel_tol=1e-12)
 a <- calc(2); b <- calc(16)
 # Radial trapezoid error is independent; isolate polar slice refinement directly.
 expect_lt(abs(b$slices[[1]]$outer$estimate-(2/3+1/48)),abs(a$slices[[1]]$outer$estimate-(2/3+1/48)))
 expect_false(a$tolerance_met)
 expect_gt(a$slice_error,0)
 expect_equal(a$error_estimate,a$outer_error+a$slice_error)
})
test_that("spherical validation and zero widths are explicit", {
 z <- function(t)0; one <- function(t)1; bad <- function(...)stop("unexpected")
 expect_equal(vc_spherical_integrate(bad,c(0,0),bad,bad,bad,bad)$estimate,0)
 expect_equal(vc_spherical_integrate(bad,c(0,1),z,z,bad,bad)$estimate,0)
 expect_equal(vc_spherical_integrate(bad,c(0,1),z,one,one,one)$estimate,0)
 expect_error(vc_spherical_integrate(sum,c(0,7),z,one,z,one),"span")
 expect_error(vc_spherical_integrate(sum,c(0,1),function(t)-.1,one,z,one),"Polar")
 expect_error(vc_spherical_integrate(sum,c(0,1),z,function(t)pi+.1,z,one),"Polar")
 expect_error(vc_spherical_integrate(sum,c(0,1),one,z,z,one),"lower")
 expect_error(vc_spherical_integrate(sum,c(0,1),z,one,function(v)-1,one),"nonnegative")
 expect_error(vc_spherical_integrate(sum,c(0,1),z,one,one,z),"lower")
 expect_error(vc_spherical_integrate(sum,c(0,1),z,one,function(v)NA_real_,one))
 expect_error(vc_spherical_integrate(function(v)1/sqrt(sum(v^2)),c(0,1),z,one,z,one))
 expect_error(vc_spherical_integrate(function(v)c(1,2),c(0,1),z,one,z,one))
 expect_error(vc_spherical_integrate(sum,c(0,0),z,one,z,one,n_middle=3))
 expect_error(vc_spherical_integrate(sum,c(0,0),z,one,z,one,abs_tol=0,rel_tol=0))
 expect_error(vc_spherical_integrate(NULL,c(0,0),z,one,z,one))
 expect_error(vc_spherical_integrate(sum,c(0,1),z,one,one,function(v)1e200))
})
test_that("spherical graphics retain physical endpoints and shell core", {
 skip_if_not_installed("plotly")
 r <- vc_spherical_integrate(sum,c(0,pi/2),function(t)0,function(t)pi,
   function(v)1,function(v)2,n_outer=2,n_middle=2,n_inner=2)
 p <- plotly::plotly_build(vc_plot_spherical_integral(r))
 a <- p$x$data[[1]]
 xyz <- cbind(as.numeric(a$x),as.numeric(a$y),as.numeric(a$z))
 expect_equal(xyz[1,],c(0,0,1))
 expect_equal(xyz[2,],c(0,0,2))
 expect_true(all(is.na(xyz[3,])))
 finite <- rowSums(is.finite(xyz))==3
 radius <- sqrt(rowSums(xyz[finite,,drop=FALSE]^2))
 expect_true(all(abs(radius-1)<1e-12 | abs(radius-2)<1e-12))
 # Independent physical formula for every stored fiber.
 expected <- do.call(rbind,lapply(seq_along(r$slices),function(i) {
   d <- r$slices[[i]]$sections
   ph <- as.vector(rbind(d$outer,d$outer,NA_real_));rho <- as.vector(rbind(d$lower,d$upper,NA_real_))
   si <- sin(ph);si[which(ph==0|ph==pi)] <- 0
   cbind(rho*si*cos(r$outer_nodes[i]),rho*si*sin(r$outer_nodes[i]),rho*cos(ph))
 }))
 expect_equal(xyz,expected[-nrow(expected),,drop=FALSE])
 expect_equal(a$type,"scatter3d")
 expect_error(vc_plot_spherical_integral(list()),"Invalid")
 z <- vc_spherical_integrate(sum,c(0,0),identity,identity,sum,sum)
 expect_error(vc_plot_spherical_integral(z),"No stored")
})

