test_that("spherical shells obey the stored angle convention and preserve seams", {
 s <- vc_solid_spherical(function(t)0,function(t)pi,function(v)1,function(v)2,
   c(0,2*pi),n=4,n_u=4,n_v=2)
 expect_equal(names(s$faces),c("theta_min","theta_max","polar_lower","polar_upper","radial_inner","radial_outer"))
 expect_equal(s$parameter,vc_solid_cartesian(function(t)0,function(t)pi,function(v)1,function(v)2,
   c(0,2*pi),n=4,n_u=4,n_v=2))
 for(name in c("radial_inner","radial_outer")) {
 f <- s$faces[[name]]
 rho <- if(name=="radial_inner")1 else 2
 expect_equal(sqrt(f$x^2+f$y^2+f$z^2),matrix(rho,5,5),tolerance=1e-12)
 expect_true(all(f$x[c(1,5),]==0 & f$y[c(1,5),]==0))
 expect_equal(f$z[1,],rep(rho,5))
 expect_equal(f$z[5,],rep(-rho,5))
 }
 expect_equal(s$faces$radial_outer$x[3,1],2)
 expect_equal(s$faces$radial_outer$y[3,2],2)
 expect_equal(s$faces$theta_min$z[1,],s$faces$radial_inner$z[,1])
 expect_equal(s$faces$theta_max$x[3,],s$faces$radial_outer$x[,5])
 expect_equal(s$faces$polar_lower$z[1,],s$faces$radial_inner$z[1,])
 expect_true(s$full_turn)
})
test_that("angular sectors custom partitions and degeneracies remain explicit", {
 p <- vc_partition(breaks=c(0,.2,.5,1))
 s <- vc_solid_spherical(function(t).1*t,function(t)1+.1*t,function(v)0,function(v)1+v[2],
   partition=p,n_u=2,n_v=2)
 expect_equal(s$parameter$base$partition,p)
 expect_equal(s$parameter$base$sections$lower,.1*p$breaks)
 expect_true(all(s$faces$radial_inner$x==0 & s$faces$radial_inner$y==0 & s$faces$radial_inner$z==0))
 raw <- s$parameter$faces$z_upper;f <- s$faces$radial_outer
 expect_equal(f$x,raw$z*sin(raw$y)*cos(raw$x),tolerance=1e-12)
 expect_equal(f$y,raw$z*sin(raw$y)*sin(raw$x),tolerance=1e-12)
 expect_equal(f$z,raw$z*cos(raw$y),tolerance=1e-12)
 expect_match(s$note,"degeneracies")
 expect_false(s$full_turn)
})
test_that("spherical sampled ordering and angular bounds are validated", {
 a <- function(t)0;b <- function(t)1
 expect_error(vc_solid_spherical(a,b,a,b,c(0,7)),"span")
 expect_error(vc_solid_spherical(a,b,a,b,c(1,0)))
 expect_error(vc_solid_spherical(function(t)-.1,b,a,b,c(0,1)),"polar")
 expect_error(vc_solid_spherical(a,function(t)pi+.1,a,b,c(0,1)),"polar")
 expect_error(vc_solid_spherical(b,a,a,b,c(0,1)),"lower")
 expect_error(vc_solid_spherical(a,b,function(v)-1,b,c(0,1)),"nonnegative")
 expect_error(vc_solid_spherical(a,b,b,a,c(0,1)),"bottom")
 expect_error(vc_solid_spherical(a,b,a,function(v)NA_real_,c(0,1)))
 expect_error(vc_solid_spherical(a,b,a,b,c(0,1),n_v=0))
 expect_error(vc_solid_spherical(a,b,a,b,partition=list()))
})
test_that("shell rendering defaults to radial patches and never resamples", {
 skip_if_not_installed("plotly")
 calls <- 0
 top <- function(v){calls <<- calls+1;2}
 s <- vc_solid_spherical(function(t)0,function(t)pi,function(v)1,top,c(0,2*pi),n=4,n_u=4,n_v=2)
 before <- calls
 p <- plotly::plotly_build(vc_plot_solid_spherical(s,"both"))
 expect_equal(calls,before)
 expect_equal(length(p$x$data),4)
 for(i in seq_len(2)) {
 name <- c("radial_inner","radial_outer")[i];t <- p$x$data[[2*i-1]]
 expect_equal(unname(t$x),s$faces[[name]]$x,ignore_attr=TRUE)
 expect_equal(unname(t$y),s$faces[[name]]$y,ignore_attr=TRUE)
 expect_equal(unname(t$z),s$faces[[name]]$z,ignore_attr=TRUE)
 }
 expect_equal(length(plotly::plotly_build(vc_plot_solid_spherical(s,faces=names(s$faces)))$x$data),6)
 expect_error(vc_plot_solid_spherical(list()),"Invalid")
 expect_error(vc_plot_solid_spherical(s,faces="missing"),"selection")
})

