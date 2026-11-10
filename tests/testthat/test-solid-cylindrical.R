test_that("cylindrical patches map stored parameter grids to physical coordinates", {
 s <- vc_solid_cylindrical(function(t)1,function(t)2,function(v)v[1]-.2*v[2],
   function(v)2+v[1]+.3*v[2],c(0,pi/2),n=4,n_u=3,n_v=2)
 expect_equal(names(s$faces),c("theta_min","theta_max","radial_inner","radial_outer","z_lower","z_upper"))
 expect_equal(s$parameter,vc_solid_cartesian(function(t)1,function(t)2,function(v)v[1]-.2*v[2],
   function(v)2+v[1]+.3*v[2],c(0,pi/2),n=4,n_u=3,n_v=2))
 for(i in seq_len(6)) {
 f <- s$faces[[i]];raw <- s$parameter$faces[[i]]
 expect_equal(dim(f$x),dim(f$z))
 expect_equal(sqrt(f$x^2+f$y^2),raw$y,tolerance=1e-12)
 expect_equal(f$z,raw$z)
 }
 expect_equal(s$faces$theta_min$x[1,],seq(1,2,length.out=4))
 expect_true(all(s$faces$theta_min$y==0))
 expect_equal(s$faces$theta_max$y[1,],seq(1,2,length.out=4))
 expect_equal(s$faces$radial_inner$z[1,],s$faces$z_lower$z[1,])
 expect_equal(s$faces$theta_min$z[3,],s$faces$z_upper$z[,1])
 expect_false(s$full_turn)
})
test_that("axis degeneracies and angular seam patches remain explicit", {
 s <- vc_solid_cylindrical(function(t)0,function(t)1,function(v)0,function(v)2,
   c(0,2*pi),n=4,n_u=2,n_v=2)
 expect_true(s$full_turn)
 expect_true(all(s$faces$radial_inner$x==0 & s$faces$radial_inner$y==0))
 expect_equal(s$faces$theta_min$x,s$faces$theta_max$x,tolerance=1e-12)
 expect_equal(s$faces$theta_min$y,s$faces$theta_max$y,tolerance=1e-12)
 # A full-turn interval alone does not establish periodic boundary agreement.
 q <- vc_solid_cylindrical(function(t)0,function(t)1+t,function(v)0,function(v)1,
   c(0,2*pi),n=4,n_u=2,n_v=2)
 expect_gt(max(q$faces$theta_max$x),max(q$faces$theta_min$x))
 expect_match(q$note,"seams")
 p <- vc_partition(breaks=c(0,.1,.4,1))
 q <- vc_solid_cylindrical(function(t)1,function(t)2,function(v)0,function(v)1,partition=p,n_u=1,n_v=1)
 expect_equal(q$parameter$base$partition,p)
})
test_that("cylindrical geometry rejects invalid bounds without swaps", {
 a <- function(t)0;b <- function(t)1
 expect_error(vc_solid_cylindrical(a,b,a,b,c(0,7)),"span")
 expect_error(vc_solid_cylindrical(a,b,a,b,c(1,0)))
 expect_error(vc_solid_cylindrical(function(t)-1,b,a,b,c(0,1)),"nonnegative")
 expect_error(vc_solid_cylindrical(b,a,a,b,c(0,1)),"lower")
 expect_error(vc_solid_cylindrical(a,b,b,a,c(0,1)),"bottom")
 expect_error(vc_solid_cylindrical(a,b,a,function(v)Inf,c(0,1)))
 expect_error(vc_solid_cylindrical(a,b,a,b,partition=list()))
 expect_error(vc_solid_cylindrical(a,b,a,b,c(0,1),n_v=0))
})
test_that("cylindrical rendering uses stored selected patches only", {
 skip_if_not_installed("plotly")
 calls <- 0
 top <- function(v){calls <<- calls+1;2+v[2]}
 s <- vc_solid_cylindrical(function(t)1,function(t)2,function(v)0,top,c(0,2*pi),n=4,n_u=2,n_v=2)
 before <- calls
 f <- c("radial_inner","radial_outer","z_lower","z_upper")
 p <- plotly::plotly_build(vc_plot_solid_cylindrical(s,"both",faces=f,wire_step=99))
 expect_equal(calls,before)
 expect_equal(length(p$x$data),8)
 for(i in seq_along(f)) {
 t <- p$x$data[[2*i-1]]
 expect_equal(unname(t$x),s$faces[[f[i]]]$x,ignore_attr=TRUE)
 expect_equal(unname(t$y),s$faces[[f[i]]]$y,ignore_attr=TRUE)
 expect_equal(unname(t$z),s$faces[[f[i]]]$z,ignore_attr=TRUE)
 }
 expect_error(vc_plot_solid_cylindrical(list()),"Invalid")
 expect_error(vc_plot_solid_cylindrical(s,faces="invalid"),"selection")
})
