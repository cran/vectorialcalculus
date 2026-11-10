test_that("3D stationary graph preserves y-row/x-column coordinates and stored height", {
 skip_if_not_installed("plotly")
 f <- function(v)v[1]^2+2*v[2]^2
 r <- vc_stationary_points(f,c(.1,.2),assume_c2=TRUE,warn=FALSE)
 before <- serialize(r,NULL)
 p <- plotly::plotly_build(vc_plot_stationary_surface(r,f,xlim=c(-1,2),ylim=c(-2,3),n=5))
 expect_identical(serialize(r,NULL),before)
 s <- p$x$data[[1]]
 expect_equal(s$z,outer(seq(-2,3,length.out=5),seq(-1,2,length.out=5),function(y,x)x^2+2*y^2),ignore_attr=TRUE)
 expect_equal(as.numeric(s$x),seq(-1,2,length.out=5))
 expect_equal(as.numeric(s$y),seq(-2,3,length.out=5))
 t <- Filter(function(t)identical(t$name,"Conditional minimum"),p$x$data)
 expect_length(t,1)
 expect_equal(as.numeric(t[[1]]$z),r$table$value)
 expect_equal(as.numeric(t[[1]]$x),as.numeric(r$points[,1]))
 expect_match(t[[1]]$text,"gradient residual")
})
test_that("saddle and degenerate candidates retain conditional labels", {
 skip_if_not_installed("plotly")
 for(case in c("saddle","inconclusive")) {
 f <- if(case=="saddle")function(v)v[1]^2-v[2]^2 else function(v)v[1]^2+v[2]^4
 r <- vc_stationary_points(f,c(0,0),assume_c2=TRUE,warn=FALSE)
 p <- plotly::plotly_build(vc_plot_stationary_surface(r,f,n=5))
 expect_true(any(vapply(p$x$data,function(z)identical(z$name,paste("Conditional",case)),logical(1))))
 expect_equal(r$table$classification,case)
 }
})
test_that("surface domain gaps and invalid arguments stay explicit", {
 skip_if_not_installed("plotly")
 f <- function(v){if(v[1] < -.5)stop("domain");v[1]^2+2*v[2]^2}
 r <- vc_stationary_points(f,c(0,0),assume_c2=TRUE,warn=FALSE)
 p <- plotly::plotly_build(vc_plot_stationary_surface(r,f,xlim=c(-1,1),ylim=c(-1,1),n=5))
 expect_true(any(!is.finite(p$x$data[[1]]$z)))
 expect_false(p$x$data[[1]]$connectgaps)
 expect_error(vc_plot_stationary_surface(list(),f),"two-dimensional")
 expect_error(vc_plot_stationary_surface(r,function(v)1+sum(v^2)),"match")
 expect_error(vc_plot_stationary_surface(r,f,n=1))
 expect_error(vc_plot_stationary_surface(r,f,opacity=2))
 expect_error(vc_plot_stationary_surface(r,f,colors=c(a="red")))
 expect_error(vc_plot_stationary_surface(r,f,xlim=c(1,0)))
})

