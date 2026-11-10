test_that("product partitions reuse validated scalar partitions", {
 p <- vc_partition2d(xbreaks=c(0,.2,1),ybreaks=c(-1,0,2))
 expect_equal(p$x,vc_partition(breaks=c(0,.2,1)))
 expect_equal(p$y,vc_partition(breaks=c(-1,0,2)))
 expect_equal(p$cells$area,c(.2,.8,.4,1.6))
 expect_equal(p$cells$ix,c(1L,2L,1L,2L))
 expect_equal(sum(p$cells$area),3)
 expect_equal(p$mesh,sqrt(.8^2+2^2))
 expect_error(vc_partition2d(c(1,0),c(0,1)))
 expect_error(vc_partition2d(xbreaks=c(0,1,1),ybreaks=c(0,1)))
 expect_error(vc_partition2d(c(0,1),c(0,1),nx=1.5))
 expect_error(vc_partition2d(c(0,1e200),c(0,1e200)),"area")
})
test_that("signed tags and sampled statistics have analytic rectangular references", {
 f <- function(v)v[1]+2*v[2]
 for(m in c("midpoint","lower_left","upper_right","sample_min","sample_max","sample_mean")) {
 r <- vc_riemann2d(f,c(0,2),c(0,1),nx=2,ny=1,method=m)
 expected <- switch(m,midpoint=4,lower_left=1,upper_right=7,sample_min=1,sample_max=7,sample_mean=4)
 expect_equal(r$estimate,expected)
 expect_equal(sum(r$table$contribution),r$estimate)
 expect_equal(r$table$contribution,r$table$value*r$table$area)
 }
 p <- vc_partition2d(xbreaks=c(0,.2,1),ybreaks=c(0,.7,2))
 expect_equal(vc_riemann2d(function(v)-2,partition=p)$estimate,-4)
 expect_equal(vc_riemann2d(f,partition=p)$estimate,5)
 # Reconstruct breaks rather than trust a modified cell table.
 p$cells$area <- 999
 expect_equal(vc_riemann2d(function(v)1,partition=p)$estimate,2)
 # Corners miss an interior maximum; sampled extrema are not Darboux bounds.
 g <- function(v)sin(pi*v[1])*sin(pi*v[2])
 a <- vc_riemann2d(g,c(0,1),c(0,1),nx=1,ny=1,method="sample_max",sample_n=1)
 expect_lt(a$estimate,1e-20)
 expect_equal(vc_riemann2d(g,c(0,1),c(0,1),nx=1,ny=1,method="sample_max",sample_n=2)$estimate,1)
 expect_match(a$note,"Darboux")
})
test_that("restricted domains use an explicit zero extension not cell clipping", {
 f <- function(v)if(v[2]>v[1])stop("outside") else 1
 r <- vc_riemann2d(f,c(0,1),c(0,1),nx=2,ny=2,lower=function(x)0,upper=identity)
 expect_equal(r$estimate,.75)
 expect_equal(r$region,"sampled_indicator_on_box")
 expect_equal(r$table$area,rep(.25,4))
 expect_true(any(!vapply(r$samples,function(q)all(q$inside),logical(1))))
 z <- vc_riemann2d_convergence(f,c(0,1),c(0,1),n_values=c(2,4,8),lower=function(x)0,upper=identity,reference=.5)
 expect_true(all(diff(z$table$abs_error)<0))
 expect_true(all(diff(z$table$mesh)<0))
 expect_equal(z$table$change,c(NA,diff(z$table$estimate)))
 expect_error(vc_riemann2d(sum,c(0,1),c(0,1),lower=identity),"both")
 expect_error(vc_riemann2d(sum,c(0,1),c(0,1),lower=function(x)1,upper=function(x)0),"lower")
 expect_error(vc_riemann2d(sum,c(0,1),c(0,1),lower=function(x)-1,upper=function(x)1),"box")
})
test_that("refinement and invalid evaluation failures are explicit", {
 r <- vc_riemann2d_convergence(function(v)sum(v^2),c(0,1),c(0,1),reference=2/3)
 expect_true(all(diff(r$table$abs_error)<0))
 expect_equal(r$table$cells,r$table$n^2)
 expect_error(vc_riemann2d(function(v)Inf,c(0,1),c(0,1)))
 expect_error(vc_riemann2d(function(v)c(1,2),c(0,1),c(0,1)))
 expect_error(vc_riemann2d(sum,c(0,1),c(0,1),sample_n=0))
 expect_error(vc_riemann2d(sum,c(0,1),c(0,1),partition=list()))
 expect_error(vc_riemann2d_convergence(sum,c(0,1),c(0,1),n_values=c(2,2)))
 expect_error(vc_riemann2d_convergence(sum,c(0,1),c(0,1),reference=Inf))
})
test_that("mesh vertices preserve each stored signed contribution", {
 skip_if_not_installed("plotly")
 r <- vc_riemann2d(function(v)v[1]-2*v[2],c(0,2),c(0,1),nx=2,ny=2)
 p <- plotly::plotly_build(vc_plot_riemann2d(r))
 expect_equal(length(p$x$data),4)
 for(i in seq_len(4)) {
 t <- p$x$data[[i]];z <- r$table[i,]
 expect_equal(range(as.numeric(t$x)),c(z$x0,z$x1))
 expect_equal(range(as.numeric(t$y)),c(z$y0,z$y1))
 expect_equal(as.numeric(t$z),c(rep(0,4),rep(z$value,4)))
 expect_equal(length(t$i),12)
 expect_equal((max(t$x)-min(t$x))*(max(t$y)-min(t$y))*t$z[5],z$contribution)
 }
 expect_error(vc_plot_riemann2d(list()),"Invalid")
})
