test_that("nested solids and all physical coordinate permutations are correct", {
 z <- function(t)0; one <- function(t)1
 iz <- function(v)0; io <- function(v)1
 tetra <- vc_triple_integrate(function(v)1,c(0,1),z,function(x)1-x,iz,function(v)1-sum(v),n_outer=2,n_middle=2,n_inner=2)
 expect_equal(tetra$estimate,1/6)
 expect_true(tetra$tolerance_met)
 for(o in c("zyx","zxy","yzx","yxz","xzy","xyz")) {
   r <- vc_triple_integrate(function(v)v[1]+2*v[2]+3*v[3],c(0,1),z,function(t)2,iz,function(v)3,
     order=o,n_outer=2,n_middle=2,n_inner=2)
   idx <- match(rev(strsplit(o,"")[[1]]),c("x","y","z"))
   means <- numeric(3); means[idx] <- c(.5,1,1.5)
   expect_equal(r$estimate,6*sum(means*c(1,2,3)))
 }
 expect_equal(length(tetra$slices),5)
 expect_equal(tetra$slices[[1]],vc_double_integrate(function(v)1,c(0,1),iz,function(s)1-s,n_outer=2,n_inner=2))
 expect_equal(vc_triple_integrate(function(v)-1,c(1,0),z,one,iz,io,n_outer=2,n_middle=2,n_inner=2)$estimate,1)
 expect_equal(vc_triple_integrate(function(v)prod(v),c(0,1),z,one,iz,io,n_outer=2,n_middle=2,n_inner=2)$estimate,1/8)
})
test_that("triple uncertainty is retained and propagated", {
 calc <- function(n)vc_triple_integrate(function(v)exp(sum(v)),c(0,1),function(t)0,function(t)1,
   function(v)0,function(v)1,method="trapezoid",n_outer=n,n_middle=n,n_inner=n,abs_tol=1e-12,rel_tol=1e-12)
 a <- calc(2); b <- calc(4)
 expect_lt(abs(b$estimate-(exp(1)-1)^3),abs(a$estimate-(exp(1)-1)^3))
 expect_false(a$tolerance_met)
 expect_gt(a$slice_error,0)
 ids <- match(a$outer$table$x,a$outer_nodes)
 expect_equal(a$slice_error,sum(abs(a$outer$table$weight)*vapply(a$slices[ids],function(z)z$error_estimate,numeric(1))))
 expect_equal(a$error_estimate,a$outer_error+a$slice_error)
 expect_true(all(vapply(a$slices,function(z)length(z$inner)>0,logical(1))))
})
test_that("empty solids and invalid inputs are explicit", {
 bad <- function(...)stop("unexpected")
 z <- function(t)0; one <- function(t)1; iz <- function(v)0; io <- function(v)1
 expect_equal(vc_triple_integrate(bad,c(0,0),bad,bad,bad,bad)$estimate,0)
 expect_equal(vc_triple_integrate(bad,c(0,1),z,z,bad,bad)$estimate,0)
 expect_equal(vc_triple_integrate(bad,c(0,1),z,one,iz,iz)$estimate,0)
 expect_error(vc_triple_integrate(sum,c(0,1),one,z,iz,io),"Middle lower")
 expect_error(vc_triple_integrate(sum,c(0,1),z,one,io,iz),"lower")
 expect_error(vc_triple_integrate(sum,c(0,1),z,one,function(v)NA_real_,io))
 expect_error(vc_triple_integrate(function(v)Inf,c(0,1),z,one,iz,io))
 expect_error(vc_triple_integrate(function(v)c(1,2),c(0,1),z,one,iz,io))
 expect_error(vc_triple_integrate(sum,c(0,Inf),z,one,iz,io))
 for(n in c(0,1,3,2.5,Inf,NA)) expect_error(vc_triple_integrate(sum,c(0,0),z,one,iz,io,n_middle=n))
 expect_error(vc_triple_integrate(sum,c(0,0),z,one,iz,io,order="xxz"))
 expect_error(vc_triple_integrate(sum,c(0,0),z,one,iz,io,abs_tol=0,rel_tol=0))
 expect_error(vc_triple_integrate(NULL,c(0,0),z,one,iz,io))
})
test_that("triple graphics retain physical fiber endpoints and projection labels", {
 skip_if_not_installed("plotly")
 for(o in c("zyx","xyz")) {
 r <- vc_triple_integrate(sum,c(0,1),function(t)0,function(t)2,function(v)0,function(v)3,
   order=o,n_outer=2,n_middle=2,n_inner=2)
 p <- plotly::plotly_build(vc_plot_triple_integral(r))
 a <- p$x$data[[1]]
 coords <- cbind(as.numeric(a$x),as.numeric(a$y),as.numeric(a$z))
 idx <- match(rev(strsplit(o,"")[[1]]),c("x","y","z"))
 expect_equal(coords[1,],c(0,0,0))
 expected <- numeric(3); expected[idx] <- c(0,0,3)
 expect_equal(coords[2,],expected)
 expect_true(all(is.na(coords[3,])))
 expect_equal(a$type,"scatter3d")
 q <- plotly::plotly_build(vc_plot_triple_integral(r,"projection"))
 expect_equal(q$x$layout$xaxis$title,rev(strsplit(o,"")[[1]])[1])
 expect_equal(q$x$layout$yaxis$title,rev(strsplit(o,"")[[1]])[2])
 }
 expect_error(vc_plot_triple_integral(list()),"Invalid")
 r <- vc_triple_integrate(sum,c(0,0),identity,identity,sum,sum)
 expect_error(vc_plot_triple_integral(r),"No stored")
 r <- vc_triple_integrate(sum,c(0,1),identity,identity,sum,sum,n_outer=2)
 expect_error(vc_plot_triple_integral(r),"No inner")
})

