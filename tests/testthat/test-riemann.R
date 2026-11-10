test_that("regular/custom partitions and signed rules agree with analytic sums", {
  p <- vc_partition(breaks=c(0,.2,1,2))
  expect_equal(p$cells$width,c(.2,.8,1))
  expect_equal(p$mesh,1)
  for(method in c("left","right","midpoint")) {
    a <- vc_riemann(function(x)3,c(0,2),method=method)
    expect_equal(a$estimate,6)
    expect_equal(sum(a$table$contribution),a$estimate)
  }
  expect_equal((identity |> vc_riemann(interval=c(-1,1)))$estimate,0,tolerance=1e-14)
  expect_equal(vc_riemann(identity,partition=p)$estimate,sum(p$cells$midpoint*p$cells$width))
  for(method in c("left","right","midpoint")) {
    a <- vc_riemann_convergence(function(x)x^2,c(0,2),method=method,reference=8/3)
    expect_true(all(diff(a$table$abs_error)<0))
  }
  a <- vc_riemann_convergence(sin,c(0,pi),reference=2)
  expect_true(all(diff(a$table$abs_error)<0))
  expect_true(all(is.na(vc_riemann_convergence(sin,c(0,1))$table$error)))
  # Interior extrema make endpoint sums unsuitable as Darboux bounds.
  a <- vc_riemann(function(x)sin(pi*x),c(0,1),n=1,method="left")
  expect_equal(a$estimate,0)
  expect_match(a$note,"not a Darboux")
})
test_that("invalid partitions and function samples fail explicitly", {
  for(b in list(c(0,0,1),c(0,Inf),c(1,0),NA_real_,c(-1e308,1e308)))
    expect_error(vc_partition(breaks=b))
  expect_error(vc_partition(c(0,1),n=1.5))
  expect_error(vc_partition(c(0,2),breaks=c(0,1)))
  expect_error(vc_riemann(function(x)NA_real_,c(0,1)))
  expect_error(vc_riemann_convergence(sin,c(0,1),n_values=c(4,2)))
})
test_that("rectangle geometry and animation encode signed contributions", {
  skip_if_not_installed("plotly")
  a <- vc_riemann(identity,c(-1,1),n=2)
  p <- plotly::plotly_build(vc_plot_riemann(a,identity))
  expect_length(p$x$data,3)
  expect_true(any(p$x$data[[2]]$y<0,na.rm=TRUE))
  expect_true(any(p$x$data[[3]]$y>0,na.rm=TRUE))
  a <- vc_riemann_convergence(function(x)x^2,c(0,2),n_values=c(2,4,8))
  p <- plotly::plotly_build(vc_plot_riemann(a,function(x)x^2))
  expect_length(p$x$frames,3)
})
