test_that("phase plots preserve candidate positions and derivative evidence", {
  skip_if_not_installed("plotly")
  f <- function(x)x^3-3*x
  a <- vc_critical_points(f,c(-2,2),n=51)
  b <- plotly::plotly_build(vc_plot_critical_points(a,f,n=51))
  expect_s3_class(b,"plotly")
  expect_equal(as.numeric(b$x$data[[2]]$x),a$points$x)
  expect_equal(as.numeric(b$x$data[[2]]$y),a$points$value)
  expect_equal(as.numeric(b$x$data[[3]]$y),a$scan$derivative)
  expect_match(b$x$layout$title$text,"Numerical evidence")
  f <- function(x)x^3
  a <- vc_concavity(f,c(-1,1),n=51)
  b <- plotly::plotly_build(a |> vc_plot_concavity(f=f,n=51))
  expect_equal(as.numeric(b$x$data[[2]]$x),a$inflections$x)
  expect_equal(as.numeric(b$x$data[[3]]$y),a$scan$second_derivative)
  expect_match(paste(b$x$data[[3]]$text,collapse=" "),"concave_up")
  up <- b$x$data[[which(vapply(b$x$data,function(t) identical(t$name,"Concave up (sampled)"),logical(1)))]]
  down <- b$x$data[[which(vapply(b$x$data,function(t) identical(t$name,"Concave down (sampled)"),logical(1)))]]
  expect_true(all(up$y[is.finite(up$y)]>0))
  expect_true(all(down$y[is.finite(down$y)]<0))
  expect_error(vc_plot_concavity(a,sin),"does not match")
  expect_error(vc_plot_concavity(a,f,n=1),"n")
  expect_error(vc_plot_concavity(a,f,colors="red"),"colors")
  expect_error(vc_plot_concavity(a,f,title=NA_character_),"labels")
})

test_that("Newton plots animate real steps and handle zero/single steps", {
  skip_if_not_installed("plotly")
  f <- function(x)x^2-2
  a <- vc_newton(f,1)
  b <- plotly::plotly_build(a |> vc_plot_newton(f=f))
  expect_length(b$x$frames,nrow(a$iterations)-1L)
  for(i in seq_along(b$x$frames)) {
    frame <- b$x$frames[[i]]$data[[1]]
    expect_equal(as.numeric(frame$x),a$iterations$x[c(i,i+1)])
    expect_equal(as.numeric(frame$y),c(a$iterations$fx[i],0))
  }
  b <- plotly::plotly_build(vc_plot_newton(a,f,animate=FALSE))
  expect_length(b$x$frames,0)
  g <- function(x)x-1
  b <- plotly::plotly_build(vc_plot_newton(vc_newton(g,0,df=function(x)1),g))
  expect_length(b$x$frames,0)
  g <- function(x)x^2
  b <- plotly::plotly_build(vc_plot_newton(vc_newton(g,0),g))
  expect_length(b$x$frames,0)
  expect_match(b$x$layout$title$text,"converged")
  expect_error(vc_plot_newton(a,f,animate=NA),"animate")
  expect_error(vc_plot_newton(a,f,frame_ms=0),"frame_ms")
  expect_error(vc_plot_newton(list(),f),"wrong analysis")
})
