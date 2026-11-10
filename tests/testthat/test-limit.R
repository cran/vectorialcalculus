test_that("finite limit probes distinguish removable, one-sided and unbounded behavior", {
  f <- function(x) { if (x==0) stop("do not evaluate x0"); sin(x)/x }
  a <- f |> vc_limit(x0=0)
  expect_identical(a$status,"compatible_finite")
  expect_equal(a$estimate,1,tolerance=1e-6)
  expect_false(any(a$samples$x==0))
  for(center in c(0,1)) {
    a <- vc_limit(identity,center)
    expect_identical(a$status,"compatible_finite")
    expect_equal(a$estimate,center,tolerance=1e-6)
  }
  a <- vc_limit(function(x)abs(x)/x,0)
  expect_identical(a$status,"side_disagreement")
  expect_true(is.na(a$estimate))
  expect_equal(a$left$estimate,-1)
  expect_equal(a$right$estimate,1)
  expect_equal(vc_limit(function(x)abs(x)/x,0,side="left")$estimate,-1)
  expect_equal(vc_limit(function(x)abs(x)/x,0,side="right")$estimate,1)
  a <- vc_limit(function(x)1/x,0)
  expect_identical(a$status,"unresolved")
  expect_true(is.na(a$estimate))
  a <- vc_limit(function(x)if(x<0)NA_real_ else x,0)
  expect_identical(a$status,"unavailable")
  expect_error(vc_limit(sin,1,h=1e-30),"collapsed")
  expect_error(vc_limit(sin,0,steps=3),"steps")
  expect_error(vc_limit(sin,0,h=0),"h")
  expect_error(vc_limit(sin,NA_real_),"finite real")
})

test_that("limit plots never evaluate or label f(x0) as a limit", {
  skip_if_not_installed("plotly")
  f <- function(x) { if(x==0) stop("not defined");sin(x)/x }
  a <- vc_limit(f,0)
  b <- plotly::plotly_build(vc_plot_limit(a,f))
  expect_s3_class(b,"plotly")
  expect_match(b$x$data[[4]]$name,"not f(x0)",fixed=TRUE)
  expect_equal(as.numeric(b$x$data[[4]]$y),a$estimate)
  expect_identical(b$x$layout$xaxis2$type,"log")
  expect_identical(b$x$layout$xaxis2$autorange,"reversed")
  expect_error(vc_plot_limit(a,cos),"does not match")
  expect_error(vc_plot_limit(a,f,xlim=c(1,2)),"approach point")
  g <- function(x)NA_real_
  expect_error(vc_plot_limit(vc_limit(g,0),g),"No finite")
})
