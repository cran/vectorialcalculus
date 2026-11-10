test_that("accumulation respects analytic functions, orientation and input order", {
  x <- c(2,-1,0,1,2)
  for(method in c("adaptive","trapezoid","simpson")) {
    a <- (function(t)t^2) |> vc_accumulation(x=x,method=method,n=256)
    expect_identical(a$table$x,x)
    expect_equal(a$table$estimate,x^3/3,tolerance=1e-5)
    expect_equal(a$table$estimate[1],a$table$estimate[5])
    expect_length(a$results,length(x))
    expect_equal(vc_accumulation(identity,c(-1,0,1),lower=1,method=method)$table$estimate,
                 c(0,-.5,0),tolerance=1e-12)
    z <- vc_accumulation(function(t)stop("must not evaluate"),x=2,lower=2,method=method)
    expect_equal(z$table$estimate,0)
    expect_identical(z$table$status,"zero_interval")
  }
  expect_equal(vc_accumulation(sin,c(-pi,0,pi))$table$estimate,c(2,0,2))
  x <- c(-1,0,1)
  expect_equal(vc_accumulation(exp,x)$table$estimate,exp(x)-1)
  a <- vc_accumulation(function(t)sin(100*t),x=1,subdivisions=1)
  expect_identical(a$table$status,"adaptive_failure")
  expect_false(a$table$tolerance_met)
})

test_that("FTC compares actual numerical accumulation derivatives", {
  x <- c(-1,0,1,2)
  a <- (function(t)t^2) |> vc_fundamental(x=x)
  expect_equal(a$table$estimate,x^2,tolerance=1e-7)
  expect_true(all(a$table$status=="consistent"))
  expect_true(all(a$table$quadrature_resolved))
  expect_true(all(a$table$integration_ok))
  expect_length(a$results[[1]]$integrals,7)
  expect_s3_class(a$results[[1]]$diagnostics,"data.frame")
  for(method in c("central","forward","backward")) {
    a <- vc_fundamental(function(t)2,x=c(-1,0,1),derivative_method=method)
    expect_equal(a$table$estimate,rep(2,3),tolerance=1e-10)
    expect_true(all(a$table$status=="consistent"))
  }
  a <- vc_fundamental(sin,c(-1,1))
  expect_equal(a$table$estimate,sin(c(-1,1)),tolerance=1e-7)
  a <- vc_fundamental(exp,c(-1,0,1))
  expect_equal(a$table$estimate,exp(c(-1,0,1)),tolerance=1e-7)
  expect_true(all(a$table$status=="consistent"))
  # Simpson integrates the polynomial accumulation exactly (up to rounding).
  a <- vc_fundamental(function(t)t^2,x=1,method="simpson",n=2)
  expect_equal(a$table$estimate,1,tolerance=1e-7)
})

test_that("discontinuities, noise and domain failures never certify FTC", {
  jump <- function(t)if(t<0)-1 else 1
  a <- vc_fundamental(jump,x=0)
  expect_equal(a$table$estimate,0)
  expect_identical(a$table$status,"inconclusive")
  expect_false(a$table$derivative_reliable)
  spike <- function(t)if(t==0)100 else 1
  a <- vc_fundamental(spike,x=0)
  expect_equal(a$table$estimate,1,tolerance=1e-10)
  expect_equal(a$table$value,100)
  expect_identical(a$table$status,"mismatch")
  a <- vc_fundamental(exp,x=1,h=1e-10)
  expect_gt(a$table$quadrature_indicator,1e-8+1e-4*exp(1))
  expect_false(a$table$quadrature_resolved)
  expect_identical(a$table$status,"inconclusive")
  a <- vc_fundamental(function(t)t^2,x=1,method="trapezoid",n=2,h=.01)
  expect_false(a$table$integration_ok)
  expect_identical(a$table$status,"inconclusive")
  a <- suppressWarnings(vc_fundamental(sqrt,x=0))
  expect_null(a$results[[1]]$derivative)
  expect_true(is.na(a$table$estimate))
  expect_identical(a$table$status,"inconclusive")
  expect_match(a$table$message,"finite")
  a <- vc_fundamental(sqrt,x=0,derivative_method="forward",h=.001)
  expect_true(is.finite(a$table$estimate))
  expect_true(is.na(a$results[[1]]$derivative$side_agreement))
  a <- vc_fundamental(identity,x=1,h=1e-320)
  expect_identical(a$table$status,"inconclusive")
  expect_match(a$table$message,"collapsed")
})

test_that("invalid accumulation and FTC inputs fail before analysis", {
  for(x in list(numeric(),c(0,Inf),NA_real_,matrix(1),1i)) {
    expect_error(vc_accumulation(sin,x))
    expect_error(vc_fundamental(sin,x))
  }
  expect_error(vc_accumulation(sin,1,lower=Inf))
  expect_error(vc_accumulation(function(t)NA_real_,1))
  expect_error(vc_fundamental(function(t)NA_real_,1))
  expect_error(vc_fundamental(sin,1,h=0))
  expect_error(vc_fundamental(sin,1,atol=-1))
  expect_error(vc_fundamental(sin,1,method="simpson",n=3))
  expect_error(vc_fundamental(sin,1,subdivisions=0))
  expect_error(vc_fundamental(sin,1,abs_tol=0,rel_tol=0))
})

test_that("accumulation widgets retain signed geometry and synchronized markers", {
  skip_if_not_installed("plotly")
  f <- function(t)t^2
  a <- vc_accumulation(f,x=c(-1,0,1))
  q <- vc_fundamental(f,x=c(-1,0,1))
  p <- plotly::plotly_build(vc_plot_accumulation(a,f,fundamental=q,n=25))
  expect_length(p$x$frames,3)
  expect_true(any(vapply(p$x$data,function(z)identical(z$name,"Numerical A'(x)"),logical(1))))
  derivative <- Filter(function(z)identical(z$name,"Numerical A'(x)"),p$x$data)[[1]]
  expect_identical(derivative$mode,"markers")
  for(i in seq_along(p$x$frames)) {
    traces <- p$x$frames[[i]]$data
    pos <- Filter(function(z)identical(z$name,"positive heights"),traces)[[1]]
    neg <- Filter(function(z)identical(z$name,"negative heights"),traces)[[1]]
    mark <- Filter(function(z)identical(z$name,"Selected upper limit"),traces)[[1]]
    expect_equal(range(as.numeric(pos$x)),sort(c(0,a$table$x[i])))
    expect_equal(as.numeric(pos$y),as.numeric(pos$x)^2)
    expect_equal(as.numeric(neg$y),rep(0,length(neg$y)))
    expect_equal(as.numeric(mark$x),rep(a$table$x[i],2))
    expect_equal(as.numeric(mark$y),c(0,a$table$estimate[i]))
    expect_identical(mark$yaxis,"y2")
    expect_identical(mark$xaxis,"x")
    expect_identical(mark$mode,"lines+markers")
  }
  # Negative heights remain below zero, independent of integral orientation.
  a <- vc_accumulation(identity,x=-1)
  p <- plotly::plotly_build(vc_plot_accumulation(a,identity,n=11))
  expect_length(p$x$frames,0)
  neg <- Filter(function(z)identical(z$name,"negative heights"),p$x$data)[[1]]
  expect_equal(as.numeric(neg$y),as.numeric(neg$x))
  expect_gt(a$table$estimate,0)
  a <- vc_accumulation(f,x=0)
  expect_s3_class(vc_plot_accumulation(a,f),"plotly")
  expect_error(vc_plot_accumulation(a,f,alpha=2))
  expect_error(vc_plot_accumulation(a,f,animate=NA))
  expect_error(vc_plot_accumulation(a,f,xlim=c(1,2)))
  expect_error(vc_plot_accumulation(a,f,fundamental=q),"same lower")
  a <- vc_accumulation(f,x=c(-1,0,1))
  expect_error(vc_plot_accumulation(a,function(t)t^2+1,fundamental=q),"does not match")
})
