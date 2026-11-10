test_that("decreasing secants reuse correct endpoint geometry", {
  f <- function(x)x^2
  h <- c(.5,.25,.125)
  for(method in c("forward","backward","central")) {
    a <- f |> vc_secant_convergence(1,steps=h,method=method)
    expect_equal(a$table$slope,switch(method,forward=2+h,backward=2-h,central=rep(2,3)))
    expect_equal(a$table$slope_difference,a$table$slope-a$tangent$slope)
    expect_length(a$results,3)
    for(z in a$results) expect_equal(z$y0+z$slope*(z$x1-z$x0),z$y1)
  }
  a <- vc_secant_convergence(abs,0,method="central")
  expect_false(a$tangent$derivative$reliable)
  expect_true(all(a$table$slope==0))
  expect_equal(vc_secant_convergence(sin,0,steps=c(.1,.01))$table$slope,sin(c(.1,.01))/c(.1,.01))
  for(steps in list(numeric(),c(1,1),c(.1,.5),0,Inf,1i))
    expect_error(vc_secant_convergence(f,1,steps=steps))
  expect_error(vc_secant_convergence(f,1,steps=1e-320),"collapsed")
  expect_error(vc_secant_convergence(function(x)NA_real_,1))
})

test_that("central animation lines pass through both endpoints, not through x0", {
  skip_if_not_installed("plotly")
  f <- function(x)x^2
  a <- vc_secant_convergence(f,1,steps=c(.5,.25),method="central")
  p <- plotly::plotly_build(vc_plot_secant_convergence(a,f))
  expect_length(p$x$frames,2)
  for(i in seq_along(a$results)) {
    tr <- p$x$frames[[i]]$data
    line <- Filter(function(z)identical(z$name,"Stored secant"),tr)[[1]]
    pts <- Filter(function(z)identical(z$name,"Secant endpoints"),tr)[[1]]
    z <- a$results[[i]]
    expect_equal(as.numeric(line$y),z$y0+z$slope*(as.numeric(line$x)-z$x0))
    expect_equal(as.numeric(pts$y),as.numeric(pts$x)^2)
    expect_identical(pts$mode,"markers")
    expect_equal(z$y0+z$slope*(1-z$x0),1+a$steps[i]^2)
  }
  expect_length(plotly::plotly_build(vc_plot_secant_convergence(a,f,animate=FALSE))$x$frames,0)
  one <- vc_secant_convergence(f,1,steps=.5)
  expect_length(plotly::plotly_build(vc_plot_secant_convergence(one,f))$x$frames,0)
  expect_error(vc_plot_secant_convergence(a,function(x)x^2+1),"does not match")
  expect_error(vc_plot_secant_convergence(a,f,xlim=c(.9,1.1)),"cover")
  expect_error(vc_plot_secant_convergence(a,f,animate=NA))
})
