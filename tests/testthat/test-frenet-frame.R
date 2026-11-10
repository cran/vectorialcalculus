test_that("the helix frame is orthonormal, oriented and consistent with numerical derivatives", {
  helix <- function(t)c(cos(t),sin(t),t/2)
  velocity <- function(t)c(-sin(t),cos(t),.5)
  acceleration <- function(t)c(-cos(t),-sin(t),0)
  exact <- vc_frenet_frame(helix,0,velocity,acceleration)
  expect_identical(exact$status,"ok")
  expect_equal(exact$point,c(1,0,0))
  expect_equal(exact$tangent,c(0,2,1)/sqrt(5))
  expect_equal(exact$normal,c(-1,0,0))
  expect_equal(exact$binormal,c(0,-1,2)/sqrt(5))
  basis <- cbind(exact$tangent,exact$normal,exact$binormal)
  expect_equal(crossprod(basis),diag(3))
  expect_equal(det(basis),1)

  t <- .3
  exact <- vc_frenet_frame(helix,t,velocity,acceleration)
  numeric <- vc_frenet_frame(helix,t,h2=.002)
  expect_identical(numeric$status,"ok")
  for (name in c("tangent","normal","binormal"))
    expect_equal(numeric[[name]],exact[[name]],tolerance=1e-6)
  expect_equal(numeric$diagnostics$velocity$type,"jacobian")
  expect_length(numeric$diagnostics$acceleration,3)
  expect_equal(numeric$diagnostics$acceleration[[1]]$h,.002)

  reverse <- vc_frenet_frame(function(s)helix(-s),-t,
    function(s)-velocity(-s),function(s)acceleration(-s))
  expect_equal(reverse$tangent,-exact$tangent)
  expect_equal(reverse$normal,exact$normal)
  expect_equal(reverse$binormal,-exact$binormal)
})

test_that("missing directions and ambiguous numerical results are explicit", {
  line <- function(t)c(2+t,4-3*t,1+2*t)
  straight <- vc_frenet_frame(line,.3)
  expect_identical(straight$status,"unresolved_curvature")
  expect_equal(straight$tangent,c(1,-3,2)/sqrt(14),tolerance=1e-7)
  expect_true(all(is.na(straight$normal)))
  expect_true(all(is.na(straight$binormal)))
  stopped <- vc_frenet_frame(function(t)c(t^2,0,0),0,
    function(t)c(2*t,0,0),function(t)c(2,0,0))
  expect_identical(stopped$status,"zero_speed")
  expect_true(all(is.na(unlist(stopped[c("tangent","normal","binormal")]))))
  expect_identical(vc_frenet_frame(function(t)c(t^3,0,0),0)$status,"unresolved_speed")
  parallel <- vc_frenet_frame(function(t)c(exp(t),2*exp(t),0),0,
    function(t)c(exp(t),2*exp(t),0),function(t)c(exp(t),2*exp(t),0))
  expect_identical(parallel$status,"unresolved_curvature")
  for (scale in c(1e-300,1e300)) {
    frame <- vc_frenet_frame(function(t)scale*c(t,t^2,0),0,
      function(t)scale*c(1,2*t,0),function(t)scale*c(0,2,0))
    expect_identical(frame$status,"ok")
    expect_equal(frame$normal,c(0,1,0))
  }
  edge <- function(t) {
    if (t<0) stop("outside curve domain")
    c(t,t^2,0)
  }
  expect_error(vc_frenet_frame(edge,0),"outside curve domain")
  expect_identical(vc_frenet_frame(edge,0,function(t)c(1,2*t,0),
    function(t)c(0,2,0))$status,"ok")
})

test_that("frame input validation rejects invalid steps, parameters and derivatives", {
  curve <- function(t)c(t,t^2,0)
  expect_error(vc_frenet_frame(1,0),"curve")
  expect_error(vc_frenet_frame(curve,NA_real_),"finite real scalar")
  expect_error(vc_frenet_frame(curve,0,h=0),"h")
  expect_error(vc_frenet_frame(curve,0,h2=-1),"positive")
  expect_error(vc_frenet_frame(curve,0,velocity=1),"velocity")
  expect_error(vc_frenet_frame(curve,0,acceleration=1),"acceleration")
  expect_error(vc_frenet_frame(curve,0,acceleration=function(t)c(0,2,Inf)),"three-vector")
  expect_error(vc_frenet_frame(function(t)c(t,t),0),"three-vector")
})

test_that("frame plots use stored coordinates and only draw selected available vectors", {
  skip_if_not_installed("plotly")
  calls <- 0L
  curve <- function(t) {
    calls <<- calls+1L
    c(t,t^2,t/2)
  }
  samples <- vc_curve_sample(curve,c(1,-1),n=6)
  frame <- vc_frenet_frame(curve,0)
  previous_calls <- calls
  plot <- vc_plot_frenet_frame(samples,frame,scale=.75,components=c("binormal","tangent"))
  built <- plotly::plotly_build(plot)
  expect_identical(calls,previous_calls)
  traces <- built$x$data
  expect_identical(vapply(traces,function(z)z$name,character(1)),
    c("Sampled curve","Frame point","binormal","tangent"))
  expect_equal(as.numeric(traces[[1]]$x),samples$samples$x)
  endpoint <- vapply(c("x","y","z"),function(axis)traces[[3]][[axis]][2],numeric(1))
  expect_equal(unname(endpoint),frame$point+.75*frame$binormal)
  expect_identical(built$x$layout$scene$aspectmode,"data")
  expect_error(vc_plot_frenet_frame(samples,frame,components="speed"),"components")
  expect_error(vc_plot_frenet_frame(samples,frame,scale=0),"scale")
  expect_error(vc_plot_frenet_frame(samples,vc_frenet_frame(curve,2)),"outside")

  line <- function(t)c(t,2*t,0)
  straight <- vc_plot_frenet_frame(vc_curve_sample(line,c(0,1)),vc_frenet_frame(line,.5))
  expect_identical(vapply(plotly::plotly_build(straight)$x$data,function(z)z$name,character(1)),
    c("Sampled curve","Frame point","tangent"))
})
