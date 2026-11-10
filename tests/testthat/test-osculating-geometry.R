test_that("osculating rings have the analytic helix center, radius and tangent", {
  curve <- function(t)c(cos(t),sin(t),t/2)
  v <- function(t)c(-sin(t),cos(t),.5)
  a <- function(t)c(-cos(t),-sin(t),0)
  result <- vc_osculating_circle(curve,c(pi/2,0,pi/2),v,a,n=16)
  expect_equal(result$t,c(pi/2,0,pi/2))
  circle <- result$circles[[1]]
  expect_identical(circle$status,"ok")
  expect_equal(circle$radius,1.25)
  expect_equal(circle$curvature,.8)
  expect_equal(circle$center,c(0,-.25,pi/4))
  expect_equal(circle$ring[1,],curve(pi/2))
  expect_identical(circle$ring[1,],circle$ring[17,])
  relative <- sweep(circle$ring,2,circle$center,"-")
  expect_equal(sqrt(rowSums(relative^2)),rep(circle$radius,17))
  expect_equal(as.numeric(relative%*%circle$frame$binormal),rep(0,17))
  expect_equal(circle$ring[5,],circle$center+circle$radius*circle$frame$tangent)
  expect_equal(result$circles[[1]],result$circles[[3]])
})

test_that("ribbon matrices preserve orientation, centers and unavailable sections", {
  curve <- function(t)c(t,t^3,0)
  v <- function(t)c(1,3*t^2,0)
  a <- function(t)c(0,6*t,0)
  ribbon <- vc_osculating_ribbon(curve,c(-1,1),v,a,n=5,n_u=4)
  expect_equal(dim(ribbon$grid$x),c(4,5))
  expect_equal(ribbon$u,c(0,1/3,2/3,1))
  expect_equal(ribbon$curve$samples$x,c(-1,-.5,0,.5,1))
  expect_identical(ribbon$sections[[3]]$status,"zero_curvature")
  expect_true(all(vapply(ribbon$grid,function(axis)all(is.na(axis[,3])),logical(1))))
  for (j in c(1,2,4,5)) {
    section <- ribbon$sections[[j]]
    expect_equal(vapply(ribbon$grid,function(axis)axis[1,j],numeric(1)),
      stats::setNames(section$frame$point,c("x","y","z")))
    expect_equal(unname(vapply(ribbon$grid,function(axis)axis[4,j],numeric(1))),
      section$center)
  }
  reversed <- vc_osculating_ribbon(curve,c(1,-1),v,a,n=5,n_u=4)
  expect_equal(reversed$grid$x,ribbon$grid$x[,5:1])
  partial <- vc_osculating_ribbon(curve,c(.5,.5),v,a,n_u=3,u_max=.4)
  expect_equal(dim(partial$grid$x),c(3,1))
  section <- partial$sections[[1]]
  expect_equal(unname(vapply(partial$grid,function(axis)axis[3,1],numeric(1))),
    section$frame$point+.4*(section$center-section$frame$point))
})

test_that("degenerate and unrepresentable circles retain their anchors", {
  line <- vc_osculating_circle(function(t)c(t,0,0),.5,
    function(t)c(1,0,0),function(t)c(0,0,0))$circles[[1]]
  expect_identical(line$status,"zero_curvature")
  expect_equal(line$curvature,0)
  expect_null(line$ring)
  expect_true(is.na(line$radius))
  expect_true(all(is.na(line$center)))
  expect_equal(line$frame$point,c(.5,0,0))
  stopped <- vc_osculating_circle(function(t)c(t^2,0,0),0,
    function(t)c(2*t,0,0),function(t)c(2,0,0))$circles[[1]]
  expect_identical(stopped$status,"zero_speed")
  expect_null(stopped$ring)
  shifted <- vc_osculating_circle(function(t)c(1e20+cos(t),sin(t),0),0,
    function(t)c(-sin(t),cos(t),0),function(t)c(-cos(t),-sin(t),0))$circles[[1]]
  expect_identical(shifted$status,"unrepresentable")
  expect_null(shifted$ring)
  expect_error(vc_osculating_circle(identity,c(0,NA)),"t")
  expect_error(vc_osculating_circle(identity,0,n=2),"n")
  expect_error(vc_osculating_ribbon(identity,c(0,1),u_max=1.1),"u_max")
  expect_error(vc_osculating_ribbon(identity,c(0,1),n_u=1),"n_u")
})

test_that("osculating displays use stored geometry and leave gaps open", {
  skip_if_not_installed("plotly")
  calls <- 0L
  curve <- function(t) {
    calls <<- calls+1L
    c(t,t^3,0)
  }
  v <- function(t)c(1,3*t^2,0)
  a <- function(t)c(0,6*t,0)
  ribbon <- vc_osculating_ribbon(curve,c(-1,1),v,a,n=5,n_u=3)
  circles <- vc_osculating_circle(curve,c(.5,0),v,a,n=12)
  before <- calls
  disk <- plotly::plotly_build(vc_plot_osculating_circle(circles,ribbon$curve,
    fill="disk",show_radius=TRUE))
  expect_identical(vapply(disk$x$data,function(trace)trace$type,character(1)),
    c("scatter3d","scatter3d","mesh3d","scatter3d"))
  expect_equal(as.integer(disk$x$data[[3]]$j),1:12)
  expect_equal(as.integer(disk$x$data[[3]]$k),2:13)
  expect_equal(as.numeric(disk$x$data[[4]]$x),
    c(circles$circles[[1]]$center[1],circles$circles[[1]]$frame$point[1]))
  rings <- plotly::plotly_build(vc_plot_osculating_circle(circles))
  expect_equal(as.numeric(rings$x$data[[2]]$x),circles$circles[[1]]$ring[,1])
  expect_error(vc_plot_osculating_circle(circles,
    vc_curve_sample(function(t)c(t,t^3,0),c(.7,1))),"outside")

  built <- plotly::plotly_build(vc_plot_osculating_ribbon(ribbon))
  expect_identical(vapply(built$x$data,function(trace)trace$type,character(1)),
    c("scatter3d","surface","scatter3d"))
  expect_true(all(is.na(built$x$data[[2]]$x[,3])))
  expect_identical(built$x$data[[2]]$connectgaps,FALSE)
  expect_true(is.na(built$x$data[[3]]$x[3]))
  expect_equal(as.numeric(built$x$data[[1]]$x),ribbon$curve$samples$x)
  expect_identical(calls,before)
  one <- vc_osculating_ribbon(curve,c(.5,.5),v,a,n_u=3)
  isolated <- plotly::plotly_build(vc_plot_osculating_ribbon(one,show_centers=FALSE))
  expect_identical(isolated$x$data[[2]]$name,"Isolated sections")
  expect_length(isolated$x$data,2)
  line <- vc_osculating_ribbon(function(t)c(t,0,0),c(0,1),n=3)
  expect_length(plotly::plotly_build(vc_plot_osculating_ribbon(line))$x$data,1)
})
