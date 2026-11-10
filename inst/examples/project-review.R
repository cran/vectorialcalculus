# Open next.Rproj and Source this file. It is the single review entry point.
# The original source is read, never modified or attached as a package.
run_project_review <- function(export_dir="project-review",legacy_dir="..") {
  if(!file.exists("R/secant_convergence.R")) stop("Open next.Rproj first.")
  if(!requireNamespace("plotly",quietly=TRUE)) stop("Install plotly first.")
  env <- new.env(parent=baseenv())
  for(file in list.files("R",pattern="[.]R$",full.names=TRUE)) sys.source(file,env)
  # Reuse approved example scripts without their automatic Viewer/export calls.
  review <- new.env(parent=baseenv())
  review$interactive <- function()FALSE
  sys.source("inst/examples/visual-validation.R",review)
  sys.source("inst/examples/integral-review.R",review)
  sys.source("inst/examples/multivariable-review.R",review)
  plots <- c(stats::setNames(review$visual_review$plots,paste0("phase1-",names(review$visual_review$plots))),
             stats::setNames(review$integral_review$plots,paste0("phase2-",names(review$integral_review$plots))),
             stats::setNames(review$multivariable_review$plots,paste0("phase3-",names(review$multivariable_review$plots))))
  helix <- function(t)c(cos(t),sin(t),t/2)
  helix_forward <- env$vc_curve_sample(helix,c(0,2*pi),n=65)
  helix_reverse <- env$vc_curve_sample(helix,c(2*pi,0),n=65)
  helix_length <- env$vc_curve_length(helix,c(0,2*pi),
    velocity=function(t)c(-sin(t),cos(t),.5))
  stopifnot(abs(helix_length$estimate-pi*sqrt(5))<1e-9)
  plots[["phase4-helix"]] <- env$vc_plot_curve(helix_forward,title="Helix in increasing parameter order")
  plots[["phase4-helix_reverse"]] <- env$vc_plot_curve(helix_reverse,"lines+markers",
    title="Helix in decreasing parameter order")
  frame <- env$vc_frenet_frame(helix,pi/4)
  stopifnot(frame$status=="ok")
  plots[["phase4-helix_frenet"]] <- env$vc_plot_frenet_frame(helix_forward,frame,
    scale=.8,title="Right-handed Frenet frame on a helix")
  line <- function(t)c(t,2*t,0)
  straight <- env$vc_frenet_frame(line,.5)
  stopifnot(straight$status=="unresolved_curvature")
  plots[["phase4-line_tangent"]] <- env$vc_plot_frenet_frame(
    env$vc_curve_sample(line,c(0,1)),straight,scale=.35,
    title="Straight line with tangent only")
  measures <- env$vc_curvature_torsion(helix,pi/4,h3=.004)
  stopifnot(abs(measures$curvature-.8)<1e-6,abs(measures$torsion-.4)<1e-5)
  plots[["phase4-helix_curvature_torsion"]] <- env$vc_plot_curvature_torsion(
    helix_forward,measures,scale=.8,title="Helix curvature and signed torsion")
  circles <- env$vc_osculating_circle(helix,c(pi/2,3*pi/2))
  plots[["phase4-helix_osculating_circles"]] <- env$vc_plot_osculating_circle(
    circles,helix_forward,fill="disk",show_radius=TRUE,title="Osculating disks on a helix")
  inflection <- function(t)c(t,t^3,t/2)
  ribbon <- env$vc_osculating_ribbon(inflection,c(-.75,.75),n=25,u_max=.3,
    velocity=function(t)c(1,3*t^2,.5),acceleration=function(t)c(0,6*t,0))
  stopifnot(ribbon$sections[[13]]$status=="zero_curvature")
  plots[["phase4-ribbon_inflection"]] <- env$vc_plot_osculating_ribbon(
    ribbon,show_centers=FALSE,title="Osculating ribbon with a gap at the inflection")
  rotation <- function(p)c(-p[2],p[1],0)
  field_lower <- function(x)-1
  field_upper <- function(x)1
  field_bottom <- function(p)-.5
  field_top <- function(p).5
  field_grid <- env$vc_vector_field_sample(rotation,field_lower,field_upper,
    field_bottom,field_top,c(-1,1),n=4,n_u=4,n_v=2)
  streamline <- env$vc_streamline(rotation,c(1,0,0),c(0,1),n=16)
  plots[["phase4-rotating_field"]] <- env$vc_plot_vector_field(field_grid,streamline,
    title="Rotating field and streamline")
  planar_line <- function(t)c(t,2*t)
  scalar_line <- env$vc_line_integral(function(p)p[1],planar_line,c(-1,1),
    method="simpson",n=16,velocity=function(t)c(1,2))
  work_curve <- function(t)c(cos(t),sin(t))
  work_field <- function(p)c(-p[2],p[1])
  planar_work <- env$vc_line_integral(work_field,work_curve,c(0,2*pi),"work",
    method="simpson",n=24,velocity=function(t)c(-sin(t),cos(t)))
  plots[["phase4-scalar_line"]] <- env$vc_plot_line_integral(scalar_line,
    title="Signed scalar curtain along a planar path")
  plots[["phase4-planar_work"]] <- env$vc_plot_line_integral(planar_work,
    title="Positive circulation around a counterclockwise circle")
  ruled <- env$vc_cylindrical_surface(helix,c(0,pi),c(-.5,1),c(0,0,2),n_t=19,n_s=9)
  plane <- function(p)c(p[1],p[2],2*p[1]+3*p[2])
  patch <- env$vc_surface_sample(plane,c(0,1),c(0,2),n_u=15,n_v=11)
  plots[["phase4-ruled_surface"]] <- env$vc_plot_surface(ruled,title="Ruled surface with a helix directrix")
  plots[["phase4-surface_graph"]] <- env$vc_plot_surface(patch,title="Graph patch for surface integration")
  green_lower <- function(x)0
  green_upper <- function(x)x^2
  green <- env$vc_green(function(p)c(0,p[1]^2),green_lower,green_upper,c(0,1),
    lower_derivative=function(x)0,upper_derivative=function(x)2*x,
    rotation=function(p)2*p[1],n_boundary=8,n_outer=4,n_inner=4)
  stopifnot(green$within_tolerance,abs(green$circulation-.5)<1e-8)
  plots[["phase4-green"]] <- env$vc_plot_green(green,
    env$vc_region(green_lower,green_upper,c(0,1)),title="Green theorem on a curved strip")
  stokes_patch <- function(p)c(p,sum(p^2))
  stokes <- env$vc_stokes(function(p)c(-p[2]*p[3],p[1]*p[3],0),stokes_patch,
    c(0,1),c(0,2),jacobian=function(p)rbind(diag(2),2*p),
    curl=function(p)c(-p[1],-p[2],2*p[3]),n_boundary=8,n_outer=4,n_inner=4)
  stopifnot(stokes$within_tolerance,abs(stokes$circulation-40/3)<1e-8)
  plots[["phase4-stokes"]] <- env$vc_plot_stokes(stokes,
    env$vc_surface_sample(stokes_patch,c(0,1),c(0,2)),title="Stokes theorem on a curved patch")
  gauss <- env$vc_gauss(identity,c(0,1),c(0,1),c(0,1),
    n_outer=2,n_middle=2,n_inner=2)
  stopifnot(gauss$within_tolerance,abs(gauss$outward_flux-3)<1e-8)
  plots[["phase4-gauss"]] <- env$vc_plot_gauss(gauss,title="Gauss theorem on the unit cube")
  f <- function(x)x^2
  steps <- 2^(-(1:8))
  forward <- env$vc_secant_convergence(f,1,steps=steps)
  central <- env$vc_secant_convergence(f,1,steps=steps,method="central")
  plots$recovered_secant_forward <- env$vc_plot_secant_convergence(forward,f)
  plots$recovered_secant_central <- env$vc_plot_secant_convergence(central,f)
  differential <- env$vc_differential(f,1,.2,warn=FALSE)
  taylor <- env$vc_taylor(f,1,order=2,warn=FALSE)
  lesson <- data.frame(capability=c("Derivative at 1","Secant from 1 to 1.2","Differential dx=0.2",
    "Actual finite change","Taylor order 2 at 1.2","Integral from 0 to 2","Accumulation at 2, lower 0","Revolution volume on [0,1]"),
    value=c(env$vc_derivative(f,1,warn=FALSE)$estimate,env$vc_secant(f,1,1.2)$slope,
      differential$dy,differential$delta_y,taylor$polynomial(1.2),env$vc_integrate(f,c(0,2),method="simpson")$estimate,
      env$vc_accumulation(f,2)$table$estimate,env$vc_revolution(f,c(0,1))$estimate),
    analytical=c(2,2.2,.4,.44,1.44,8/3,8/3,pi/5))
  stopifnot(all(abs(lesson$value-lesson$analytical)<1e-7))
  for(name in c("recovered_secant_forward","recovered_secant_central"))
    stopifnot(length(plotly::plotly_build(plots[[name]])$x$frames)==length(steps))
  comparisons <- data.frame(original=character(),capability=character(),legacy=numeric(),
    rebuilt=numeric(),difference=numeric())
  findings <- data.frame(original=character(),issue=character(),observed_error=numeric())
  compare <- function(original,capability,a,b) {
    stopifnot(is.finite(a),is.finite(b),abs(a-b)<1e-8)
    comparisons <<- rbind(comparisons,data.frame(original=original,capability=capability,
      legacy=a,rebuilt=b,difference=b-a))
  }
  legacy_available <- file.exists(file.path(legacy_dir,"R/secant_tangent.R"))
  map <- data.frame(original=character(),phase=character(),recovery=character())
  if(legacy_available) {
    old <- new.env(parent=baseenv())
    # Suppress legacy automatic print only in this review environment.
    old$print <- function(x,...)invisible(x)
    for(file in c("secant_tangent.R","riemann_rectangles2d.R","riemann_sum_1d_plot.R","solid_of_revolution_y.R",
                 "gradient_scalar.R","directional_derivative3d.R","gradient_direction2d.R",
                 "total_differential_nd.R","tangent_plane3d.R","plot_surface_with_tangents.R","critical_points_nd.R",
                 "partial_derivative_surface.R","newton_raphson2d.R","lagrange_check.R","integrate_double_xy.R","integrate_double_polar.R","integrate_triple_general.R","solid_cylindrical3d.R","related_rates_grad.R","riemann_sum_2d_plot.R","riemann_prisms3d.R","xy_region.R","region_xyz0.R","solid_xyz3d.R","solid_spherical3d.R","critical_points_2d.R","curve_sample3d.R","plot_curve3d.R","arc_length3d.R",
                 "tangent3d.R","normal3d.R","binormal3d.R","frenet_frame3d.R","curvature_torsion3d.R",
                 "osculating_circle3d.R","osculating_ribbon3d.R",
                 "vector_field3d.R","streamline_and_field3d.R","divergence_field.R","curl3d.R",
                 "line_integral2d.R","line_integral_vector2d.R","line_integral3d_work.R",
                 "cylindrical_surface3d.R","surface_integral_z.R","surface_parametric_area.R"))
      sys.source(file.path(legacy_dir,"R",file),old)
    old_curve <- old$curve_sample3d(cos,sin,function(t)t/2,2*pi,0,n_samples=65)
    old_ruled <- old$cylindrical_surface3d(cos,sin,function(t)t/2,0,pi,c(-.5,1),
      c(0,0,2),n_t=19,n_s=9,plot=FALSE)
    for (axis in c("x","y","z")) {
      m <- old_ruled[[paste0(toupper(axis),"mat")]]
      stopifnot(isTRUE(all.equal(m,ruled$grid[[axis]])))
      compare("cylindrical_surface3d",paste("ruled grid",axis),m[4,7],ruled$grid[[axis]][4,7])
    }
    g <- function(x,y)2*x+3*y
    old_graph <- old$surface_integral_z(g,c(0,1),c(0,2),nx=9,ny=7,
      scalar_phi=function(x,y,z)-x,vector_F=function(x,y,z)c(1,2,3),plot=FALSE)
    new_scalar_graph <- env$vc_surface_graph_integral(function(p)g(p[1],p[2]),
      c(0,1),c(0,2),function(p)-p[1],n_outer=2,n_inner=2)
    new_flux <- env$vc_surface_graph_integral(function(p)g(p[1],p[2]),
      c(0,1),c(0,2),function(p)c(1,2,3),"flux",n_outer=2,n_inner=2)
    compare("surface_integral_z","plane signed scalar integral",
      old_graph$area_density_integral,new_scalar_graph$estimate)
    compare("surface_integral_z","plane upward flux",old_graph$flux_integral,new_flux$estimate)
    old_area <- old$surface_parametric_area(function(u,v)u,function(u,v)v,g,
      c(0,1),c(0,2),nu=9,nv=7)
    new_area <- env$vc_surface_integral(plane,c(0,1),c(0,2),n_outer=2,n_inner=2)
    compare("surface_parametric_area","affine patch area",old_area$area,new_area$estimate)
    old_scalar <- old$line_integral2d(function(x,y)x,planar_line,0,1,plot=FALSE)
    new_scalar <- env$vc_line_integral(function(p)p[1],planar_line,c(0,1),
      velocity=function(t)c(1,2))
    compare("line_integral2d","weighted planar segment",old_scalar$value,new_scalar$estimate)
    old_planar <- old$line_integral_vector2d(function(x,y)c(-y,x),work_curve,0,2*pi,
      h=1e-5,plot=FALSE)
    compare("line_integral_vector2d","unit circle circulation",old_planar$value,planar_work$estimate)
    spatial <- function(t)c(t,2*t,3*t)
    constant <- function(x,y,z)c(0,0,1)
    new_spatial <- env$vc_line_integral(function(p)c(0,0,1),spatial,c(0,1),"work",
      velocity=function(t)c(1,2,3))
    odd <- old$line_integral3d_work(constant,spatial,0,1,n_curve=101,plot=FALSE)
    even <- old$line_integral3d_work(constant,spatial,0,1,n_curve=100,plot=FALSE)
    compare("line_integral3d_work","constant work, complete Simpson grid",odd$value,new_spatial$estimate)
    stopifnot(abs(even$value-3*98/99)<1e-8)
    findings <- rbind(findings,data.frame(original="line_integral3d_work",
      issue="With an even sample count, the original drops the last path interval before Simpson integration. The reconstruction integrates the full parameter interval.",
      observed_error=abs(even$value-new_spatial$estimate)))
    fxyz <- function(x,y,z)c(-y,x,0)
    zlo <- function(x,y)-.5
    zhi <- function(x,y).5
    old_field <- old$vector_field3d(fxyz,field_lower,field_upper,zlo,zhi,-1,1,
      NX=4,NY=4,NZ=2,plot=FALSE)
    stopifnot(isTRUE(all.equal(as.matrix(old_field$points[c("x","y","z","mag")]),
      as.matrix(field_grid$samples[c("x","y","z","magnitude")]),check.attributes=FALSE)))
    compare("vector_field3d","sample magnitude",old_field$points$mag[11],field_grid$samples$magnitude[11])
    old_path <- old$streamline_and_field3d(fxyz,field_lower,field_upper,zlo,zhi,-1,1,
      NX=4,NY=4,NZ=2,p=c(1,0,0),t_final=1,step=1/16,arrows="none")
    stopifnot(isTRUE(all.equal(old_path$traj,streamline$samples,check.attributes=FALSE)))
    compare("streamline_and_field3d","final y",old_path$traj$y[17],streamline$samples$y[17])
    growth <- function(p)c(2*p[1],-p[2],3*p[3])
    compare("divergence_field","linear field trace",old$divergence_field(growth,c(1,2,3),h=.01),
      env$vc_divergence(growth,c(1,2,3),h=.01)$estimate)
    compare("curl3d","rotation z component",old$curl3d(fxyz,1,2,3,h=.01)[3],
      env$vc_curl(rotation,c(1,2,3),h=.01)$estimate[3])
    for(name in c("t","x","z"))
      compare("curve_sample3d",paste("reversed helix",name),old_curve[[name]][17],
        helix_reverse$samples[[name]][17])
    old_plot <- plotly::plotly_build(old$plot_curve3d(old_curve))$x$data[[1]]
    new_plot <- plotly::plotly_build(plots[["phase4-helix_reverse"]])$x$data[[1]]
    stopifnot(all(vapply(c("x","y","z"),function(axis)
      isTRUE(all.equal(as.numeric(old_plot[[axis]]),as.numeric(new_plot[[axis]]))),logical(1))))
    compare("plot_curve3d","reversed helix plotted y",old_plot$y[17],new_plot$y[17])
    old_length <- old$arc_length3d(cos,sin,function(t)t/2,0,2*pi,method_int="integrate")
    compare("arc_length3d","helix arc length",old_length,helix_length$estimate)
    cubic <- function(t)c(t,t^2/2,t^3/6)
    step <- 2^-8
    old_measures <- old$curvature_torsion3d(identity,function(t)t^2/2,
      function(t)t^3/6,0,h=step)
    new_measures <- env$vc_curvature_torsion(cubic,0,h=step,h2=step,h3=step)
    compare("curvature_torsion3d","cubic curvature",old_measures$kappa,new_measures$curvature)
    stopifnot(abs(old_measures$tau+1)<1e-8,abs(new_measures$torsion-1)<1e-8)
    findings <- rbind(findings,data.frame(original="curvature_torsion3d",
      issue="Original third-difference stencil reverses signed torsion. The cubic (t, t^2/2, t^3/6) gives -1 in the original and +1 in the reconstruction, matching the analytical value.",
      observed_error=abs(old_measures$tau-1)))
    parabola <- function(t)c(t,t^2,t/2)
    parabola_y <- function(t)t^2
    parabola_z <- function(t)t/2
    step <- 2^-8
    old_circle <- old$osculating_circle3d(identity,parabola_y,parabola_z,-.5,.5,0,h=step)$data
    new_circle <- env$vc_osculating_circle(parabola,0,h=step,h2=step)$circles[[1]]
    compare("osculating_circle3d","parabola radius",old_circle$radius,new_circle$radius)
    for (i in seq_len(3)) {
      axis <- c("x","y","z")[i]
      compare("osculating_circle3d",paste("parabola center",axis),
        old_circle[[paste0("c",axis)]],new_circle$center[i])
    }
    old_ribbon <- old$osculating_ribbon3d(identity,parabola_y,parabola_z,-.5,.5,
      h=step,n_t=5,n_u=3,u_max=.75,plot=TRUE)
    old_grid <- plotly::plotly_build(old_ribbon$plot)$x$data[[1]]
    new_ribbon <- env$vc_osculating_ribbon(parabola,c(-.5,.5),h=step,h2=step,
      n=5,n_u=3,u_max=.75)
    for (axis in c("x","y","z")) {
      stopifnot(isTRUE(all.equal(old_grid[[axis]],new_ribbon$grid[[axis]],tolerance=1e-8,
        check.attributes=FALSE)))
      compare("osculating_ribbon3d",paste("parabola ribbon",axis),
        old_grid[[axis]][3,4],new_ribbon$grid[[axis]][3,4])
    }
    legacy_gap <- old$osculating_ribbon3d(identity,function(t)t^3,parabola_z,-.75,.75,
      h=step,n_t=25,n_u=3)$data
    stopifnot(is.na(legacy_gap$x[13]),ribbon$curve$samples$x[13]==0,
      all(is.na(ribbon$grid$x[,13])))
    findings <- rbind(findings,data.frame(original="osculating_ribbon3d",
      issue="Original omits the curve point at an inflection along with its undefined ribbon section. The reconstruction retains the curve point and leaves only the ribbon section unavailable.",
      observed_error=NA_real_))
    step <- 2^-10
    new_frame <- env$vc_frenet_frame(parabola,0,h=step,h2=step)
    legacy_frame <- old$frenet_frame3d(identity,function(t)t^2,function(t)t/2,-1,1,0,h=step)
    for (component in c("tangent","normal","binormal")) {
      name <- paste0(component,"3d")
      single <- old[[name]](identity,function(t)t^2,function(t)t/2,-1,1,0,h=step)
      axes <- paste0(toupper(substr(component,1,1)),c("x","y","z"))
      for (i in seq_along(axes)) {
        column <- axes[i]
        compare(name,paste("parabola",column),single[[column]],new_frame[[component]][i])
        compare("frenet_frame3d",paste("parabola frame",column),
          legacy_frame[[column]],new_frame[[component]][i])
      }
    }
    sec <- old$secant_tangent(f,1,h_vals=steps,method="forward")
    compare("secant_tangent","finest forward slope",sec$derivative,forward$table$slope[length(steps)])
    plots$legacy_secant <- sec$plot
    plots$legacy_riemann_animation <- old$riemann_rectangles2d(f,0,2,n_vals=c(4,8,16,32),y0=0)
    sums <- old$riemann_sum_1d_plot(f,c(0,2),n=32)
    for(rule in c("left","right","midpoint")) {
      value <- switch(rule,left=sums$lower_sum,right=sums$upper_sum,midpoint=sums$mid_sum)
      compare("riemann_sum_1d_plot",paste("x^2 monotone",rule),value,env$vc_riemann(f,c(0,2),n=32,method=rule)$estimate)
    }
    plots$legacy_riemann_static <- sums$figure
    shifted <- function(x)2+x
    solid <- old$solid_of_revolution_y(shifted,c(0,1),a=2,deriv=function(x)1,include_end_caps=TRUE,plot=FALSE)
    rebuilt <- env$vc_revolution(shifted,c(0,1),g=function(x)2,axis_at=2)
    surface <- env$vc_surface_area(shifted,c(0,1),axis_at=2,df=function(x)1,end_caps=TRUE)
    compare("solid_of_revolution_y","shifted cone volume",solid$volume,rebuilt$estimate)
    compare("solid_of_revolution_y","lateral area",solid$surface_area_lateral,surface$lateral_area)
    compare("solid_of_revolution_y","total area with end disks",solid$surface_area_total,surface$estimate)
    plots$recovered_shifted_solid <- env$vc_plot_revolution(rebuilt,shifted,g=function(x)2,sections=c(.25,.5,.75))
    mv <- review$multivariable_review
    legacy_gradient <- old$gradient_scalar(mv$field,mv$point,plot=FALSE)
    old_partials <- old$partial_derivatives_surface(function(x,y)x^2+2*y^2,1,1,plot=FALSE)
    compare("partial_derivatives_surface","coordinate partial x",old_partials$fx,unname(mv$gradient$estimate[1]))
    compare("partial_derivatives_surface","coordinate partial y",old_partials$fy,unname(mv$gradient$estimate[2]))
    compare("partial_derivatives_surface","anchor value",old_partials$f0,mv$gradient$value)
    for(i in 1:2) compare("gradient_scalar",paste("component",i),legacy_gradient[i],unname(mv$gradient$estimate[i]))
    legacy_direction <- old$directional_derivative3d(mv$field,mv$point,c(3,4),plot=FALSE)
    compare("directional_derivative3d","direct unit-direction derivative",legacy_direction$D,mv$directional$estimate)
    plots$legacy_gradient_animation <- old$gradient_direction2d(function(x,y)x^2+2*y^2,
      x0=1,y0=1,xlim=c(-1,3),ylim=c(-1,3),n_grid=21,theta_vals=seq(0,2*pi,length.out=13))
    linear <- mv$linearization
    old_total <- old$total_differential_nd(mv$field,mv$point)
    compare("total_differential_nd","requested differential",old_total$differential(c(.1,.2)),linear$differential(c(.1,.2)))
    compare("total_differential_nd","affine local value",old_total$value+old_total$differential(c(.1,.2)),linear$model(c(1.1,1.2)))
    old_plane <- old$tangent_plane3d(function(x,y)x^2+2*y^2,mv$point,plot=TRUE,grid=21,
      x_window=2,y_window=2,z_window=6)
    compare("tangent_plane3d","plane x coefficient",old_plane$fx,unname(linear$gradient$estimate[1]))
    compare("tangent_plane3d","plane y coefficient",old_plane$fy,unname(linear$gradient$estimate[2]))
    compare("tangent_plane3d","plane value",old_plane$plane_fun(1.1,1.2),linear$model(c(1.1,1.2)))
    compare("tangent_plane3d","implicit offset",unname(old_plane$plane_coeff[4]),unname(linear$plane_coefficients[4]))
    plots$legacy_tangent_plane <- old_plane$fig
    old_surface <- old$plot_surface_with_tangents(function(x,y)x^2+2*y^2,
      x0=1,y0=1,xlim=c(-1,3),ylim=c(-1,3),n=21)
    old_built <- plotly::plotly_build(old_surface)
    for(i in 1:2) {
      line <- Filter(function(z)identical(z$name,if(i==1)"Tangent along x" else "Tangent along y"),old_built$x$data)[[1]]
      u <- if(i==1)as.numeric(line$x) else as.numeric(line$y)
      slope <- (line$z[length(line$z)]-line$z[1])/(u[length(u)]-u[1])
      compare("plot_surface_with_tangents",paste("coordinate tangent",i),slope,unname(linear$gradient$estimate[i]))
    }
    original_surface <- Filter(function(z)identical(z$type,"surface"),old_built$x$data)[[1]]
    expected_surface <- outer(original_surface$y,original_surface$x,function(y,x)x^2+2*y^2)
    grid_error <- max(abs(original_surface$z-expected_surface))
    stopifnot(grid_error>1)
    findings <- rbind(findings,data.frame(original="plot_surface_with_tangents",
      issue="Legacy expand.grid/matrix order transposes surface coordinates; the reconstructed mesh uses y rows and x columns.",observed_error=grid_error))
    plots$legacy_coordinate_tangents <- old_surface
    # Evaluate the unaltered private Hessian definition, not its optimizer.
    definitions <- Filter(function(expr)is.call(expr)&&identical(expr[[1]],as.name("<-"))&&
      identical(expr[[2]],as.name(".hess")),as.list(body(old$critical_points_nd)))
    stopifnot(length(definitions)==1)
    henv <- new.env(parent=baseenv())
    henv$n <- 2L
    henv$f <- mv$hessian_field
    henv$.steps <- function(x)rep(.0025,length(x)) # match the finest new step
    eval(definitions[[1]],henv)
    old_hessian <- henv$.hess(mv$point)
    for(ij in list(c(1,1),c(1,2),c(2,2))) compare("critical_points_nd private .hess",
      paste("direct Hessian entry",paste(ij,collapse=",")),old_hessian[ij[1],ij[2]],mv$hessian$estimate[ij[1],ij[2]])
    old_newton <- old$newton_raphson2d(function(x,y)mv$newton_field(c(x,y)),c(.8,.2),
      J=function(x,y)mv$newton_jacobian(c(x,y)),tol=1e-10,n_grid=51,
      xlim=c(-1,1),ylim=c(-1,1))
    plots$legacy_newton_system <- old_newton$plot
    for(i in 1:2) {
      compare("newton_raphson2d",paste("terminal coordinate",i),old_newton$root[i],mv$newton_reference$estimate[i])
      compare("newton_raphson2d",paste("first undamped update",i),old_newton$iterates[2,c("x","y")][[i]],
        mv$newton_reference$points[2,i])
    }
    old_contours <- plotly::plotly_build(old_newton$plot)$x$data
    grid <- seq(-1,1,length.out=51)
    contour_error <- vapply(1:2,function(i) {
      expected <- outer(grid,grid,Vectorize(function(y,x)mv$newton_field(c(x,y))[i]))
      max(abs(matrix(as.numeric(old_contours[[i]]$z),51)-expected))
    },numeric(1))
    stopifnot(max(abs(contour_error-c(1,6)))<1e-10)
    findings <- rbind(findings,data.frame(original="newton_raphson2d",
      issue="Original zero-contour matrices use x rows instead of y rows: asymmetric test component discrepancies 1 and 6; new view corrects orientation.",
      observed_error=max(contour_error)))
    old_stationary_surface <- (function() {
      had_seed <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
      if(had_seed)seed <- get(".Random.seed",envir=.GlobalEnv)
      on.exit(if(had_seed)assign(".Random.seed",seed,envir=.GlobalEnv) else
        if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv))
      old$critical_points_2d(function(x,y)x^2+2*y^2,c(-1,2),c(-2,3),
        start_n=c(1,1),n_rand=0,h=.01,maxit=50,plot=TRUE,grid_plot=c(5,5))
    })()
    old_graph <- plotly::plotly_build(old_stationary_surface$fig)$x$data[[1]]
    new_graph <- plotly::plotly_build(plots[["phase3-stationary_surface_minimum"]])$x$data[[1]]
    stopifnot(max(abs(as.numeric(old_graph$z)-as.numeric(new_graph$z)))<1e-12)
    compare("critical_points_2d","3D surface asymmetric grid height A",old_graph$z[2,4],new_graph$z[2,4])
    compare("critical_points_2d","3D surface asymmetric grid height B",old_graph$z[4,2],new_graph$z[4,2])
    old_sph <- old$solid_spherical3d(R1=function(t,p).5+.1*cos(t)*sin(p),
      R2=function(t,p)1+.2*sin(t)*sin(p),theta_range=c(0,pi/2),phi_range=c(0,pi),
      n_theta=5,n_phi=4,plot=FALSE)
    for(pair in list(c("R1_surf","radial_inner"),c("R2_surf","radial_outer")))
      for(axis in c("x","y","z")) {
        a <- old_sph[[pair[1]]][[toupper(axis)]]
        b <- mv$solid_sph_bridge$faces[[pair[2]]][[axis]]
        stopifnot(identical(dim(a),dim(b)),max(abs(as.numeric(a)-as.numeric(b)))<1e-12)
      }
    compare("solid_spherical3d","inner shell sampled x",old_sph$R1_surf$X[2,2],mv$solid_sph_bridge$faces$radial_inner$x[2,2])
    compare("solid_spherical3d","inner shell sampled z",old_sph$R1_surf$Z[3,4],mv$solid_sph_bridge$faces$radial_inner$z[3,4])
    compare("solid_spherical3d","outer shell sampled y",old_sph$R2_surf$Y[2,4],mv$solid_sph_bridge$faces$radial_outer$y[2,4])
    compare("solid_spherical3d","outer shell sampled z",old_sph$R2_surf$Z[3,3],mv$solid_sph_bridge$faces$radial_outer$z[3,3])
    old_cyl_faces <- old$solid_cylindrical3d(R1=function(t).5+.1*t,R2=function(t)1+.2*t,
      Z1=function(r,t)t-.2*r,Z2=function(r,t)2+t+.3*r,th_min=0,th_max=pi/2,
      plot=TRUE,n_theta=5,n_u=4,n_v=3,mode="faces",show_surface_grid=FALSE)
    cyl_faces <- Filter(function(z)identical(z$type,"surface"),plotly::plotly_build(old_cyl_faces$fig)$x$data)
    stopifnot(length(cyl_faces)==6)
    for(i in seq_len(6))for(axis in c("x","y","z")) {
      target <- mv$solid_cyl_sector$faces[[i]][[axis]]
      if(i<=2)target <- t(target) # Original theta patches use u rows/v columns.
      stopifnot(identical(dim(cyl_faces[[i]][[axis]]),dim(target)),
        max(abs(as.numeric(cyl_faces[[i]][[axis]])-as.numeric(target)))<1e-12)
    }
    compare("solid_cylindrical3d","theta patch sampled x",cyl_faces[[1]]$x[2,1],mv$solid_cyl_sector$faces$theta_min$x[1,2])
    compare("solid_cylindrical3d","outer wall sampled z",cyl_faces[[4]]$z[2,3],mv$solid_cyl_sector$faces$radial_outer$z[2,3])
    compare("solid_cylindrical3d","bottom sampled x",cyl_faces[[5]]$x[3,4],mv$solid_cyl_sector$faces$z_lower$x[3,4])
    compare("solid_cylindrical3d","top sampled y",cyl_faces[[6]]$y[2,3],mv$solid_cyl_sector$faces$z_upper$y[2,3])
    old_box <- old$solid_xyz3d(H1=function(x)0,H2=function(x)1,G1=function(x,y)0,G2=function(x,y)2,
      a=0,b=3,plot=FALSE,n_x=3,n_u=3,n_v=4,compute_volume=TRUE,vol_method="grid",nx_vol=5,ny_vol=5)
    old_wedge <- old$solid_xyz3d(H1=function(x)0,H2=function(x)1-x,G1=function(x,y)0,G2=function(x,y)1,
      a=0,b=1,plot=FALSE,n_x=3,n_u=3,n_v=4,compute_volume=TRUE,vol_method="grid",nx_vol=5,ny_vol=5)
    compare("solid_xyz3d","box volume via existing triple engine",old_box$volume$estimate,mv$solid_box_volume$estimate)
    compare("solid_xyz3d","wedge volume via existing triple engine",old_wedge$volume$estimate,mv$solid_wedge_volume$estimate)
    old_face <- old$solid_xyz3d(H1=function(x)0,H2=function(x)1,G1=function(x,y)0,G2=function(x,y)2,
      a=0,b=3,plot=TRUE,n_x=3,n_u=3,n_v=4,show_faces=c(FALSE,FALSE,FALSE,FALSE,FALSE,TRUE),
      show_edges=FALSE,show_surface_grid=FALSE,mode="faces")
    face <- plotly::plotly_build(old_face$fig)$x$data[[1]]
    for(axis in c("x","y","z"))stopifnot(max(abs(as.numeric(face[[axis]])-
      as.numeric(mv$solid_box$faces$z_upper[[axis]])))<1e-12)
    compare("solid_xyz3d","upper-face sampled x",face$x[1,2],mv$solid_box$faces$z_upper$x[1,2])
    compare("solid_xyz3d","upper-face sampled z",face$z[2,3],mv$solid_box$faces$z_upper$z[2,3])
    old_xy <- old$xy_region(H1=function(x)x^2,H2=function(x)2+x,a=0,b=2,D=4,
      plot=FALSE,n_curve=5,n_u=4)
    old_plate <- old$region_xyz0(H1=function(x)x^2,H2=function(x)2+x,a=0,b=2,z0=2,D=4,
      plot=FALSE,n_curve=5)
    compare("xy_region","sampled lower boundary",old_xy$y_low[2],mv$region_bridge$sections$lower[2])
    compare("xy_region","sampled upper boundary",old_xy$y_high[4],mv$region_bridge$sections$upper[4])
    compare("region_xyz0","plate boundary data",old_plate$y_high[2],mv$region_bridge$sections$upper[2])
    compare("region_xyz0","strip break",old_plate$x_part[4],mv$region_bridge$partition$breaks[4])
    old_riemann2d <- old$riemann_sum_2d_plot(function(x,y)x+2*y,c(0,2),c(0,1),nx=2,ny=1,
      show_surface=FALSE,base_plane=FALSE,show_base_grid=FALSE)
    compare("riemann_sum_2d_plot","rectangular midpoint",old_riemann2d$mid_sum,mv$riemann2d_mid$estimate)
    compare("riemann_sum_2d_plot","corner sampled minimum",old_riemann2d$lower_sum,mv$riemann2d_min$estimate)
    compare("riemann_sum_2d_plot","corner sampled maximum",old_riemann2d$upper_sum,mv$riemann2d_max$estimate)
    old_prisms <- old$riemann_prisms3d(f=function(x,y)x+2*y,f1=function(x)0,f2=function(x)1,a=0,b=2,
      N=2,M=1,plot=FALSE,sample_n=1)
    compare("riemann_prisms3d","rectangular sampled minimum",old_prisms$sum_lower,mv$riemann2d_min$estimate)
    compare("riemann_prisms3d","rectangular sampled maximum",old_prisms$sum_upper,mv$riemann2d_max$estimate)
    compare("riemann_prisms3d","rectangular sampled mean",old_prisms$sum_mean,mv$riemann2d_mean$estimate)
    old_ladder_rate <- old$related_rates_grad(function(x,y)x^2+y^2-25,c(4,3),c(-1,NA),2)
    old_time_rate <- old$related_rates_grad(function(x,t)x^2-t,c(2,4),c(NA,1),1)
    compare("related_rates_grad","ladder missing rate",old_ladder_rate$rate,mv$rates_ladder$estimate)
    compare("related_rates_grad","explicit time coordinate",old_time_rate$rate,mv$rates_time$estimate)
    old_cylinder <- old$solid_cylindrical3d(function(t)0,function(t)1,function(r,t)0,function(r,t)2,
      0,2*pi,plot=FALSE,n_theta=4,n_u=4,n_v=4,compute_volume=TRUE,vol_method="adaptive")
    old_cup <- old$solid_cylindrical3d(function(t)0,function(t)1,function(r,t)0,function(r,t)1+r,
      0,pi/2,plot=FALSE,n_theta=4,n_u=4,n_v=4,compute_volume=TRUE,vol_method="adaptive")
    compare("solid_cylindrical3d","cylinder volume",old_cylinder$volume$estimate,mv$cylindrical_cylinder$estimate)
    compare("solid_cylindrical3d","radially varying height",old_cup$volume$estimate,mv$cylindrical_cup$estimate)
    old_tetra <- old$integrate_triple_general(function(x,y,z)1,0,1,function(x)0,function(x)1-x,
      function(x,y)0,function(x,y)1-x-y,n_outer=4,n_middle=4,n_inner=4,plot_xy_domain=FALSE)
    old_box <- old$integrate_triple_general(function(x,y,z)x+2*y+3*z,0,1,function(x)0,function(x)2,
      function(x,y)0,function(x,y)3,n_outer=4,n_middle=4,n_inner=4,plot_xy_domain=FALSE)
    compare("integrate_triple_general","tetrahedron volume",old_tetra$integral_value,mv$triple_tetra$estimate)
    compare("integrate_triple_general","asymmetric box signed field",old_box$integral_value,mv$triple_box$estimate)
    polar_warnings <- character()
    old_polar <- withCallingHandlers(old$integrate_double_polar(function(x,y)1,0,2*pi,
      function(t)0,function(t)1,n_theta=4,n_r=4,plot_domain=FALSE),
      warning=function(w) {
        polar_warnings <<- c(polar_warnings,conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    stopifnot(length(polar_warnings)>0,abs(old_polar$integral_value-49*pi/36)<1e-10)
    findings <- rbind(findings,data.frame(original="integrate_double_polar",
      issue="Simpson coefficients have n+2 entries: recycled vectors; unit disk n=4 gives 49*pi/36 instead of pi",
      observed_error=abs(old_polar$integral_value-pi)))
    old_double1 <- old$integrate_double_xy(function(x,y)x+2*y,0,1,function(x)0,
      function(x)x,n_outer=32,n_inner=32,plot_domain=FALSE)
    old_double2 <- old$integrate_double_xy(function(x,y)x+2*y,0,1,function(y)y,
      function(y)1,region_type="type2",n_outer=32,n_inner=32,plot_domain=FALSE)
    compare("integrate_double_xy","type I triangle signed integral",old_double1$integral_value,mv$double_type1$estimate)
    compare("integrate_double_xy","type II same triangle signed integral",old_double2$integral_value,mv$double_type2$estimate)
    old_lagrange <- old$lagrange_check(mv$lagrange_objective,mv$lagrange_constraint,c(.5,.5),h=.0025)
    old_lagrange_max <- old$lagrange_check(function(v)-sum(v^2),mv$lagrange_constraint,c(.5,.5),h=.0025)
    old_lagrange_invalid <- old$lagrange_check(mv$lagrange_objective,mv$lagrange_constraint,c(0,0),h=.0025)
    compare("lagrange_check","regular minimum multiplier",old_lagrange$lambda,unname(mv$lagrange_minimum$lambda))
    compare("lagrange_check","regular maximum multiplier",old_lagrange_max$lambda,unname(mv$lagrange_maximum$lambda))
    compare("lagrange_check","regular minimum stationarity residual",old_lagrange$ok_stationarity,
      max(abs(mv$lagrange_minimum$stationarity_residual)))
    for(ij in list(c(1,1),c(1,2),c(2,2)))compare("lagrange_check",
      paste("Lagrangian Hessian entry",paste(ij,collapse=",")),old_lagrange$H_L[ij[1],ij[2]],
      mv$lagrange_minimum$lagrangian_hessian[ij[1],ij[2]])
    stopifnot(old_lagrange_max$clasificacion=="indeterminado",
      old_lagrange_invalid$clasificacion=="minimo",old_lagrange_invalid$ok_feasible==1,
      env$vc_lagrange(mv$lagrange_objective,mv$lagrange_constraint,c(0,0),warn=FALSE)$classification=="infeasible")
    findings <- rbind(findings,data.frame(original="lagrange_check",
      issue="Original leading-minor rule reports indeterminate for the known regular maximum of -(x^2+y^2) on x+y=1; reduced tangent test recovers a conditional maximum.",observed_error=NA_real_),
      data.frame(original="lagrange_check",
      issue="Original reports a minimum at (0,0) for x^2+y^2 constrained by x+y=1, despite constraint residual 1; new assessment rejects infeasible candidates.",observed_error=1))
    lines <- readLines("development/CAPABILITIES.md",warn=FALSE)
    lines <- lines[grepl("^\\| [A-Za-z_0-9]+ \\| (I|II|III|IV) \\|",lines)]
    if(length(lines)) {
      cells <- lapply(strsplit(lines,"|",fixed=TRUE),trimws)
      map <- data.frame(original=vapply(cells,function(z)z[2],character(1)),
        phase=vapply(cells,function(z)z[3],character(1)),recovery=vapply(cells,function(z)z[4],character(1)))
      exports <- readLines(file.path(legacy_dir,"NAMESPACE"),warn=FALSE)
      exports <- sub("^export\\((.*)\\)$","\\1",exports[grepl("^export\\(",exports)])
      stopifnot(nrow(map)==47,setequal(map$original,exports))
    }
  }
  dir.create(export_dir,recursive=TRUE,showWarnings=FALSE)
  for(pair in list(c("PHASE-III-AUDIT.html","phase-III-audit.html"),
                   c("PHASE-III-AUDIT.csv","phase-III-audit.csv"),
                   c("phase-III-spherical-comparison.csv","phase-III-spherical-comparison.csv"))) {
    audit_file <- file.path("development",pair[1])
    if(file.exists(audit_file)) file.copy(audit_file,file.path(export_dir,pair[2]),overwrite=TRUE)
  }
  assets <- basename(tempfile("assets-"))
  for(name in names(plots)) htmlwidgets::saveWidget(plots[[name]],file.path(export_dir,paste0(name,".html")),
    selfcontained=FALSE,libdir=assets)
  if(legacy_available) {
    paired <- function(file,left,right,label,note="")writeLines(c(
      '<!doctype html><html lang="en"><meta charset="utf-8">',paste0('<title>',label,'</title><h1>',label,'</h1>'),
      '<p>Left: original 1.0.5; right: reconstruction. Playback controls are independent. Compare the same mathematical procedure, not identical APIs.</p>',
      if(nzchar(note))paste0('<p>',note,'</p>'),
      '<div style="display:flex;flex-wrap:wrap;gap:1%;">',
      paste0('<iframe title="Original" src="',left,'.html" style="width:49%;min-width:400px;height:850px;border:1px solid #aaa"></iframe>'),
      paste0('<iframe title="Reconstruction" src="',right,'.html" style="width:49%;min-width:400px;height:850px;border:1px solid #aaa"></iframe>'),
      '</div></html>'),file.path(export_dir,file))
    paired("compare-secants.html","legacy_secant","recovered_secant_forward","Secants: original and recovered")
    paired("compare-riemann.html","legacy_riemann_animation","phase2-riemann","Riemann rectangles: original and recovered")
    paired("compare-gradients.html","legacy_gradient_animation","phase3-gradient_animation","Gradient directions: original and recovered")
    paired("compare-planes.html","legacy_tangent_plane","phase3-linearization","Local model and tangent plane: original and recovered")
    paired("compare-coordinate-tangents.html","legacy_coordinate_tangents","phase3-linearization",
      "Coordinate tangent recovery and corrected surface mesh",
      "Audit finding: the original surface grid interchanges coordinates on this asymmetric example. Tangent slopes are compared numerically; the reconstruction corrects the surface ordering instead of copying it.")
    paired("compare-newton-systems.html","legacy_newton_system","phase3-newton_system",
      "Newton systems: recurrence recovery and corrected zero contours",
      "The original contour grid interchanges coordinates on this asymmetric system. Recurrence values match with the same analytic Jacobian and undamped steps; the reconstruction corrects contour orientation. Numerical-Jacobian safeguards and optional damping are separate enhancements.")
  }
  utils::write.csv(comparisons,file.path(export_dir,"legacy-comparison.csv"),row.names=FALSE)
  utils::write.csv(map,file.path(export_dir,"legacy-map.csv"),row.names=FALSE)
  utils::write.csv(lesson,file.path(export_dir,"connected-lesson.csv"),row.names=FALSE)
  utils::write.csv(findings,file.path(export_dir,"legacy-findings.csv"),row.names=FALSE)
  manuals <- c("../../../checks/multivariable-F/vectorialcalculus-manual.pdf",
               "../../../checks/multivariable-E/vectorialcalculus-manual.pdf",
               "../../../checks/multivariable-D/vectorialcalculus-manual.pdf",
               "../../../checks/multivariable-C/vectorialcalculus-manual.pdf",
               "../../../checks/multivariable-B/vectorialcalculus-manual.pdf",
               "../../../checks/multivariable-A/vectorialcalculus-manual.pdf",
               "../../../checks/continuity/vectorialcalculus-manual.pdf")
  manual <- manuals[file.exists(manuals)][1]
  manual_ready <- !is.na(manual)&&file.exists(manual)
  if(manual_ready) stopifnot(file.copy(manual,file.path(export_dir,"manual.pdf"),overwrite=TRUE))
  links <- function(names)paste0('<li><a href="',names,'.html">',names,'</a></li>')
  writeLines(c('<!doctype html><html lang="en"><meta charset="utf-8"><title>vectorialcalculus: continuity</title>',
    '<h1>vectorialcalculus: 1.0.5 to 2.0</h1>',
    '<p>One mathematical project, rebuilt progressively. Original functions are preserved; calculation and graphics now have separate interfaces.</p>',
    '<p>Automated analytical/geometry checks do not replace human visual review. No continuity or differentiability proof is implied.</p>',
    '<h2>One function connects the curriculum: f(x)=x^2</h2>',
    '<p>Move from secants to derivatives, differentials and Taylor; then from Riemann sums to integrals, accumulation and revolution. These are related views of the same function.</p>',
    '<table border="1" cellpadding="5"><tr><th>Procedure</th><th>Computed</th><th>Analytical</th></tr>',
    paste0('<tr><td>',lesson$capability,'</td><td>',format(lesson$value,digits=8),'</td><td>',format(lesson$analytical,digits=8),'</td></tr>'),'</table>',
    '<h2>Recovered procedures and original examples</h2><ul>',
    if(legacy_available)c('<li><a href="compare-secants.html">Compare original/recovered secants side by side</a></li>',
                         '<li><a href="compare-riemann.html">Compare original/recovered Riemann rectangles side by side</a></li>',
                         '<li><a href="compare-gradients.html">Compare original/recovered gradient directions side by side</a></li>',
                         '<li><a href="compare-planes.html">Compare original/recovered tangent planes side by side</a></li>',
                         '<li><a href="compare-coordinate-tangents.html">Coordinate tangents and the surface-grid audit finding</a></li>',
                         '<li><a href="compare-newton-systems.html">Newton system recurrence and corrected zero contours</a></li>'),
    links(names(plots)[!grepl("^phase[1234]-",names(plots))]),'</ul>',
    '<p>Central secants in the reconstruction pass through both actual endpoints. Legacy central lines used a different anchor.</p>',
    '<p>Legacy lower/upper endpoint labels are not general Darboux bounds. Comparison uses monotone x^2. Shifted-axis regions explicitly specify their second boundary.</p>',
    '<h2>Phase I: differential calculus</h2><ul>',links(names(plots)[grepl("^phase1-",names(plots))]),'</ul>',
    '<h2>Phase II: integral calculus</h2><ul>',links(names(plots)[grepl("^phase2-",names(plots))]),'</ul>',
    '<h2>Phase III: partials, gradients, directional restrictions and local models</h2>',
    '<p>The same one-variable derivative engine now follows coordinate slices and straight paths on f(x,y)=x^2+2*y^2. At (1,1), gradient=(2,4), and direction (3/5,4/5) gives 4.4.</p><ul>',
    links(names(plots)[grepl("^phase3-",names(plots))]),'</ul>',
    '<p>Partial derivatives do not prove total differentiability. Directional values are computed along their own restrictions; gradient dot products are conditional comparisons.</p>',
    '<p>The candidate total differential predicts change 1 for increment (0.1,0.2) at (1,1); the actual change is 1.09. Halving the increment reduces the quadratic remainder by a factor of four. A separate counterexample retains a nonvanishing normalized remainder.</p>',
    '<p>Derivative matrices retain initial, half and quarter step estimates. Jacobian rows identify outputs and columns identify inputs. Hessian mixed entries use direct rectangles; symmetry is imposed, not a smoothness proof. The large-offset example marks unresolved entries.</p>',
    '<p>Stationary searches show a quadratic minimum, a saddle and a degenerate inconclusive test. Extremum labels require caller-declared C2 regularity and passing numerical checks; they remain conditional. Colored history shows objective evaluations, not Newton iterations or a descent flow. No boundary extrema or complete root enumeration is claimed.</p>',
    '<p>Newton system views show accepted iterates on two zero level sets, optional step damping with rejected trials, and a singular-Jacobian failure. The numerical engine reuses vc_jacobian, solves linear systems without explicit inverses and never treats small steps alone as convergence. Original recurrence comparisons use the same analytic Jacobian without damping; an original contour-orientation defect is explicitly recorded and corrected.</p>',
    '<p>Lagrange views assess supplied candidates: minimum on a line, maximum on a circle and a singular constraint. Feasibility, stationarity and constraint rank precede conditional reduced-Hessian classification. The original multipliers and Lagrangian matrices are compared directly; its missed maximum and infeasible minimum labels are recorded as audit findings. Arrowheads are static annotations; legend toggles affect traces.</p>',
    '<p>The original private Hessian formula is compared directly without running its critical-point optimizer or random starts.</p>',
    '<h2>Phase IV: curves, vector fields, integrals and theorems</h2><ul>',
    links(names(plots)[grepl("^phase4-",names(plots))]),'</ul>',
    '<p>The helix examples cover traversal, Frenet directions, curvature 0.8, signed torsion 0.4 and osculating disks. The straight line has only a tangent. The ribbon leaves a gap at an inflection and retains the curve point.</p>',
    '<h2>Traceability</h2><ul><li><a href="connected-lesson.csv">Connected x^2 lesson</a></li><li><a href="legacy-comparison.csv">Numerical comparisons with original source</a></li>',
    '<li><a href="phase-III-audit.html">Phase III coverage audit and remaining work</a></li>',
    '<li><a href="legacy-map.csv">All 47 original exports and recovery destinations</a></li>',
    '<li><a href="legacy-findings.csv">Explicit legacy audit findings (not copied into reconstruction)</a></li></ul>',
    if(manual_ready)'<p><a href="manual.pdf">Reference manual (PDF)</a></p>' else '<p>Reference manual PDF has not been copied into this export.</p>',
    if(!legacy_available)'<p>Original source not found: legacy comparison was not executed. Supply legacy_dir pointing to the original repository.</p>' else
      '<p>Original source was loaded in an isolated environment. It was not modified, installed or attached over the reconstructed package.</p>',
    '<p>Keep this index, all HTML files and the shared asset directory together.</p></html>'),file.path(export_dir,"index.html"))
  invisible(list(plots=plots,comparisons=comparisons,map=map,legacy_executed=legacy_available,
    lesson=lesson,findings=findings,index=normalizePath(file.path(export_dir,"index.html")),forward=forward,central=central))
}

project_review <- run_project_review()
message("Review entry point: ",project_review$index)
# browseURL(project_review$index)
# project_review$comparisons
# print(project_review$plots$recovered_secant_forward)

