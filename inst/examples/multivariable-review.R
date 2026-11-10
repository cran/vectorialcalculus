# Source from next.Rproj. Reuses scalar engines through explicit restrictions.
if(!file.exists("R/multivariable.R")) stop("Open next.Rproj first.")
multivariable_review <- local({
  env <- new.env(parent=baseenv())
  for(file in list.files("R",pattern="[.]R$",full.names=TRUE)) sys.source(file,env)
  f <- function(v)v[1]^2+2*v[2]^2
  point <- c(1,1)
  gradient <- env$vc_gradient(f,point)
  angles <- seq(0,2*pi,length.out=13)
  # Axis-aligned angles use exact zeros to avoid roundoff-sized active components.
  axes <- list(c(1,0),c(0,1),c(-1,0),c(0,-1),c(1,0))
  directions <- lapply(seq_along(angles),function(i) {
    u <- if((i-1)%%3==0)axes[[1+(i-1)%/%3]] else c(cos(angles[i]),sin(angles[i]))
    env$vc_directional(f,point,u,gradient=gradient)
  })
  directional <- env$vc_directional(f,point,c(3,4),gradient=gradient)
  partial <- env$vc_partial(f,point,2)
  scalar_slice <- function(t)f(c(1,t))
  scalar_derivative <- env$vc_derivative(scalar_slice,1,warn=FALSE)
  stopifnot(max(abs(gradient$estimate-c(2,4)))<1e-6,abs(directional$estimate-4.4)<1e-6,
    abs(partial$estimate-scalar_derivative$estimate)<1e-12,
    all(vapply(directions,function(d)abs(d$estimate-sum(c(2,4)*d$direction))<1e-6,logical(1))))
  f3 <- function(v)sum(v^2)
  gradient3 <- env$vc_gradient(f3,c(1,1,1))
  linearization <- env$vc_linearization(f,point,dx=c(.1,.2),gradient=gradient)
  counter_field <- function(v)if(all(v==0))0 else v[1]^3/sum(v^2)
  counter_model <- env$vc_linearization(counter_field,c(0,0),dx=c(.1,.1))
  vector_field <- function(v)c(field=unname(v[1]^2+2*v[2]^2),product=unname(v[1]*v[2]))
  jacobian <- env$vc_jacobian(vector_field,point)
  rectangular_field <- function(v)c(first=unname(v[1]^2+v[2]),second=unname(v[2]*v[3]))
  rectangular_jacobian <- env$vc_jacobian(rectangular_field,c(x=1,y=2,z=3))
  hessian_field <- function(v)v[1]^2+v[1]*v[2]+2*v[2]^2
  hessian <- env$vc_hessian(hessian_field,point,h=.01)
  hessian_offset <- env$vc_hessian(function(v)1e12+hessian_field(v),point,warn=FALSE)
  stationary_minimum <- env$vc_stationary_points(f,rbind(c(.1,.2),c(-.2,.1)),
    assume_c2=TRUE,warn=FALSE)
  saddle_field <- function(v)v[1]^2-v[2]^2
  stationary_saddle <- env$vc_stationary_points(saddle_field,c(.1,.2),assume_c2=TRUE,warn=FALSE)
  degenerate_field <- function(v)v[1]^2+v[2]^4
  stationary_inconclusive <- env$vc_stationary_points(degenerate_field,c(0,0),assume_c2=TRUE,warn=FALSE)
  newton_field <- function(v)c(v[1]^2+2*v[2]^2-1,v[1]-2*v[2])
  newton_jacobian <- function(v)rbind(c(2*v[1],4*v[2]),c(1,-2))
  newton_system <- env$vc_newton_system(newton_field,c(.8,.2))
  newton_reference <- env$vc_newton_system(newton_field,c(.8,.2),
    jacobian=newton_jacobian,backtrack=FALSE,f_tol=1e-10)
  damped_field <- function(v)c(v[1]^2-2,v[2])
  damped_system <- env$vc_newton_system(damped_field,c(.1,0),
    jacobian=function(v)diag(c(2*v[1],1)))
  singular_field <- function(v)c(v[1]^2+1,v[2])
  singular_system <- env$vc_newton_system(singular_field,c(0,1),
    jacobian=function(v)diag(c(2*v[1],1)))
  lagrange_objective <- function(v)sum(v^2)
  lagrange_constraint <- function(v)sum(v)-1
  lagrange_minimum <- env$vc_lagrange(lagrange_objective,lagrange_constraint,c(.5,.5),
    assume_c2=TRUE,hessian_h=.01,warn=FALSE)
  lagrange_maximum <- env$vc_lagrange(function(v)-sum(v^2),lagrange_constraint,c(.5,.5),
    assume_c2=TRUE,hessian_h=.01,warn=FALSE)
  circle_constraint <- function(v)sum(v^2)-1
  lagrange_circle <- env$vc_lagrange(function(v)v[1],circle_constraint,c(1,0),
    assume_c2=TRUE,hessian_h=.01,warn=FALSE)
  singular_constraint <- function(v)sum(v^2)
  lagrange_singular <- env$vc_lagrange(function(v)v[1],singular_constraint,c(0,0),
    assume_c2=TRUE,warn=FALSE)
  stopifnot(lagrange_minimum$classification=="minimum",abs(lagrange_minimum$lambda-1)<1e-7,
    lagrange_maximum$classification=="maximum",abs(lagrange_maximum$lambda+1)<1e-7,
    lagrange_circle$classification=="maximum",lagrange_singular$classification=="inconclusive")
  stopifnot(newton_system$converged,max(abs(newton_system$estimate-c(2,1)/sqrt(6)))<1e-8,
    damped_system$converged,damped_system$iterations$damping[2]==1/8,
    singular_system$status=="jacobian_singular",!singular_system$converged)
  stopifnot(nrow(stationary_minimum$points)==1,
    max(abs(stationary_minimum$points))<1e-5,
    stationary_minimum$table$classification=="minimum",
    stationary_saddle$table$classification=="saddle",
    stationary_inconclusive$table$classification=="inconclusive")
  stopifnot(max(abs(jacobian$estimate[1,]-gradient$estimate))<1e-8,
    max(abs(hessian$estimate-matrix(c(2,1,1,4),2)))<1e-7,!hessian_offset$all_reliable)
  stopifnot(abs(linearization$differential(c(.1,.2))-1)<1e-8,
    abs(linearization$model(c(1.1,1.2))-4)<1e-8,
    max(abs(linearization$comparison$remainder-c(.09,.0225,.005625)))<1e-8,
    max(abs(counter_model$comparison$normalized_remainder-1/(2*sqrt(2))))<1e-8)
  solid_sph_bridge <- env$vc_solid_spherical(function(t)0,function(t)pi,
    function(v).5+.1*cos(v[1])*sin(v[2]),function(v)1+.2*sin(v[1])*sin(v[2]),
    c(0,pi/2),n=4,n_u=3,n_v=2)
  solid_sph_shell <- env$vc_solid_spherical(function(t)0,function(t)pi,function(v)1,function(v)2,
    c(0,2*pi),n=16,n_u=12,n_v=3)
  solid_sph_ball <- env$vc_solid_spherical(function(t)0,function(t)pi,function(v)0,function(v)1,
    c(0,2*pi),n=16,n_u=12,n_v=2)
  solid_sph_sector <- env$vc_solid_spherical(function(t).1*t,function(t)1+.1*t,
    function(v).3,function(v)1+.1*sin(v[2]),c(0,pi/2),n=6,n_u=4,n_v=3)
  solid_cyl_sector <- env$vc_solid_cylindrical(function(t).5+.1*t,function(t)1+.2*t,
    function(v)v[1]-.2*v[2],function(v)2+v[1]+.3*v[2],c(0,pi/2),n=4,n_u=3,n_v=2)
  solid_cyl_annulus <- env$vc_solid_cylindrical(function(t)1,function(t)2,
    function(v)0,function(v)1,c(0,2*pi),n=20,n_u=3,n_v=3)
  solid_cyl_cup <- env$vc_solid_cylindrical(function(t)0,function(t)1+.2*cos(t),
    function(v)0,function(v)1+.5*v[2],c(0,2*pi),n=20,n_u=4,n_v=4)
  solid_box <- env$vc_solid_cartesian(function(x)0,function(x)1,function(v)0,function(v)2,
    c(0,3),n=2,n_u=2,n_v=3)
  solid_wedge <- env$vc_solid_cartesian(function(x)0,function(x)1-x,function(v)0,function(v)1,
    c(0,1),n=4,n_u=3,n_v=3)
  solid_curved <- env$vc_solid_cartesian(function(x)0,function(x)1+x,function(v)v[1]-v[2],
    function(v)2+v[1]+2*v[2],c(0,2),n=5,n_u=4,n_v=4)
  solid_box_volume <- env$vc_triple_integrate(function(v)1,c(0,3),function(x)0,function(x)1,
    function(v)0,function(v)2,n_outer=2,n_middle=2,n_inner=2)
  solid_wedge_volume <- env$vc_triple_integrate(function(v)1,c(0,1),function(x)0,function(x)1-x,
    function(v)0,function(v)1,n_outer=2,n_middle=2,n_inner=2)
  stopifnot(abs(solid_box_volume$estimate-6)<1e-10,abs(solid_wedge_volume$estimate-.5)<1e-10)
  region_bridge <- env$vc_region(function(x)x^2,function(x)2+x,c(0,2),n=4,n_u=3)
  region_triangle <- env$vc_region(function(x)0,identity,c(0,1),n=8,n_u=4)
  region_curved <- env$vc_region(function(x)x^2,sqrt,c(0,1),n=12,n_u=8)
  region_custom <- env$vc_region(function(x)-1,function(x)1+x,
    partition=env$vc_partition(breaks=c(0,.1,.4,1)),n_u=3)
  riemann2d_field <- function(v)v[1]+2*v[2]
  riemann2d_mid <- env$vc_riemann2d(riemann2d_field,c(0,2),c(0,1),nx=2,ny=1)
  riemann2d_min <- env$vc_riemann2d(riemann2d_field,c(0,2),c(0,1),nx=2,ny=1,method="sample_min")
  riemann2d_max <- env$vc_riemann2d(riemann2d_field,c(0,2),c(0,1),nx=2,ny=1,method="sample_max")
  riemann2d_mean <- env$vc_riemann2d(riemann2d_field,c(0,2),c(0,1),nx=2,ny=1,method="sample_mean")
  riemann2d_signed <- env$vc_riemann2d(function(v)v[1]-2*v[2],c(0,2),c(0,1),nx=4,ny=4)
  riemann2d_region <- env$vc_riemann2d(function(v)1,c(0,1),c(0,1),nx=6,ny=6,lower=function(x)0,upper=identity)
  riemann2d_refinement <- env$vc_riemann2d_convergence(function(v)sum(v^2),c(0,1),c(0,1),reference=2/3)
  stopifnot(abs(riemann2d_mid$estimate-4)<1e-10,abs(riemann2d_min$estimate-1)<1e-10,
    abs(riemann2d_max$estimate-7)<1e-10,all(diff(riemann2d_refinement$table$abs_error)<0))
  rates_constraint <- function(v)sum(v^2)-25
  rates_ladder <- env$vc_related_rates(rates_constraint,c(4,3),c(-1,NA),TRUE,warn=FALSE)
  rates_time_constraint <- function(v)v[1]^2-v[2]
  rates_time <- env$vc_related_rates(rates_time_constraint,c(2,4),c(NA,1),TRUE,warn=FALSE)
  rates_unresolved_constraint <- function(v)v[1]+1e-12*v[2]
  rates_unresolved <- env$vc_related_rates(rates_unresolved_constraint,c(0,0),c(1,NA),TRUE,warn=FALSE)
  stopifnot(rates_ladder$status=="solved",abs(rates_ladder$estimate-4/3)<1e-7,
    rates_time$status=="solved",abs(rates_time$estimate-.25)<1e-7,
    rates_unresolved$status=="denominator_unresolved")
  spherical_sphere <- env$vc_spherical_integrate(function(v)1,c(0,2*pi),function(t)0,
    function(t)pi,function(v)0,function(v)1,n_outer=2,n_middle=16,n_inner=2)
  spherical_hemisphere <- env$vc_spherical_integrate(function(v)v[3],c(0,2*pi),function(t)0,
    function(t)pi/2,function(v)0,function(v)1,n_outer=2,n_middle=16,n_inner=2)
  spherical_shell <- env$vc_spherical_integrate(function(v)1,c(0,2*pi),function(t)0,
    function(t)pi,function(v)1,function(v)2,n_outer=4,n_middle=16,n_inner=2)
  stopifnot(abs(spherical_sphere$estimate-4*pi/3)<3e-6,
    abs(spherical_hemisphere$estimate-pi/4)<2e-6,abs(spherical_shell$estimate-28*pi/3)<2e-5)
  cylindrical_cylinder <- env$vc_cylindrical_integrate(function(v)1,c(0,2*pi),function(t)0,
    function(t)1,function(v)0,function(v)2,n_outer=4,n_middle=2,n_inner=2)
  cylindrical_cup <- env$vc_cylindrical_integrate(function(v)1,c(0,pi/2),function(t)0,
    function(t)1,function(v)0,function(v)1+v[2],n_outer=4,n_middle=2,n_inner=2)
  cylindrical_annulus <- env$vc_cylindrical_integrate(function(v)1,c(0,2*pi),function(t)1,
    function(t)2,function(v)0,function(v)1,n_outer=8,n_middle=2,n_inner=2)
  triple_tetra <- env$vc_triple_integrate(function(v)1,c(0,1),function(t)0,function(t)1-t,
    function(v)0,function(v)1-sum(v),n_outer=2,n_middle=2,n_inner=2)
  triple_box <- env$vc_triple_integrate(function(v)v[1]+2*v[2]+3*v[3],c(0,1),function(t)0,
    function(t)2,function(v)0,function(v)3,n_outer=2,n_middle=2,n_inner=2)
  stopifnot(abs(triple_tetra$estimate-1/6)<1e-10,abs(triple_box$estimate-42)<1e-10)
  polar_disk <- env$vc_polar_integrate(function(v)1,c(0,2*pi),function(t)0,function(t)1)
  polar_sector <- env$vc_polar_integrate(function(v)v[1],c(0,pi/2),function(t)0,function(t)1)
  polar_annulus <- env$vc_polar_integrate(function(v)1,c(0,2*pi),function(t)1,function(t)2)
  stopifnot(abs(polar_disk$estimate-pi)<1e-10,abs(polar_annulus$estimate-3*pi)<1e-10)
  double_field <- function(v)v[1]+2*v[2]
  double_type1 <- env$vc_double_integrate(double_field,c(0,1),function(x)0,function(x)x)
  double_type2 <- env$vc_double_integrate(double_field,c(0,1),function(y)y,function(y)1,order="xy")
  double_curved <- env$vc_double_integrate(function(v)1,c(0,1),function(x)x^2,function(x)sqrt(x))
  stopifnot(abs(double_type1$estimate-2/3)<1e-10,abs(double_type2$estimate-2/3)<1e-10)
  plots <- list(stationary_surface_minimum=env$vc_plot_stationary_surface(stationary_minimum,f,xlim=c(-1,2),ylim=c(-2,3),n=5),
    stationary_surface_saddle=env$vc_plot_stationary_surface(stationary_saddle,saddle_field),
    stationary_surface_degenerate=env$vc_plot_stationary_surface(stationary_inconclusive,degenerate_field),solid_sph_shell=env$vc_plot_solid_spherical(solid_sph_shell,"both"),
    solid_sph_ball=env$vc_plot_solid_spherical(solid_sph_ball,faces="radial_outer"),
    solid_sph_sector=env$vc_plot_solid_spherical(solid_sph_sector,"both",faces=names(solid_sph_sector$faces)),
    solid_cyl_sector=env$vc_plot_solid_cylindrical(solid_cyl_sector,"both"),
    solid_cyl_annulus=env$vc_plot_solid_cylindrical(solid_cyl_annulus,
      faces=c("radial_inner","radial_outer","z_lower","z_upper")),
    solid_cyl_cup=env$vc_plot_solid_cylindrical(solid_cyl_cup,"wireframe",
      faces=c("radial_outer","z_lower","z_upper")),solid_box=env$vc_plot_solid_cartesian(solid_box,"both"),
    solid_wedge=env$vc_plot_solid_cartesian(solid_wedge),
    solid_curved_wire=env$vc_plot_solid_cartesian(solid_curved,"wireframe"),region_triangle=env$vc_plot_region(region_triangle),
    region_curved_plate=env$vc_plot_region(region_curved,"3d",plane_z=2),
    region_custom=env$vc_plot_region(region_custom),riemann2d_mid=env$vc_plot_riemann2d(riemann2d_mid),
    riemann2d_signed=env$vc_plot_riemann2d(riemann2d_signed),
    riemann2d_region=env$vc_plot_riemann2d(riemann2d_region),rates_ladder=env$vc_plot_related_rates(rates_ladder,rates_constraint,vector_scale=.2),
    rates_time=env$vc_plot_related_rates(rates_time,rates_time_constraint,vector_scale=.2),
    rates_unresolved=env$vc_plot_related_rates(rates_unresolved,rates_unresolved_constraint),spherical_sphere=env$vc_plot_spherical_integral(spherical_sphere),
    spherical_hemisphere=env$vc_plot_spherical_integral(spherical_hemisphere),
    spherical_shell=env$vc_plot_spherical_integral(spherical_shell),cylindrical_cylinder=env$vc_plot_cylindrical_integral(cylindrical_cylinder),
    cylindrical_cup=env$vc_plot_cylindrical_integral(cylindrical_cup),
    cylindrical_annulus_base=env$vc_plot_cylindrical_integral(cylindrical_annulus,"projection"),triple_tetra=env$vc_plot_triple_integral(triple_tetra),
    triple_projection=env$vc_plot_triple_integral(triple_tetra,"projection"),
    triple_box=env$vc_plot_triple_integral(triple_box),polar_disk=env$vc_plot_polar_integral(polar_disk),
    polar_sector=env$vc_plot_polar_integral(polar_sector),
    polar_annulus=env$vc_plot_polar_integral(polar_annulus),double_type1=env$vc_plot_double_integral(double_type1),
    double_type2=env$vc_plot_double_integral(double_type2),
    double_curved=env$vc_plot_double_integral(double_curved),
    gradient_animation=env$vc_plot_gradient(gradient,f,directions=directions,
      xlim=c(-1,3),ylim=c(-1,3),vector_scale=.3),
    directional_surface=env$vc_plot_directional(directional,f),
    gradient3=env$vc_plot_gradient(gradient3,f3),
    linearization=env$vc_plot_linearization(linearization,f,show_normal=TRUE,xlim=c(-1,3),ylim=c(-1,3)),
    linearization_counterexample=env$vc_plot_linearization(counter_model,counter_field,
      title="Counterexample: stable partials do not establish a total differential"),
    jacobian_matrix=env$vc_plot_matrix(jacobian),
    jacobian_rectangular=env$vc_plot_matrix(rectangular_jacobian),
    hessian_matrix=env$vc_plot_matrix(hessian),
    hessian_offset=env$vc_plot_matrix(hessian_offset,title="Large additive offset: inspect unresolved entries"),
    stationary_minimum=env$vc_plot_stationary(stationary_minimum,f,xlim=c(-1,1),ylim=c(-1,1)),
    stationary_saddle=env$vc_plot_stationary(stationary_saddle,saddle_field,xlim=c(-1,1),ylim=c(-1,1)),
    stationary_inconclusive=env$vc_plot_stationary(stationary_inconclusive,degenerate_field,
      xlim=c(-1,1),ylim=c(-1,1),show_history=FALSE),
    newton_system=env$vc_plot_newton_system(newton_system,newton_field,xlim=c(-1,1),ylim=c(-1,1)),
    newton_system_damped=env$vc_plot_newton_system(damped_system,damped_field,
      xlim=c(-.2,2),ylim=c(-1,1),show_rejected=TRUE),
    newton_system_singular=env$vc_plot_newton_system(singular_system,singular_field,
      xlim=c(-1,1),ylim=c(-1,1)),
    lagrange_minimum=env$vc_plot_lagrange(lagrange_minimum,lagrange_objective,lagrange_constraint,
      xlim=c(-.5,1.5),ylim=c(-.5,1.5),show_weighted=TRUE),
    lagrange_circle=env$vc_plot_lagrange(lagrange_circle,function(v)v[1],circle_constraint,
      xlim=c(-1.5,2),ylim=c(-1.3,1.3),show_weighted=TRUE),
    lagrange_singular=env$vc_plot_lagrange(lagrange_singular,function(v)v[1],singular_constraint,
      xlim=c(-1,1),ylim=c(-1,1)))
  built <- plotly::plotly_build(plots$gradient_animation)
  stopifnot(length(built$x$frames)==length(directions))
  stopifnot(length(plotly::plotly_build(plots$linearization)$x$frames)==3)
  list(solid_sph_bridge=solid_sph_bridge,solid_cyl_sector=solid_cyl_sector,solid_box=solid_box,solid_box_volume=solid_box_volume,solid_wedge_volume=solid_wedge_volume,region_bridge=region_bridge,riemann2d_mid=riemann2d_mid,riemann2d_min=riemann2d_min,riemann2d_max=riemann2d_max,
    riemann2d_mean=riemann2d_mean,rates_ladder=rates_ladder,rates_time=rates_time,cylindrical_cylinder=cylindrical_cylinder,cylindrical_cup=cylindrical_cup,triple_tetra=triple_tetra,triple_box=triple_box,double_type1=double_type1,double_type2=double_type2,
    field=f,point=point,gradient=gradient,directional=directional,partial=partial,
    scalar_derivative=scalar_derivative,directions=directions,linearization=linearization,
    counter_model=counter_model,jacobian=jacobian,rectangular_jacobian=rectangular_jacobian,
    hessian_field=hessian_field,hessian=hessian,hessian_offset=hessian_offset,
    stationary_minimum=stationary_minimum,stationary_saddle=stationary_saddle,
    stationary_inconclusive=stationary_inconclusive,newton_field=newton_field,
    newton_jacobian=newton_jacobian,newton_system=newton_system,
    newton_reference=newton_reference,lagrange_objective=lagrange_objective,
    lagrange_constraint=lagrange_constraint,lagrange_minimum=lagrange_minimum,
    lagrange_maximum=lagrange_maximum,plots=plots)
})
if(interactive()) print(multivariable_review$plots$gradient_animation)
# print(multivariable_review$plots$directional_surface)
# Use project-review.R to export everything in the common old/new menu.
