# Open next.Rproj, install plotly if needed, then Source this file.
# Geometry/analytical checks are not a human visual inspection.
if (!file.exists("R/riemann.R")) stop("Open next.Rproj before sourcing this script.")
integral_review <- local({
  env <- new.env(parent=baseenv())
  for(file in list.files("R",pattern="[.]R$",full.names=TRUE)) sys.source(file,envir=env)
  f <- function(x)x^2
  a <- env$vc_riemann_convergence(f,c(0,2),reference=8/3)
  stopifnot(all(diff(a$table$abs_error)<0))
  p <- env$vc_plot_riemann(a,f)
  built <- plotly::plotly_build(p)
  stopifnot(length(built$x$frames)==length(a$results))
  comparison <- env$vc_integral_compare(f,c(0,2),n=4,reference=8/3)
  stopifnot(abs(comparison$results$simpson$estimate-8/3)<1e-12,
            abs(comparison$results$adaptive$estimate-8/3)<1e-12)
  plots <- lapply(comparison$results,function(z)env$vc_plot_integral(z,f))
  for(m in c("trapezoid","simpson")) {
    b <- plotly::plotly_build(plots[[m]])
    d <- comparison$results[[m]]$table
    stopifnot(isTRUE(all.equal(as.numeric(b$x$data[[2]]$x),d$x)),
              isTRUE(all.equal(as.numeric(b$x$data[[2]]$y),d$value)))
    if(m=="simpson") stopifnot(max(abs(b$x$data[[1]]$y-b$x$data[[1]]$x^2))<1e-12)
  }
  limits <- seq(-1,2,length.out=13)
  accumulation <- env$vc_accumulation(f,x=limits)
  fundamental <- env$vc_fundamental(f,x=limits)
  stopifnot(max(abs(accumulation$table$estimate-limits^3/3))<1e-12,
            max(abs(fundamental$table$estimate-limits^2))<1e-7)
  ap <- env$vc_plot_accumulation(accumulation,f,fundamental=fundamental)
  bp <- plotly::plotly_build(ap)
  stopifnot(length(bp$x$frames)==length(limits))
  signed <- env$vc_accumulation(identity,x=seq(-1,1,length.out=9))
  stopifnot(max(abs(signed$table$estimate-signed$table$x^2/2))<1e-12)
  sp <- env$vc_plot_accumulation(signed,identity)
  area <- env$vc_area(identity,c(-1,1),breaks=0)
  cone <- env$vc_revolution(identity,c(0,1))
  sphere_f <- function(t)sqrt(1-t^2)
  sphere <- env$vc_revolution(sphere_f,c(-1,1))
  shell_f <- function(t)1-t
  shells <- env$vc_revolution(shell_f,c(0,1),axis="vertical")
  outer_f <- function(t)2
  inner_f <- function(t)1
  annulus <- env$vc_revolution(outer_f,c(0,1),g=inner_f)
  arc_result <- env$vc_arc_length(identity,c(0,1))
  surface <- env$vc_surface_area(identity,c(0,1),df=function(t)1,end_caps=TRUE)
  stopifnot(abs(area$estimate-1)<1e-12,abs(area$signed_estimate)<1e-12,
    abs(cone$estimate-pi/3)<1e-12,abs(sphere$estimate-4*pi/3)<1e-10,
    abs(shells$estimate-pi/3)<1e-12,abs(annulus$estimate-3*pi)<1e-12,
    abs(arc_result$estimate-sqrt(2))<1e-8,abs(surface$estimate-pi*(sqrt(2)+1))<1e-12)
  applications <- list(area=area,cone=cone,sphere=sphere,shells=shells,
                       annulus=annulus,arc_length=arc_result,surface_area=surface)
  app_plots <- list(area=env$vc_plot_area(area,identity),
    cone=env$vc_plot_revolution(cone,identity,sections=seq(.1,1,length.out=7)),
    sphere=env$vc_plot_revolution(sphere,sphere_f,sections=seq(-.9,.9,length.out=7)),
    shells=env$vc_plot_revolution(shells,shell_f,sections=seq(.1,1,length.out=7)),
    annulus=env$vc_plot_revolution(annulus,outer_f,g=inner_f,sections=c(.25,.5,.75)))
  cb <- plotly::plotly_build(app_plots$cone)
  ob <- Filter(function(z)identical(z$name,"Outer radial boundary"),cb$x$data)[[1]]
  stopifnot(length(cb$x$frames)==7,max(abs(ob$y^2+ob$z^2-ob$x^2))<1e-12)
  list(result=a,plot=p,comparison=comparison,accumulation=accumulation,
       fundamental=fundamental,signed=signed,applications=applications,
       plots=c(list(riemann=p),plots,list(accumulation=ap,signed_accumulation=sp),app_plots))
})
if(interactive()) print(integral_review$plot)
# Examine additional widgets one at a time in RStudio:
# print(integral_review$plots$trapezoid)
# print(integral_review$plots$simpson)
# print(integral_review$plots$adaptive)
# print(integral_review$plots$accumulation)
# print(integral_review$plots$signed_accumulation)
# integral_review$comparison$table
# integral_review$fundamental$table
# print(integral_review$plots$area)
# print(integral_review$plots$cone)
# print(integral_review$plots$sphere)
# print(integral_review$plots$shells)
# print(integral_review$plots$annulus)
# integral_review$applications$surface_area

# Optional export of all eleven widgets. Keep the generated asset directories.
# export_integral_review()
# browseURL(normalizePath("integral-review/index.html"))
export_integral_review <- function(directory="integral-review") {
  dir.create(directory,recursive=TRUE,showWarnings=FALSE)
  for(name in names(integral_review$plots)) {
    htmlwidgets::saveWidget(integral_review$plots[[name]],
      file.path(directory,paste0(name,".html")),selfcontained=FALSE,
      libdir=basename(tempfile("assets-")))
  }
  links <- paste0('<li><a href="',names(integral_review$plots),'.html">',
                  names(integral_review$plots),'</a></li>')
  writeLines(c('<!doctype html><html lang="en"><meta charset="utf-8">',
    '<title>Integral review</title><h1>Numerical integration</h1>',
    '<p>Quadrature: x^2 on [0,2], exact integral 8/3. Accumulation from zero: x^3/3 for x^2; x^2/2 for identity. Limits below zero reverse orientation.</p>',
    '<p>Automated geometry checks do not replace human visual review.</p>',
    '<p>Applications: identity on [-1,1] has geometric area 1 and signed integral 0; unit cone V=pi/3; unit sphere V=4*pi/3; shells give the same cone volume; annular cylinder V=3*pi.</p>',
    '<ul>',links,'</ul></html>'),file.path(directory,"index.html"))
  invisible(normalizePath(file.path(directory,"index.html")))
}
