# Source this file from RStudio (R >= 4.1) with plotly installed.
# Use the R/ files of this next/ directory when available. Otherwise verify
# that the installed package is version >= 2.0.0 and exports the expected API;
# the legacy package cannot silently satisfy this script.
# After sourcing, open each Viewer plot explicitly, as needed:
# print(visual_review$plots$secant)
# print(visual_review$plots$tangent)
# print(visual_review$plots$differential)
# print(visual_review$plots$taylor)
# print(visual_review$plots$critical_points)
# print(visual_review$plots$concavity)
# print(visual_review$plots$newton)
# print(visual_review$plots$limit)
# In RStudio, sourcing also exports visual-review/index.html under getwd().
# browseURL(visual_review$files[["index"]]) opens the exported menu.
# run_visual_validation(package_dir = "/path/to/next", display = TRUE,
#                       export_dir = "visual-review")
# Review zoom, hover, legend toggles, labels and negative/zero increments.
# Assertions below verify geometry; they do not replace human visual review.

.vc_review_script <- local({
  paths <- vapply(sys.frames(), function(frame) {
    if (is.character(frame$ofile) && length(frame$ofile)) frame$ofile[1] else NA_character_
  }, character(1))
  paths <- paths[!is.na(paths)]
  paths <- paths[basename(paths) == "visual-validation.R"]
  if (length(paths)) normalizePath(paths[length(paths)], winslash="/", mustWork=TRUE) else NA_character_
})

run_visual_validation <- function(package_dir = NULL, display = interactive(),
                                  export_dir = NULL) {
  if (!requireNamespace("plotly", quietly = TRUE)) stop("Install plotly first.")
  if (is.null(package_dir)) {
    candidates <- c(".", "next", if (!is.na(.vc_review_script))
      dirname(dirname(dirname(.vc_review_script))))
    matches <- candidates[file.exists(file.path(candidates, "R", "vc_derivative.R"))]
    if (length(matches)) package_dir <- matches[1]
  }
  if (is.null(package_dir)) {
    needed <- c("vc_derivative","vc_secant","vc_tangent","vc_differential",
                "vc_plot_secant","vc_plot_tangent","vc_plot_differential",
                "vc_taylor","vc_plot_taylor","vc_critical_points","vc_plot_critical_points",
                "vc_concavity","vc_plot_concavity","vc_newton","vc_plot_newton",
                "vc_limit","vc_plot_limit")
    if (!requireNamespace("vectorialcalculus", quietly=TRUE) ||
        utils::packageVersion("vectorialcalculus") < "2.0.0" ||
        !all(needed %in% getNamespaceExports("vectorialcalculus")))
      stop("Supply package_dir pointing to next/ or install the reconstruction.")
    env <- asNamespace("vectorialcalculus")
  } else {
    if (!file.exists(file.path(package_dir,"R","vc_derivative.R")))
      stop("package_dir must point to the reconstruction source next/.")
    env <- new.env(parent = baseenv())
    for (file in list.files(file.path(package_dir, "R"), pattern = "[.]R$", full.names = TRUE))
      sys.source(file, envir = env)
  }
  quad <- function(x) x^2
  critical_function <- function(x) x^3-3*x
  concavity_function <- function(x) x^3
  root_function <- function(x) x^2-2
  limit_function <- function(x) sin(x)/x
  results <- list(secant = env$vc_secant(quad, 1, 2),
                  tangent = env$vc_tangent(sin, 0),
                  differential = env$vc_differential(quad, 1, 0.2),
                  taylor = env$vc_taylor(exp, 0, order=2, x=seq(-1,1,length.out=21)),
                  critical_points=env$vc_critical_points(critical_function,c(-2,2),n=81),
                  concavity=env$vc_concavity(concavity_function,c(-1,1),n=81),
                  newton=env$vc_newton(root_function,1),limit=env$vc_limit(limit_function,0))
  plots <- list(secant = env$vc_plot_secant(results$secant, quad, xlim = c(0,3)),
                tangent = env$vc_plot_tangent(results$tangent, sin, xlim = c(-2,2)),
                differential = env$vc_plot_differential(results$differential, quad,
                                                        xlim = c(0.5,1.5)),
                taylor = env$vc_plot_taylor(results$taylor, exp, xlim=c(-1,1)),
                critical_points=env$vc_plot_critical_points(results$critical_points,critical_function),
                concavity=env$vc_plot_concavity(results$concavity,concavity_function),
                newton=env$vc_plot_newton(results$newton,root_function),
                limit=env$vc_plot_limit(results$limit,limit_function))
  built <- lapply(plots, plotly::plotly_build)
  near <- function(a,b,tol = 1e-7) stopifnot(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol)))
  # Analytic references: secant y=3*x-2; tangent to sin at zero y=x;
  # for x^2 at 1 and dx=.2: dy=.4, delta_y=.44, remainder=.04.
  near(results$secant$slope, 3)
  near(built$secant$x$data[[2]]$y, 3*built$secant$x$data[[2]]$x - 2)
  near(built$secant$x$data[[3]]$x, c(1,2))
  near(built$secant$x$data[[3]]$y, c(1,4))
  near(results$tangent$slope, 1)
  near(built$tangent$x$data[[2]]$y, built$tangent$x$data[[2]]$x)
  near(built$tangent$x$data[[3]]$y, 0)
  near(results$differential$dy, .4)
  near(results$differential$delta_y, .44)
  near(results$differential$remainder, .04)
  near(built$differential$x$data[[4]]$x, c(1,1.2))
  near(built$differential$x$data[[5]]$y, c(1,1.4))
  near(built$differential$x$data[[6]]$y, c(1,1.44))
  for (name in c("secant","tangent","differential")) {
    curve <- built[[name]]$x$data[[1]]
    f <- if (name == "tangent") sin else quad
    near(curve$y, f(curve$x))
    stopifnot(identical(curve$connectgaps, FALSE))
  }
  # Taylor reference at zero: P2(x)=1+x+x^2/2; coefficients are numerical.
  tb <- built$taylor$x$data
  tx <- as.numeric(tb[[1]]$x)
  near(results$taylor$coefficients, c(1,1,.5), tol=1e-6)
  near(tb[[1]]$y, exp(tx))
  near(tb[[2]]$y, 1+tx+tx^2/2, tol=1e-6)
  near(tb[[2]]$y, results$taylor$polynomial(tx))
  near(tb[[3]]$x, 0)
  near(tb[[3]]$y, 1)
  near(tb[[4]]$y, exp(tx)-results$taylor$polynomial(tx))
  near(tb[[5]]$y, abs(exp(tx)-results$taylor$polynomial(tx)))
  stopifnot(length(tb)==5L, built$taylor$x$layout$xaxis$rangeslider$visible,
            all(vapply(tb,function(t) identical(t$xaxis,"x"),logical(1))))
  near(results$critical_points$points$x,c(-1,1),tol=1e-5)
  near(built$critical_points$x$data[[2]]$y,c(2,-2))
  near(built$concavity$x$data[[2]]$x,0)
  near(results$newton$estimate,sqrt(2))
  stopifnot(length(built$newton$x$frames)==nrow(results$newton$iterations)-1L)
  near(results$limit$estimate,1,tol=1e-6)
  near(built$limit$x$data[[4]]$y,results$limit$estimate)
  # Additional signed and degenerate increment geometry.
  for (dx in c(-.2, 0)) {
    a <- env$vc_differential(quad, 1, dx)
    b <- plotly::plotly_build(env$vc_plot_differential(a, quad))
    near(b$x$data[[5]]$y, c(1,1+2*dx))
    near(b$x$data[[6]]$y, c(1,(1+dx)^2))
  }
  message("PASS: analytic references, curves, lines, anchors and increments agree.")
  files <- character()
  if (!is.null(export_dir)) {
    dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
    export_dir <- normalizePath(export_dir,winslash="/",mustWork=TRUE)
    files <- stats::setNames(file.path(export_dir,paste0(names(plots),".html")),names(plots))
    for (name in names(plots)) {
      # Fresh asset directories make repeated exports safe without deleting
      # existing user files or hiding filesystem warnings.
      assets <- basename(tempfile(paste0(name,"_files-"), tmpdir=export_dir))
      htmlwidgets::saveWidget(plots[[name]], file.path(export_dir,paste0(name,".html")),
                              selfcontained=FALSE, libdir=assets)
    }
    index <- file.path(export_dir,"index.html")
    writeLines(c('<!doctype html><html lang="en"><meta charset="utf-8">',
      '<title>vectorialcalculus graphical review</title><body>',
      '<h1>Interactive graphical review</h1>',
      '<p>Geometry checked automatically; human visual review is still required.</p><ul>',
      '<li><a href="secant.html">Secant: x^2, slope 3 through (1,1) and (2,4)</a></li>',
      '<li><a href="tangent.html">Tangent: sin(x) at zero, approximate slope 1</a></li>',
      '<li><a href="differential.html">Differential: x^2 at 1, dx=.2, dy=.4, delta_y=.44</a></li>',
      '<li><a href="taylor.html">Taylor: exp(x) at zero, P2 approximately 1+x+x^2/2</a></li>',
      '<li><a href="critical_points.html">Critical points: x^3-3x, maximum at -1 and minimum at 1</a></li>',
      '<li><a href="concavity.html">Concavity: x^3, smooth inflection candidate at 0</a></li>',
      '<li><a href="newton.html">Newton: x^2-2, root candidate sqrt(2), computed tangent steps</a></li>',
      '<li><a href="limit.html">Limit: sin(x)/x approaching 0, compatible value 1</a></li>',
      '</ul><p>Keep each HTML file together with its accompanying asset directories.</p></body></html>'),
      index, useBytes=TRUE)
    files <- c(index=index,files)
    message("HTML review menu: ",index)
  }
  if (isTRUE(display)) for (p in plots) print(p)
  invisible(list(results = results, plots = plots, files = files, verified = TRUE))
}

visual_review <- run_visual_validation(export_dir =
  if (interactive()) file.path(getwd(),"visual-review") else NULL)
