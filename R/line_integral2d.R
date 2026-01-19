#' Line integral of a scalar field along a planar curve, with optional 3D visualization
#'
#' @description
#' Computes a numerical line integral of a scalar field along a parametric
#' curve in the plane. The function integrates the quantity
#' \code{f(r(t)) * speed(t)}, where \code{r(t)} is the parametric curve
#' and \code{speed(t)} is the length of its derivative. Optionally, it
#' produces a 3D visualization that includes a surface representing
#' \code{z = f(x, y)}, the curve itself in the plane, a lifted version of
#' the curve showing \code{z = f(x(t), y(t))}, and a vertical curtain
#' over the curve.
#'
#' @details
#' The function evaluates the scalar field along the curve and computes an
#' approximation of the derivative of \code{r(t)} using centered finite
#' differences. The integral can be computed either through a built-in
#' adaptive integration routine or via a composite Simpson rule with a user
#' specified number of subintervals.
#'
#' For visualization, the function:
#' \itemize{
#'   \item builds a rectangular surface \code{z = f(x, y)} using only the
#'         endpoints of the curve,
#'   \item plots the curve in the plane,
#'   \item plots a lifted copy of the curve where the height is given by
#'         the scalar field,
#'   \item optionally constructs a vertical curtain over the curve by
#'         extruding the height values.
#' }
#'
#' @param f A scalar field, given as \code{function(x, y)} returning a
#' numeric value.
#' @param r A parametric curve, given as \code{function(t)} that returns a
#' numeric vector of length two, interpreted as \code{c(x, y)}.
#' @param a,b Numeric scalars giving the parameter limits, with \code{b > a}.
#' @param plot Logical; if \code{TRUE}, produces a visualization using
#' \pkg{plotly}.
#' @param n_curve Number of sample points along the curve for plotting.
#' @param n_curtain_v Number of subdivisions in the vertical direction for
#' the curtain.
#' @param n_surf_x,n_surf_y Grid resolution for sampling the surface
#' \code{z = f(x, y)}.
#' @param colorscale Color scale for the surface. It may be a \pkg{plotly}
#' scale name, a single color, or a vector of colors defining a gradient.
#' @param surface_opacity Numeric opacity for the surface, between 0 and 1.
#' @param show_surface_grid Logical; if \code{TRUE}, overlays grid lines on
#' the surface.
#' @param surface_grid_color Character string giving the color of surface
#' grid lines.
#' @param surface_grid_width Numeric width of surface grid lines.
#' @param curtain Numeric value between 0 and 1 giving the opacity of the
#' vertical curtain. A value of zero disables the curtain.
#' @param curtain_colorscale Color scale for the curtain, using the same
#' formats accepted by \code{colorscale}.
#' @param curve_color Color for the curve drawn at height zero.
#' @param curve_width Width of the curve drawn at height zero.
#' @param lifted_color Color for the lifted curve whose height is
#' \code{f(x(t), y(t))}.
#' @param lifted_width Width of the lifted curve.
#' @param h Step size used for centered finite differences when computing
#' the derivative of \code{r(t)}. If \code{NULL}, a default is chosen
#' automatically.
#' @param method Integration method; may be \code{"adaptive"} (using
#' \code{stats::integrate}) or \code{"simpson"} (using a composite Simpson
#' rule).
#' @param n_simpson Number of subintervals for the Simpson method. This
#' value is automatically adjusted to be even.
#' @param scene List configuring the 3D scene in \pkg{plotly}.
#' @param bg List with background color settings for the figure, with entries
#' such as \code{paper} and \code{plot}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{value}: the numeric value of the line integral,
#'   \item \code{fig}: a \pkg{plotly} figure when \code{plot = TRUE}, or
#'         \code{NULL} otherwise.
#' }
#'
#' @examples
#' f <- function(x, y) x^2 + y^2
#' r <- function(t) c(t*cos(t), t*sin(3*t))
#' line_integral2d(f, r, a = 0, b = 2*pi, plot = FALSE)
#'
#' @export
line_integral2d <- function(
    f, r, a, b,
    plot = TRUE,
    n_curve = 400,
    n_curtain_v = 40,
    n_surf_x = 80, n_surf_y = 80,
    colorscale = "Viridis",
    surface_opacity = 0.65,
    show_surface_grid  = TRUE,
    surface_grid_color = "rgba(80,80,80,0.25)",
    surface_grid_width = 1,
    curtain = 0.4,
    curtain_colorscale = c("white", "#2a9d8f"),
    curve_color = "black",
    curve_width = 3,
    lifted_color = "red",
    lifted_width = 2,
    h = NULL,
    method = c("adaptive", "simpson"),
    n_simpson = 1000,
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z")
    ),
    bg = list(paper = "white", plot = "white")
) {
  method <- match.arg(method)

  # --- basic checks -----------------------------------------------------
  if (!is.function(f)) stop("'f' must be function(x, y) -> scalar.", call. = FALSE)
  if (!is.function(r)) stop("'r' must be function(t) -> c(x, y).", call. = FALSE)
  if (!is.numeric(a) || !is.numeric(b) ||
      length(a) != 1L || length(b) != 1L ||
      !is.finite(a) || !is.finite(b) || b <= a) {
    stop("'a' and 'b' must be finite scalars with b > a.", call. = FALSE)
  }
  if (!is.numeric(n_curve) || n_curve < 10) {
    stop("'n_curve' must be a numeric value >= 10.", call. = FALSE)
  }
  n_curve <- as.integer(n_curve)

  if (!is.null(h)) {
    if (!is.numeric(h) || length(h) != 1L || !is.finite(h) || h <= 0) {
      stop("'h' must be a positive numeric scalar or NULL.", call. = FALSE)
    }
  }

  # --- helpers to validate outputs -------------------------------------
  eval_r <- function(t) {
    out <- r(t)
    if (!is.numeric(out) || length(out) != 2L || any(!is.finite(out))) {
      stop("r(t) must return a finite numeric vector of length 2.", call. = FALSE)
    }
    as.numeric(out)
  }

  eval_f <- function(x, y) {
    out <- f(x, y)
    if (!is.numeric(out) || length(out) != 1L || !is.finite(out)) {
      stop("f(x, y) must return a finite numeric scalar.", call. = FALSE)
    }
    as.numeric(out)
  }

  # --- color helpers ----------------------------------------------------
  is_color <- function(x) {
    if (!is.character(x) || length(x) != 1L) return(FALSE)
    if (grepl("^rgba?\\(", x, ignore.case = TRUE)) return(TRUE)
    if (grepl("^#([0-9A-Fa-f]{6}|[0-9A-Fa-f]{8})$", x)) return(TRUE)
    ok <- TRUE
    tryCatch(grDevices::col2rgb(x), error = function(...) ok <<- FALSE)
    ok
  }

  as_rgba <- function(col, alpha = NULL) {
    if (grepl("^rgba?\\(", col, ignore.case = TRUE)) return(col)
    rgb <- grDevices::col2rgb(col) / 255
    a_  <- if (is.null(alpha)) 1 else max(0, min(1, alpha))
    sprintf("rgba(%g,%g,%g,%g)", 255 * rgb[1], 255 * rgb[2], 255 * rgb[3], a_)
  }

  as_colorscale <- function(x, alpha = NULL) {
    if (is.null(x)) return(NULL)
    if (is.list(x) && length(x) >= 2L && is.numeric(x[[1]][[1]])) return(x)
    if (is.character(x) && length(x) == 1L) {
      if (is_color(x)) {
        ccol <- as_rgba(x, alpha)
        return(list(list(0, ccol), list(1, ccol)))
      }
      return(x)
    }
    if (is.character(x) && length(x) > 1L) {
      cols <- vapply(x, as_rgba, character(1), alpha = alpha)
      pos  <- seq(0, 1, length.out = length(cols))
      return(lapply(seq_along(cols), function(i) list(pos[i], cols[i])))
    }
    stop("Unrecognized 'colorscale' format.", call. = FALSE)
  }

  # --- numerical derivative of r(t) -------------------------------------
  if (is.null(h)) h <- (b - a) * 1e-4 + 1e-6

  rprime <- function(t) {
    tt1 <- pmax(a, pmin(b, t - h))
    tt2 <- pmax(a, pmin(b, t + h))
    (eval_r(tt2) - eval_r(tt1)) / (tt2 - tt1)
  }

  speed <- function(t) {
    v <- rprime(t)
    sqrt(sum(v * v))
  }

  # --- integrand: f(r(t)) * ||r'(t)|| -----------------------------------
  integrand <- function(t) {
    p  <- eval_r(t)
    fx <- eval_f(p[1], p[2])
    fx * speed(t)
  }

  integrand_vec <- Vectorize(integrand, "t")

  # --- integral value ---------------------------------------------------
  value <- if (method == "adaptive") {
    stats::integrate(integrand_vec, lower = a, upper = b, rel.tol = 1e-6)$value
  } else {
    n <- as.integer(n_simpson)
    if (n %% 2L == 1L) n <- n + 1L
    tt <- seq(a, b, length.out = n + 1L)
    yy <- integrand_vec(tt)
    hS <- (b - a) / n
    hS * (yy[1] + yy[length(yy)] +
            4 * sum(yy[seq(2, n, by = 2L)]) +
            2 * sum(yy[seq(3, n - 1L, by = 2L)])) / 3
  }

  # --- plotting data ----------------------------------------------------
  fig <- NULL
  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("To plot you need 'plotly' installed.", call. = FALSE)
    } else {
      t_curve <- seq(a, b, length.out = n_curve)
      Rxy     <- t(vapply(t_curve, eval_r, numeric(2)))
      cx <- Rxy[, 1]
      cy <- Rxy[, 2]
      cz0 <- rep(0, length(t_curve))

      cf <- vapply(
        seq_along(t_curve),
        function(i) eval_f(cx[i], cy[i]),
        numeric(1)
      )
      cz <- cf

      # surface z = f(x, y) over rectangle joining r(a) and r(b)
      ra <- eval_r(a)
      rb <- eval_r(b)
      xs <- seq(ra[1], rb[1], length.out = n_surf_x)
      ys <- seq(ra[2], rb[2], length.out = n_surf_y)
      X  <- matrix(rep(xs, each = n_surf_y), nrow = n_surf_y)
      Y  <- matrix(rep(ys, times = n_surf_x), nrow = n_surf_y)
      Z  <- matrix(NA_real_, nrow = n_surf_y, ncol = n_surf_x)
      for (j in seq_len(n_surf_x)) {
        Z[, j] <- vapply(ys, function(y) eval_f(xs[j], y), numeric(1))
      }

      plt <- plotly::plot_ly()

      # curtain over the curve
      if (is.numeric(curtain) && curtain > 0) {
        vseq <- seq(0, 1, length.out = n_curtain_v)
        Xc   <- matrix(rep(cx, each = n_curtain_v), nrow = n_curtain_v)
        Yc   <- matrix(rep(cy, each = n_curtain_v), nrow = n_curtain_v)
        Zc   <- matrix(rep(cz, each = n_curtain_v), nrow = n_curtain_v) *
          matrix(rep(vseq, times = n_curve), nrow = n_curtain_v)

        plt <- plt |>
          plotly::add_surface(
            x = Xc, y = Yc, z = Zc,
            colorscale = as_colorscale(curtain_colorscale),
            showscale  = FALSE,
            opacity    = curtain
          )
      }

      # optional grid on surface z = f(x, y)
      contours_arg <- if (isTRUE(show_surface_grid)) list(
        x = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
        y = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
        z = list(show = FALSE)
      ) else NULL

      # base surface
      plt <- plt |>
        plotly::add_surface(
          x = X, y = Y, z = Z,
          colorscale = as_colorscale(colorscale),
          showscale  = FALSE,
          opacity    = surface_opacity,
          contours   = contours_arg
        )

      # curve at z = 0 and lifted curve
      plt <- plt |>
        plotly::add_trace(
          x = cx, y = cy, z = cz0,
          type = "scatter3d", mode = "lines",
          line = list(color = curve_color, width = curve_width),
          showlegend = FALSE
        ) |>
        plotly::add_trace(
          x = cx, y = cy, z = cz,
          type = "scatter3d", mode = "lines",
          line = list(color = lifted_color, width = lifted_width),
          showlegend = FALSE
        ) |>
        plotly::layout(
          title = sprintf("Line integral (value ~ %.6g)", value),
          scene = scene,
          paper_bgcolor = bg$paper,
          plot_bgcolor  = bg$plot
        )

      fig <- plt
      print(fig)
    }
  }

  list(value = value, fig = fig)
}
