#' Directional derivative in any dimension, with optional 2D visualization
#'
#' @description
#' Computes a numerical directional derivative of a multivariate function
#' at a given point, along a specified direction. The derivative is
#' approximated using centered finite differences. If the dimension is two
#' and \code{plot = TRUE}, the function displays a local visualization of
#' the surface defined by \code{z = f(x, y)}, including the evaluation point,
#' the direction, and the curve traced along that direction.
#'
#' @details
#' The function accepts two types of input for the function \code{f}:
#' \itemize{
#'   \item a function of several numeric arguments, for example
#'         \code{f(x, y, z, ...)}, or
#'   \item a function that takes a single numeric vector, such as
#'         \code{f(x_vec)}.
#' }
#'
#' At the evaluation point \code{x0}, the function:
#' \itemize{
#'   \item normalizes the direction vector \code{v} to obtain a unit
#'         direction,
#'   \item computes forward and backward perturbations along this unit
#'         direction,
#'   \item evaluates the function at those perturbed points,
#'   \item estimates the directional derivative using a centered
#'         finite-difference formula.
#' }
#'
#' In two dimensions, if \code{plot = TRUE}, the function builds a small
#' rectangular window around \code{x0} and evaluates the function over a
#' fine grid to produce a strip of the surface. It then overlays:
#' \itemize{
#'   \item the base point,
#'   \item the selected direction,
#'   \item the trajectory of the directional curve,
#'   \item (optionally) a surface grid and other plot elements.
#' }
#'
#' @param f A function returning a numeric scalar. It may be defined either
#' as \code{f(x, y, ...)} with several numeric arguments, or as a function
#' of a numeric vector, \code{f(x_vec)}.
#' @param x0 Numeric vector giving the evaluation point. Its length
#' determines the dimension of the problem.
#' @param v Numeric vector giving the direction along which the directional
#' derivative is computed. Must be nonzero and have the same length as
#' \code{x0}.
#' @param h Numeric step size used for centered finite-difference
#' approximations.
#' @param plot Logical; if \code{TRUE} and \code{length(x0) == 2}, displays a
#' 3D visualization of the local surface using \pkg{plotly}.
#' @param x_window,y_window Numeric half-widths defining the size of the
#' rectangular window around the evaluation point in the 2D case.
#' @param n_s Integer giving the number of samples along the direction line
#' in the 2D visualization.
#' @param n_r Integer giving the number of samples across the strip in the
#' 2D visualization.
#' @param show_strip Logical; if \code{TRUE}, draws a surface strip around
#' the evaluation point.
#' @param strip_colorscale Character string specifying a \pkg{plotly}
#' colorscale for the strip.
#' @param strip_opacity Numeric value between 0 and 1 determining the
#' opacity of the strip.
#' @param show_surface_grid Logical; if \code{TRUE}, overlays a grid onto
#' the surface strip.
#' @param surface_grid_color Character string giving the grid line color.
#' @param surface_grid_width Numeric value giving the grid line width.
#' @param curve_line,point_marker,u_line,w_line Lists with \pkg{plotly}
#' style options for the directional curve, evaluation point, and auxiliary
#' lines.
#' @param t_range Numeric vector of length two giving the parameter range
#' for the directional curve in the 2D visualization.
#' @param scene,bg Lists specifying 3D scene and background options when
#' plotting.
#' @param tol Numeric tolerance used for detecting numerical degeneracies.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{D}}{The directional derivative at \code{x0} along the
#'         normalized direction.}
#'   \item{\code{v_hat}}{The normalized direction vector.}
#'   \item{\code{fx}, \code{fy}}{Centered partial derivatives in two
#'         dimensions (only returned when \code{length(x0) == 2}).}
#'   \item{\code{fig}}{A \pkg{plotly} visualization when \code{plot = TRUE},
#'         otherwise \code{NULL}.}
#' }
#'
#' @examples
#' # General n-dimensional usage:
#' f3 <- function(x, y, z) x^2 + y^2 + z
#' directional_derivative3d(f3, x0 = c(1, 0, 0), v = c(1, 1, 0))
#'
#' # Two-dimensional example without plotting (fast, no plotly required):
#' f2 <- function(x, y) x^2 + y^2
#' directional_derivative3d(f2, x0 = c(0, 0), v = c(1, 2), plot = FALSE)
#'
#' @export
directional_derivative3d <- function(
    f,
    x0,
    v,
    h = 1e-6,
    plot = FALSE,
    x_window = 2, y_window = 2,
    n_s = 180, n_r = 50,
    show_strip = TRUE,
    strip_colorscale = "Blues",
    strip_opacity = 0.30,
    show_surface_grid = TRUE,
    surface_grid_color = "rgba(60,80,200,0.25)",
    surface_grid_width = 1,
    curve_line   = list(color = "red",   width = 1),
    point_marker = list(color = "black", size = 3, symbol = "circle"),
    u_line = list(color = "black", width = 1.0),
    w_line = list(color = "black", width = 0.5),
    t_range = c(-2, 2),
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z")
    ),
    bg = list(paper = "white", plot = "white"),
    tol = 1e-12
) {
  # basic checks
  stopifnot(is.numeric(h), length(h) == 1L, is.finite(h), h > 0)
  if (!is.numeric(x0) || any(!is.finite(x0))) {
    stop("'x0' must be finite numeric.", call. = FALSE)
  }
  if (!is.numeric(v) || any(!is.finite(v)) || length(v) != length(x0)) {
    stop("'v' must be finite numeric and same length as 'x0'.", call. = FALSE)
  }
  nv <- sqrt(sum(v * v))
  if (nv < tol) {
    stop("'v' cannot be the zero vector.", call. = FALSE)
  }
  v_hat <- as.numeric(v / nv)

  # robust evaluator: accepts f(x, y, ...) or f(c(...))
  eval_f <- function(xx) {
    xx <- as.numeric(xx)
    out <- tryCatch(
      do.call(f, as.list(xx)),
      error = function(e) f(xx)
    )
    if (!is.numeric(out) || length(out) != 1L || !is.finite(out)) {
      stop("f(...) did not return a finite numeric scalar.", call. = FALSE)
    }
    out
  }

  # directional derivative
  g <- function(t) eval_f(x0 + t * v_hat)
  D <- (g(h) - g(-h)) / (2 * h)

  fx <- fy <- NA_real_
  fig <- NULL

  # 2D case: partial derivatives and optional plot
  if (length(x0) == 2L) {
    a <- x0[1]
    b <- x0[2]
    fx <- (eval_f(c(a + h, b)) - eval_f(c(a - h, b))) / (2 * h)
    fy <- (eval_f(c(a, b + h)) - eval_f(c(a, b - h))) / (2 * h)

    if (isTRUE(plot)) {
      if (!requireNamespace("plotly", quietly = TRUE)) {
        warning("To plot you need 'plotly' installed.", call. = FALSE)
      } else {
        # base line in x if direction is not vertical, else in y
        xs_line <- seq(a - x_window, a + x_window, length.out = n_s)

        build_strip <- function() {
          if (abs(v_hat[1]) >= 1e-14) {
            slope  <- v_hat[2] / v_hat[1]

            # line in the xy-plane centered at (a, b)
            y_center <- b + slope * (xs_line - a)

            r_seq <- seq(0, 1, length.out = n_r)

            # grid for the strip
            Xmat <- matrix(rep(xs_line, each = n_r), nrow = n_r)

            # interpolate between top edge (b + y_window) and the center line
            Ymat <- r_seq %o% rep(b + y_window, length(xs_line)) +
              (1 - r_seq) %o% y_center
          } else {
            y_line <- seq(b - y_window, b + y_window, length.out = n_s)
            r_seq  <- seq(0, 1, length.out = n_r)
            Ymat <- matrix(rep(y_line, each = n_r), nrow = n_r)
            Xmat <- matrix(
              rep(r_seq * (a + x_window) + (1 - r_seq) * a,
                  times = length(y_line)),
              nrow = n_r
            )
          }
          Zmat <- matrix(
            mapply(function(xx, yy) eval_f(c(xx, yy)),
                   as.numeric(Xmat), as.numeric(Ymat)),
            nrow = nrow(Xmat), ncol = ncol(Xmat)
          )
          list(X = Xmat, Y = Ymat, Z = Zmat)
        }

        plt <- plotly::plot_ly()
        if (isTRUE(show_strip)) {
          strip <- build_strip()
          contours_arg <- if (isTRUE(show_surface_grid)) {
            list(
              x = list(show = TRUE, color = surface_grid_color,
                       width = surface_grid_width),
              y = list(show = TRUE, color = surface_grid_color,
                       width = surface_grid_width),
              z = list(show = FALSE)
            )
          } else {
            NULL
          }
          plt <- plt |>
            plotly::add_surface(
              x = strip$X, y = strip$Y, z = strip$Z,
              colorscale = strip_colorscale, showscale = FALSE,
              opacity = strip_opacity, contours = contours_arg
            )
        }

        f0 <- eval_f(c(a, b))
        p0 <- c(a, b, f0)

        # vectors in R3
        u3 <- c(v_hat[1], v_hat[2], D)
        nu3 <- sqrt(sum(u3 * u3))
        if (nu3 < tol) {
          nu3 <- 1
        }
        u <- u3 / nu3
        w <- c(u[1], u[2], 0)

        add_seg <- function(plt_obj, p, v_vec, line_opt) {
          plotly::add_trace(
            plt_obj,
            x = c(p[1], p[1] + v_vec[1]),
            y = c(p[2], p[2] + v_vec[2]),
            z = c(p[3], p[3] + v_vec[3]),
            type = "scatter3d", mode = "lines",
            line = line_opt,
            showlegend = FALSE, hoverinfo = "none"
          )
        }

        plt <- plt |>
          plotly::add_trace(
            x = p0[1], y = p0[2], z = p0[3],
            type = "scatter3d", mode = "markers",
            marker = point_marker,
            showlegend = FALSE, hoverinfo = "none"
          )
        plt <- add_seg(plt, p0, u, u_line)
        plt <- add_seg(plt, p0, w, w_line)

        # directional curve (x(t), y(t), z(t) = f(x(t), y(t)))
        ts <- seq(t_range[1], t_range[2], length.out = 400)
        x_dir <- a + v_hat[1] * ts
        y_dir <- b + v_hat[2] * ts
        z_dir <- vapply(
          seq_along(ts),
          function(i) eval_f(c(x_dir[i], y_dir[i])),
          numeric(1)
        )

        plt <- plt |>
          plotly::add_trace(
            x = x_dir, y = y_dir, z = z_dir,
            type = "scatter3d", mode = "lines",
            line = curve_line,
            showlegend = FALSE, hoverinfo = "none"
          ) |>
          plotly::layout(
            title = sprintf("Directional derivative D = %.6g", D),
            scene = scene,
            paper_bgcolor = bg$paper,
            plot_bgcolor  = bg$plot
          )

        fig <- plt
        print(plt)
      }
    }
  }

  out <- list(D = D, v_hat = v_hat, fx = fx, fy = fy, fig = fig)
  class(out) <- c("directional_derivative", class(out))
  out
}

