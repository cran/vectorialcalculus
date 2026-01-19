# Internal utility: centered finite difference for f'(t)
# Not exported.
.d_central <- function(f, t, h) {
  (f(t + h) - f(t - h)) / (2 * h)
}

#' Numeric arc length of a 3D parametric curve
#'
#' Computes a numerical approximation to the arc length of the parametric
#' curve \eqn{(X(t), Y(t), Z(t))} on the interval \eqn{[a, b]} by integrating
#' the speed \eqn{\sqrt{(dx/dt)^2 + (dy/dt)^2 + (dz/dt)^2}}.
#'
#' Derivatives are approximated by centered finite differences and the
#' integral is computed either by Romberg integration (via \pkg{pracma})
#' or by \code{\link[stats]{integrate}} from base R. Optionally, the
#' curve can be visualized with [plot_curve3d()].
#'
#' @param X,Y,Z Functions of one variable \code{t} defining the parametric
#'   curve coordinates.
#' @param a,b Numeric parameter limits for \code{t}.
#' @param h Numeric step size for centered finite differences used to approximate
#'   the derivatives \eqn{dX/dt}, \eqn{dY/dt}, and \eqn{dZ/dt}.
#' @param method_int Character string. Either \code{"romberg"} (requires
#'   the \pkg{pracma} package) or \code{"integrate"} (base R).
#' @param n_samples Integer. Number of sample points used when plotting
#'   the curve (if \code{plot = TRUE}).
#' @param plot Logical. If \code{TRUE}, the function also produces a 3D
#'   visualization of the curve using [plot_curve3d()].
#' @param plot_mode Character string passed to [plot_curve3d()] as the
#'   \code{mode} argument.
#' @param plot_line List with line styling options passed to [plot_curve3d()].
#' @param plot_marker Optional list with marker styling options passed
#'   to [plot_curve3d()], or \code{NULL}.
#' @param plot_title Optional title for the plot. If \code{NULL}, a title
#'   including the estimated arc length is generated.
#' @param plot_scene List specifying 3D axes and options, passed to
#'   [plot_curve3d()].
#' @param plot_bg List with background colors, passed to [plot_curve3d()].
#'
#' @return
#' A single numeric value: the approximated arc length of the curve
#' on the interval \eqn{[a, b]}.
#'
#' @seealso [curve_sample3d()], [plot_curve3d()]
#'
#' @examples
#' X <- function(t) t^2 * cos(t)
#' Y <- function(t) t^3 * sin(3 * t)
#' Z <- function(t) t
#' arc_length3d(X, Y, Z, 0, 2 * pi)
#'
#' # \donttest{
#' # if (requireNamespace("plotly", quietly = TRUE)) {
#' #   arc_length3d(
#' #     X, Y, Z, 0, 2 * pi,
#' #     plot = TRUE,
#' #     plot_line = list(color = "red", width = 3),
#' #     n_samples = 300
#' #   )
#' # }
#' # }
#'
#' @importFrom stats integrate
#' @export
arc_length3d <- function(
    X, Y, Z, a, b,
    h = 1e-6,
    method_int = c("romberg", "integrate"),
    n_samples = 400,
    plot = FALSE,
    plot_mode  = "lines",
    plot_line  = list(color = "blue", width = 3, dash = "solid"),
    plot_marker = NULL,
    plot_title  = NULL,
    plot_scene  = list(
      xaxis = list(title = "x(t)"),
      yaxis = list(title = "y(t)"),
      zaxis = list(title = "z(t)")
    ),
    plot_bg = list(paper = "white", plot = "white")
) {
  method_int <- match.arg(method_int)

  speed <- function(t) {
    dx <- .d_central(X, t, h)
    dy <- .d_central(Y, t, h)
    dz <- .d_central(Z, t, h)
    sqrt(dx * dx + dy * dy + dz * dz)
  }

  length_val <- switch(
    method_int,
    romberg = {
      if (!requireNamespace("pracma", quietly = TRUE)) {
        stop("For method_int = 'romberg' you need to install the 'pracma' package.")
      }
      pracma::romberg(speed, a, b)$value
    },
    integrate = stats::integrate(speed, lower = a, upper = b, rel.tol = 1e-8)$value
  )

  if (isTRUE(plot)) {
    data <- curve_sample3d(X, Y, Z, a, b, n_samples = n_samples)
    plot_curve3d(
      data,
      mode   = plot_mode,
      line   = plot_line,
      marker = plot_marker,
      title  = if (is.null(plot_title)) {
        paste0("Arc length =", signif(length_val, 6))
      } else {
        plot_title
      },
      scene = plot_scene,
      bg    = plot_bg
    ) |> print()
  }

  length_val
}
