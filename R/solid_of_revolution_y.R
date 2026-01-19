#' Solid of revolution around a horizontal line
#'
#' Construct a three-dimensional surface for the solid obtained by rotating
#' the graph of a function \code{f(x)} around the line \code{y = a} on a
#' finite interval, and compute its volume and surface areas.
#'
#' @param f Function of one numeric argument \code{x} that returns a numeric value.
#' @param xlim Numeric vector of length two with the limits of the x interval
#'   \code{c(xmin, xmax)}, with \code{xmax > xmin}.
#' @param a Numeric scalar that gives the horizontal axis of rotation.
#' @param nx Integer number of grid points along the x direction (for plotting).
#' @param nt Integer number of grid points along the angular direction (for plotting).
#' @param deriv Optional function that returns the derivative of \code{f}. If \code{NULL}
#'   a numeric derivative is used.
#' @param h Optional numeric step for the numeric derivative.
#' @param include_end_caps Logical value. If \code{TRUE} the area of the circular
#'   cross sections at the ends of the interval is added to the total area.
#' @param plot Logical; if \code{TRUE} and \pkg{plotly} is available, the solid is drawn.
#' @param colors List with optional entries \code{surface}, \code{axis} and \code{curve} that
#'   control the colours used in the plot.
#' @param opacity Numeric value between 0 and 1 that controls the surface
#'   opacity in the plot.
#' @param show_axis Logical value indicating whether the axis of rotation
#'   is drawn.
#' @param show_profile_curve Logical value indicating whether the generating
#'   curve is drawn on the surface.
#' @param curve_thetas Numeric vector of angles (in radians) where profile curves are
#'   drawn.
#' @param curve_width Numeric line width for the profile curves.
#' @param curve_opacity Numeric value between 0 and 1 for the profile curves.
#' @param scene List of \pkg{plotly} scene options used in \code{plotly::layout()}.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{volume}: numeric value of the volume.
#'   \item \code{surface_area_lateral}: numeric value of the lateral area.
#'   \item \code{surface_area_total}: numeric value of the total area.
#'   \item \code{figure}: \pkg{plotly} object with the three-dimensional plot
#'     if \code{plot = TRUE} and \pkg{plotly} is available; otherwise \code{NULL}.
#' }
#'
#' @examples
#' f <- function(x) sqrt(x)
#' solid_of_revolution_y(f, xlim = c(0, 4), a = 0, plot = FALSE)
#'
#' @export
solid_of_revolution_y <- function(
    f,
    xlim,
    a,
    nx = 120L,
    nt = 120L,
    deriv = NULL,
    h = NULL,
    include_end_caps = FALSE,
    plot = TRUE,
    colors = list(surface = "steelblue", axis = "black", curve = "firebrick"),
    opacity = 0.9,
    show_axis = TRUE,
    show_profile_curve = TRUE,
    curve_thetas = 0,
    curve_width = 4,
    curve_opacity = 1,
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z")
    )
){
  stopifnot(is.function(f))
  if (length(xlim) != 2L || !is.finite(xlim[1]) || !is.finite(xlim[2]) || xlim[2] <= xlim[1]) {
    stop("'xlim' must be c(xmin, xmax) with finite numbers and xmax > xmin.", call. = FALSE)
  }
  if (!is.numeric(a) || length(a) != 1L || !is.finite(a)) {
    stop("'a' must be a finite numeric scalar.", call. = FALSE)
  }

  nx <- as.integer(nx)
  nt <- as.integer(nt)
  if (!is.finite(nx) || !is.finite(nt) || nx < 2L || nt < 2L) {
    stop("'nx' and 'nt' must be integer values >= 2.", call. = FALSE)
  }

  # Vectorized wrappers (in case f/deriv are not vectorized)
  f_vec <- function(x) vapply(x, f, numeric(1))
  if (is.null(deriv)) {
    d1 <- function(x) {
      hx0 <- if (is.null(h)) 1e-5 else as.numeric(h)
      g <- function(x0) {
        hx <- if (is.null(h)) 1e-5 * (1 + abs(x0)) else hx0
        (f(x0 + hx) - f(x0 - hx)) / (2 * hx)
      }
      vapply(x, g, numeric(1))
    }
    deriv_vec <- d1
  } else {
    stopifnot(is.function(deriv))
    deriv_vec <- function(x) vapply(x, deriv, numeric(1))
  }

  # --- Integrands
  vol_integrand <- function(x) {
    r <- f_vec(x) - a
    (r * r) * pi
  }
  area_integrand <- function(x) {
    r  <- abs(f_vec(x) - a)
    fp <- deriv_vec(x)
    2 * pi * r * sqrt(1 + fp * fp)
  }

  V <- try(stats::integrate(vol_integrand, lower = xlim[1], upper = xlim[2])$value, silent = TRUE)
  if (inherits(V, "try-error")) {
    stop("Volume integral failed; check that f(x) is integrable on 'xlim'.", call. = FALSE)
  }

  S_lat <- try(stats::integrate(area_integrand, lower = xlim[1], upper = xlim[2])$value, silent = TRUE)
  if (inherits(S_lat, "try-error")) {
    stop("Surface-area integral failed; check f and its derivative on 'xlim'.", call. = FALSE)
  }

  S_tot <- S_lat
  if (isTRUE(include_end_caps)) {
    r0 <- f_vec(xlim[1]) - a
    r1 <- f_vec(xlim[2]) - a
    S_tot <- S_lat + pi * (r0 * r0 + r1 * r1)
  }

  # --- Plot
  fig <- NULL
  if (isTRUE(plot)) {
    have_plotly <- requireNamespace("plotly", quietly = TRUE)
    if (!have_plotly) {
      warning("Package 'plotly' not found; returning numeric results only.", call. = FALSE)
    } else {
      xs     <- seq(xlim[1], xlim[2], length.out = nx)
      thetas <- seq(0, 2 * pi, length.out = nt)
      r      <- f_vec(xs) - a

      # Parametric surface matrices (nt x nx)
      X <- matrix(rep(xs, each = nt), nrow = nt, ncol = nx)
      Y <- matrix(a, nrow = nt, ncol = nx) + outer(cos(thetas), r, `*`)
      Z <- outer(sin(thetas), r, `*`)

      p <- plotly::plot_ly()
      p <- plotly::add_surface(
        p,
        x = X, y = Y, z = Z,
        showscale = FALSE,
        opacity = opacity,
        colorscale = list(
          list(0, colors$surface %||% "steelblue"),
          list(1, colors$surface %||% "steelblue")
        ),
        name = "solid"
      )

      if (isTRUE(show_axis)) {
        p <- plotly::add_trace(
          p,
          x = xlim, y = c(a, a), z = c(0, 0),
          type = "scatter3d", mode = "lines",
          line = list(color = colors$axis %||% "black", width = 4),
          name = "axis y = a",
          showlegend = TRUE,
          inherit = FALSE
        )
      }

      # Profile curve(s) mapped onto the surface
      if (isTRUE(show_profile_curve) && length(curve_thetas) > 0) {
        xsc <- seq(xlim[1], xlim[2], length.out = max(400L, nx))
        rc  <- f_vec(xsc) - a
        for (th in as.numeric(curve_thetas)) {
          yc <- a + rc * cos(th)
          zc <-       rc * sin(th)
          p <- plotly::add_trace(
            p, x = xsc, y = yc, z = zc,
            type = "scatter3d", mode = "lines",
            line = list(color = colors$curve %||% "firebrick", width = curve_width),
            opacity = curve_opacity,
            name = sprintf("profile f(x) at angle %.2f rad", th),
            showlegend = TRUE,
            inherit = FALSE
          )
        }
      }

      p <- plotly::layout(
        p,
        title = "Solid of revolution about y = a",
        scene = scene
      )
      fig <- p
      print(fig)
    }
  }

  list(
    volume = V,
    surface_area_lateral = S_lat,
    surface_area_total = S_tot,
    figure = fig
  )
}
