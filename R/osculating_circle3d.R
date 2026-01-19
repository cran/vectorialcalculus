#' Osculating discs and circles of a spatial curve
#'
#' @description
#' For a three-dimensional parametric curve, this function constructs
#' numerical approximations to the osculating circles (and associated
#' discs) at a set of parameter values. At each requested point on the
#' curve, it approximates the Frenet frame and the curvature, and then
#' uses this information to define the center and radius of the local
#' osculating circle. Optionally, it can display these circles or discs
#' in an interactive 3D visualization using \pkg{plotly}.
#'
#' @details
#' For each parameter value in \code{t_points}, the function:
#' \itemize{
#'   \item evaluates the curve and approximates its first and second
#'         derivatives,
#'   \item constructs approximate tangent, normal and binormal directions,
#'   \item estimates the curvature from the derivative information,
#'   \item defines the center of the osculating circle by moving from the
#'         curve point along the normal direction by a distance equal to
#'         the reciprocal of the curvature,
#'   \item records the corresponding radius as that same reciprocal
#'         quantity.
#' }
#'
#' Depending on the value of \code{fill}, the function either:
#' \itemize{
#'   \item builds a filled disc that lies in the osculating plane and is
#'         bounded by the osculating circle, or
#'   \item draws only the circumference corresponding to that circle.
#' }
#'
#' A regular sampling of angles around the osculating circle is used to
#' generate the discrete representation. For filled discs, radial
#' subdivisions are added to obtain a surface mesh. The resulting objects
#' can be combined with a sampled version of the base curve and additional
#' elements such as radius segments.
#'
#' @param X,Y,Z Functions of \code{t} returning the coordinate components
#' of the curve.
#' @param a,b Numeric endpoints of the parameter interval.
#' @param t_points Numeric vector of parameter values at which osculating
#' circles or discs are constructed.
#' @param h Step size for centered finite-difference approximations.
#' @param plot Logical; if \code{TRUE}, creates a 3D visualization using
#' \pkg{plotly}.
#' @param n_samples Number of sample points used to draw the base curve
#' when \code{show_curve = TRUE}.
#' @param fill Character; either \code{"disk"} for a filled surface or
#' \code{"ring"} for the circumference only.
#' @param ru Number of radial subdivisions when drawing a filled disc.
#' @param rv Number of angular subdivisions; also used as the number of
#' points on each ring.
#' @param colorscale Character string giving the \pkg{plotly} colorscale
#' used for the discs.
#' @param opacity Numeric value between 0 and 1 controlling the opacity of
#' the discs when \code{fill = "disk"}.
#' @param ring_line List with style options for the ring when
#' \code{fill = "ring"}.
#' @param show_curve,show_points Logical values indicating whether the base
#' curve and the corresponding points on the curve should be displayed.
#' @param curve_line,point_marker Lists with \pkg{plotly} style options for
#' the base curve and the points.
#' @param show_radius Logical; if \code{TRUE}, draws a radius segment from
#' the center of each osculating circle to its boundary.
#' @param radius_phase Angle, in radians, that determines the direction of
#' the displayed radius.
#' @param radius_line List with \pkg{plotly} style options for the radius
#' segment.
#' @param scene List with 3D scene settings for the \pkg{plotly} figure.
#' @param bg List defining background colors for the figure, typically with
#' entries \code{paper} and \code{plot}.
#' @param lighting List with lighting options for \code{add_surface} when
#' \code{fill = "disk"}.
#' @param tol Numeric tolerance used in derivative-based checks and to
#' detect degenerate cases in which curvature or frame vectors cannot be
#' computed reliably.
#'
#' @return
#' A list with two components:
#' \describe{
#'   \item{\code{data}}{A tibble with columns
#'     \code{t, x, y, z, kappa, cx, cy, cz, radius,
#'     Tx, Ty, Tz, Nx, Ny, Nz, Bx, By, Bz}, containing the parameter
#'     values, the curve coordinates, the numerical curvature, the centers
#'     and radii of the osculating circles, and the associated Frenet frame
#'     vectors.}
#'   \item{\code{plot}}{A \pkg{plotly} object when \code{plot = TRUE},
#'     otherwise \code{NULL}.}
#' }
#'
#' @examples
#' X <- function(t) cos(t)
#' Y <- function(t) sin(t)
#' Z <- function(t) 0.2 * t
#' osculating_circle3d(
#'   X, Y, Z,
#'   a = 0, b = 6 * pi,
#'   t_points = c(pi, 2 * pi),
#'   plot = FALSE
#' )
#'
#' @export
osculating_circle3d <- function(
    X, Y, Z,
    a, b,
    t_points,
    h = 1e-4,
    plot = FALSE,
    n_samples = 400,
    fill = c("disk", "ring"),
    ru = 24,
    rv = 72,
    colorscale = "Reds",
    opacity = 0.6,
    ring_line = list(color = "red", width = 4, dash = "solid"),
    show_curve = TRUE,
    show_points = TRUE,
    curve_line = list(color = "blue", width = 2, dash = "solid"),
    point_marker = list(color = "black", size = 3, symbol = "circle"),
    show_radius = FALSE,
    radius_phase = 0,
    radius_line = list(color = "orange", width = 5, dash = "solid"),
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x(t)"),
      yaxis = list(title = "y(t)"),
      zaxis = list(title = "z(t)")
    ),
    bg = list(paper = "white", plot = "white"),
    lighting = list(
      ambient = 1, diffuse = 0.15, specular = 0,
      roughness = 1, fresnel = 0
    ),
    tol = 1e-10
) {
  fill <- match.arg(fill)

  if (!is.numeric(t_points) || any(!is.finite(t_points))) {
    stop("'t_points' must be finite numeric values.", call. = FALSE)
  }
  if (any(t_points < a | t_points > b)) {
    stop("All 't_points' must lie within [a, b].", call. = FALSE)
  }
  if (b < a) stop("'b' must be >= 'a'.", call. = FALSE)

  d1 <- function(f, t, h) (f(t + h) - f(t - h)) / (2 * h)
  d2 <- function(f, t, h) (f(t + h) - 2 * f(t) + f(t - h)) / (h * h)
  r  <- function(t) c(X(t), Y(t), Z(t))
  dot   <- function(a, b) sum(a * b)
  norm  <- function(a) sqrt(dot(a, a))
  cross <- function(a, b) c(
    a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1]
  )

  compute_one <- function(t0) {
    r0 <- r(t0)
    r1 <- c(d1(X, t0, h), d1(Y, t0, h), d1(Z, t0, h))
    r2 <- c(d2(X, t0, h), d2(Y, t0, h), d2(Z, t0, h))
    n1 <- norm(r1)

    if (n1 < tol) {
      return(c(
        t = t0, r0,
        kappa = NA_real_,
        rep(NA_real_, 3),     # cx, cy, cz
        radius = NA_real_,
        rep(NA_real_, 9)      # T, N, B
      ))
    }

    T <- r1 / n1
    c12 <- cross(r1, r2)
    n12 <- norm(c12)

    if (n12 < tol) {
      # curvature numerically ~ 0; radius infinite, center undefined
      return(c(
        t = t0, r0,
        kappa = 0,
        rep(NA_real_, 3),     # cx, cy, cz
        radius = Inf,
        T,
        rep(NA_real_, 6)      # N, B
      ))
    }

    B <- c12 / n12
    N <- cross(B, T)
    N <- N / norm(N)

    kappa <- n12 / (n1^3)
    C <- r0 + N / kappa
    R <- 1 / kappa

    c(
      t = t0, r0,
      kappa = kappa,
      C,
      radius = R,
      T,
      N,
      B
    )
  }

  M <- t(vapply(
    t_points,
    compute_one,
    numeric(1 + 3 + 1 + 3 + 1 + 3 + 3 + 3)
  ))
  colnames(M) <- c(
    "t", "x", "y", "z",
    "kappa",
    "cx", "cy", "cz",
    "radius",
    "Tx", "Ty", "Tz",
    "Nx", "Ny", "Nz",
    "Bx", "By", "Bz"
  )
  out <- tibble::as_tibble(M)

  fig <- NULL

  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("You need 'plotly' installed for plotting.", call. = FALSE)
    } else {
      curve_df <- NULL
      if (isTRUE(show_curve)) {
        ts <- seq(a, b, length.out = n_samples)
        curve_df <- tibble::tibble(
          x = X(ts),
          y = Y(ts),
          z = Z(ts)
        )
      }

      plt <- if (!is.null(curve_df)) {
        plotly::plot_ly(
          curve_df,
          x = ~x, y = ~y, z = ~z,
          type = "scatter3d", mode = "lines",
          line = curve_line,
          hoverinfo = "none",
          showlegend = FALSE
        )
      } else {
        plotly::plot_ly()
      }

      if (isTRUE(show_points)) {
        plt <- plt |>
          plotly::add_trace(
            x = out$x, y = out$y, z = out$z,
            type = "scatter3d", mode = "markers",
            marker = point_marker,
            showlegend = FALSE,
            hoverinfo = "none"
          )
      }

      V <- seq(0, 2 * pi, length.out = rv)
      cosV <- cos(V)
      sinV <- sin(V)

      for (i in seq_len(nrow(out))) {
        if (!is.finite(out$kappa[i]) || out$kappa[i] <= tol) next

        C <- c(out$cx[i], out$cy[i], out$cz[i])
        R <- out$radius[i]
        T <- c(out$Tx[i], out$Ty[i], out$Tz[i])
        N <- c(out$Nx[i], out$Ny[i], out$Nz[i])

        if (fill == "disk") {
          U <- seq(0, 1, length.out = ru)
          Xmat <- matrix(NA_real_, ru, rv)
          Ymat <- matrix(NA_real_, ru, rv)
          Zmat <- matrix(NA_real_, ru, rv)

          for (j in seq_len(rv)) {
            dir <- -N * cosV[j] + T * sinV[j]
            Xmat[, j] <- C[1] + R * U * dir[1]
            Ymat[, j] <- C[2] + R * U * dir[2]
            Zmat[, j] <- C[3] + R * U * dir[3]
          }

          plt <- plt |>
            plotly::add_surface(
              x = Xmat, y = Ymat, z = Zmat,
              colorscale = colorscale,
              showscale = FALSE,
              opacity = opacity,
              lighting = lighting
            )
        } else {
          x_ring <- C[1] + R * (-N[1] * cosV + T[1] * sinV)
          y_ring <- C[2] + R * (-N[2] * cosV + T[2] * sinV)
          z_ring <- C[3] + R * (-N[3] * cosV + T[3] * sinV)

          x_ring <- c(x_ring, x_ring[1])
          y_ring <- c(y_ring, y_ring[1])
          z_ring <- c(z_ring, z_ring[1])

          plt <- plt |>
            plotly::add_trace(
              x = x_ring, y = y_ring, z = z_ring,
              type = "scatter3d", mode = "lines",
              line = ring_line,
              showlegend = FALSE,
              hoverinfo = "none"
            )
        }

        if (isTRUE(show_radius)) {
          v0 <- radius_phase
          dir0 <- -N * cos(v0) + T * sin(v0)
          p_end <- C + R * dir0

          plt <- plt |>
            plotly::add_trace(
              x = c(C[1], p_end[1]),
              y = c(C[2], p_end[2]),
              z = c(C[3], p_end[3]),
              type = "scatter3d", mode = "lines",
              line = radius_line,
              showlegend = FALSE,
              hoverinfo = "none"
            )
        }
      }

      plt <- plt |>
        plotly::layout(
          title = if (fill == "disk") "Osculating Discs" else "Osculating Circles",
          scene = scene,
          paper_bgcolor = bg$paper,
          plot_bgcolor = bg$plot
        )

      fig <- plt
      print(fig)
    }
  }

  list(
    data = out,
    plot = fig
  )
}
