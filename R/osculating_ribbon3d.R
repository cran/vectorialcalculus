#' Osculating ribbon along a 3D parametric curve
#'
#' @description
#' Constructs a narrow ribbon that follows a three-dimensional parametric
#' curve. The ribbon is based on the Frenet-Serret frame of the curve,
#' computed numerically along a set of sample points. The ribbon extends a
#' small distance in the normal direction of the curve, producing a thin
#' band that helps visualize how the curve bends and twists in space.
#'
#' @details
#' The function samples the curve at \code{n_t} points and computes the
#' numerical tangent, normal and binormal directions using finite-difference
#' approximations of the derivatives. At each sampled point, a short segment
#' is taken in the normal direction to define the width of the ribbon. These
#' segments are interpolated across the curve and subdivided according to
#' \code{n_u} to produce a mesh that represents the ribbon surface.
#'
#' Optionally, the function can display the ribbon in an interactive 3D plot
#' using \pkg{plotly}. The base curve, the centerline, and optional grid
#' lines on the ribbon surface can be shown or hidden independently.
#'
#' @param X,Y,Z Functions returning the coordinate components of the curve as
#' functions of the parameter \code{t}.
#' @param a,b Numeric values giving the endpoints of the parameter interval.
#' @param h Step size used in the finite-difference approximations.
#' @param plot Logical; if \code{TRUE}, produces a 3D visualization of the
#' ribbon using \pkg{plotly}.
#' @param n_t Number of sample points along the curve.
#' @param n_u Number of subdivisions across the width of the ribbon.
#' @param u_max Half-width of the ribbon, measured in units of the normal
#' direction.
#' @param colorscale Character string giving the \pkg{plotly} colorscale used
#' for the ribbon surface.
#' @param opacity Numeric value between 0 and 1 controlling the transparency
#' of the ribbon.
#' @param show_curve Logical; if \code{TRUE}, draws the base curve.
#' @param show_centers Logical; if \code{TRUE}, draws the centerline joining
#' the midpoints of the ribbon cross-sections.
#' @param curve_line List with \pkg{plotly} style options for the base curve.
#' @param centers_line List with \pkg{plotly} style options for the
#' centerline.
#' @param show_surface_grid Logical; if \code{TRUE}, draws a grid on the
#' surface of the ribbon.
#' @param surface_grid_color Character string giving the color of the grid
#' lines on the ribbon.
#' @param surface_grid_width Numeric value giving the width of the surface
#' grid lines.
#' @param show_axis_grid Logical; if \code{TRUE}, displays gridlines on the
#' coordinate axes in the \pkg{plotly} scene.
#' @param scene List with 3D scene settings for the \pkg{plotly} figure.
#' @param bg List defining the background colors of the figure, typically
#' with entries \code{paper} and \code{plot}.
#' @param lighting List with lighting settings for the surface in
#' \pkg{plotly}.
#' @param tol Numeric tolerance used to detect numerical instabilities when
#' computing the derivative-based frame vectors.
#'
#' @return
#' A list with two components:
#' \describe{
#'   \item{\code{data}}{A tibble containing the sampled parameter
#'     values, the coordinates of the curve, and the corresponding tangent,
#'     normal and binormal directions.}
#'   \item{\code{plot}}{A \pkg{plotly} object if \code{plot = TRUE}, otherwise
#'     \code{NULL}.}
#' }
#'
#' @examples
#' X <- function(t) cos(t)
#' Y <- function(t) sin(t)
#' Z <- function(t) 0.2 * t
#' osculating_ribbon3d(X, Y, Z, a = 0, b = 4*pi, plot = FALSE)
#'
#' @export
osculating_ribbon3d <- function(
    X, Y, Z,
    a, b,
    h = 1e-4,
    plot = FALSE,
    n_t = 400,
    n_u = 25,
    u_max = 1,
    colorscale = "Blues",
    opacity = 0.35,
    show_curve = TRUE,
    show_centers = TRUE,
    curve_line   = list(color = "black", width = 2),
    centers_line = list(color = "red",   width = 2),
    show_surface_grid = TRUE,
    surface_grid_color = "rgba(60,80,200,0.25)",
    surface_grid_width = 1,
    show_axis_grid = FALSE,
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
  if (b < a) stop("'b' must be >= 'a'.", call. = FALSE)
  if (!(u_max > 0 && u_max <= 1)) {
    stop("'u_max' must satisfy 0 < u_max <= 1.", call. = FALSE)
  }

  # Derivatives and utilities
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

  t_seq <- seq(a, b, length.out = n_t)
  u_seq <- seq(0, u_max, length.out = n_u)

  # Surface matrices (n_u x n_t)
  Xmat <- matrix(NA_real_, n_u, n_t)
  Ymat <- matrix(NA_real_, n_u, n_t)
  Zmat <- matrix(NA_real_, n_u, n_t)

  # For curve and centers (for traces)
  curve_df   <- data.frame(
    t = t_seq,
    x = NA_real_, y = NA_real_, z = NA_real_
  )
  centers_df <- data.frame(
    t = t_seq,
    cx = NA_real_, cy = NA_real_, cz = NA_real_
  )

  # To store Frenet frame
  Tmat <- matrix(NA_real_, n_t, 3)
  Nmat <- matrix(NA_real_, n_t, 3)
  Bmat <- matrix(NA_real_, n_t, 3)

  for (j in seq_along(t_seq)) {
    t0 <- t_seq[j]
    r0 <- r(t0)
    r1 <- c(d1(X, t0, h), d1(Y, t0, h), d1(Z, t0, h))
    r2 <- c(d2(X, t0, h), d2(Y, t0, h), d2(Z, t0, h))

    n1 <- norm(r1)
    if (n1 < tol) next

    T <- r1 / n1
    c12 <- cross(r1, r2)
    n12 <- norm(c12)
    if (n12 < tol) next

    B <- c12 / n12
    N <- cross(B, T)
    N <- N / norm(N)

    kappa <- n12 / (n1^3)
    if (kappa <= tol) next

    R <- 1 / kappa
    C <- r0 + N * R

    # Save curve point and center
    curve_df$x[j] <- r0[1]
    curve_df$y[j] <- r0[2]
    curve_df$z[j] <- r0[3]

    centers_df$cx[j] <- C[1]
    centers_df$cy[j] <- C[2]
    centers_df$cz[j] <- C[3]

    # Save Frenet frame
    Tmat[j, ] <- T
    Nmat[j, ] <- N
    Bmat[j, ] <- B

    # Ribbon surface column: u in [0, u_max] along N
    Xmat[, j] <- r0[1] + u_seq * R * N[1]
    Ymat[, j] <- r0[2] + u_seq * R * N[2]
    Zmat[, j] <- r0[3] + u_seq * R * N[3]
  }

  fig <- NULL

  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("To plot you need 'plotly' installed.", call. = FALSE)
    } else {
      contours_arg <- if (isTRUE(show_surface_grid)) list(
        x = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
        y = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
        z = list(show = FALSE)
      ) else NULL

      plt <- plotly::plot_ly() |>
        plotly::add_surface(
          x = Xmat, y = Ymat, z = Zmat,
          colorscale = colorscale,
          showscale  = FALSE,
          opacity    = opacity,
          lighting   = lighting,
          contours   = contours_arg
        )

      if (isTRUE(show_curve) && any(is.finite(curve_df$x))) {
        plt <- plt |>
          plotly::add_trace(
            data = curve_df,
            x = ~x, y = ~y, z = ~z,
            type = "scatter3d", mode = "lines",
            line = curve_line,
            hoverinfo = "none",
            showlegend = FALSE
          )
      }

      if (isTRUE(show_centers) && any(is.finite(centers_df$cx))) {
        plt <- plt |>
          plotly::add_trace(
            x = centers_df$cx,
            y = centers_df$cy,
            z = centers_df$cz,
            type = "scatter3d", mode = "lines",
            line = centers_line,
            hoverinfo = "none",
            showlegend = FALSE
          )
      }

      scene_final <- scene
      for (ax in c("xaxis", "yaxis", "zaxis")) {
        if (is.null(scene_final[[ax]])) scene_final[[ax]] <- list()
        scene_final[[ax]]$showgrid <- isTRUE(show_axis_grid)
      }

      plt <- plt |>
        plotly::layout(
          title = sprintf("Osculating ribbon (u_max = %.2f)", u_max),
          scene = scene_final,
          paper_bgcolor = bg$paper,
          plot_bgcolor  = bg$plot
        )

      fig <- plt
      print(fig)
    }
  }

  data_out <- tibble::tibble(
    t  = t_seq,
    x  = curve_df$x,
    y  = curve_df$y,
    z  = curve_df$z,
    Tx = Tmat[, 1],
    Ty = Tmat[, 2],
    Tz = Tmat[, 3],
    Nx = Nmat[, 1],
    Ny = Nmat[, 2],
    Nz = Nmat[, 3],
    Bx = Bmat[, 1],
    By = Bmat[, 2],
    Bz = Bmat[, 3]
  )

  list(
    data = data_out,
    plot = fig
  )
}
