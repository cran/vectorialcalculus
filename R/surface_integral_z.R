#' Surface integral over a graph z = g(x, y)
#'
#' Computes numeric approximations of surface integrals over a surface given
#' in graph form \code{z = g(x, y)} on a rectangular domain in the x-y plane.
#'
#' Two types of integrals can be computed:
#' \enumerate{
#'   \item A scalar surface integral of the form \code{Integral_S phi dS}, where
#'     \code{phi(x, y, z)} is a scalar field evaluated on the surface.
#'   \item A flux (vector surface integral) of the form
#'     \code{Integral_S F dot n dS}, where \code{F(x, y, z)} is a vector field
#'     and \code{n} is the chosen unit normal direction.
#' }
#'
#' The surface is parametrized by \code{(x, y) -> (x, y, g(x, y))} over
#' a rectangular domain given by \code{xlim} and \code{ylim}. Partial
#' derivatives of \code{g} with respect to x and y are approximated by
#' finite differences on a uniform grid, and the integrals are computed
#' using a composite trapezoid rule on that grid.
#'
#' @param gfun Function of two variables \code{function(x, y)} returning
#'   the height \code{z = g(x, y)} of the surface.
#' @param xlim Numeric vector of length 2 giving the range for x,
#'   \code{c(x_min, x_max)} with \code{x_max > x_min}.
#' @param ylim Numeric vector of length 2 giving the range for y,
#'   \code{c(y_min, y_max)} with \code{y_max > y_min}.
#' @param nx Integer, number of grid points in the x direction
#'   (recommended: at least 20).
#' @param ny Integer, number of grid points in the y direction
#'   (recommended: at least 20).
#' @param scalar_phi Optional scalar field \code{function(x, y, z)}.
#'   If provided, the function computes the scalar surface integral
#'   of \code{scalar_phi} over the surface.
#' @param vector_F Optional vector field \code{function(x, y, z)} that
#'   returns a numeric vector of length 3 \code{c(Fx, Fy, Fz)}.
#'   If provided, the function computes the flux integral of \code{vector_F}
#'   across the oriented surface.
#' @param orientation Character string indicating the orientation of
#'   the normal vector, either \code{"up"} or \code{"down"}. This affects
#'   the sign of the flux integral.
#' @param plot Logical. If \code{TRUE}, returns a \pkg{plotly} surface
#'   plot colored by the available integrand (scalar field times area
#'   density or flux density).
#' @param title Character string used as the base title for the plot.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{area_density_integral}: numeric scalar with the value
#'     of the scalar surface integral (or \code{NA} if \code{scalar_phi}
#'     is not supplied).
#'   \item \code{flux_integral}: numeric scalar with the value of the
#'     flux integral (or \code{NA} if \code{vector_F} is not supplied).
#'   \item \code{plot}: a \pkg{plotly} surface object if \code{plot = TRUE},
#'     otherwise \code{NULL}.
#'   \item \code{grids}: list with matrices and grid information, including
#'     \code{X}, \code{Y}, \code{Z}, partial derivatives, and weights.
#'   \item \code{fields}: list with the scalar and/or flux integrand
#'     evaluated on the grid.
#' }
#'
#' @examples
#' \dontshow{if (interactive()) \{}
#' # Surface z = x^2 + y^2 on [-1,1] x [-1,1]
#' gfun <- function(x, y) x^2 + y^2
#'
#' # Scalar field phi(x,y,z) = 1 (surface area of the patch)
#' phi <- function(x, y, z) 1
#'
#' # Vector field F = (0, 0, 1), flux through the surface
#' Fvec <- function(x, y, z) c(0, 0, 1)
#'
#' res <- surface_integral_z(
#'   gfun,
#'   xlim = c(-1, 1),
#'   ylim = c(-1, 1),
#'   nx = 60, ny = 60,
#'   scalar_phi = phi,
#'   vector_F = Fvec,
#'   orientation = "up",
#'   plot = FALSE
#' )
#' res$area_density_integral
#' res$flux_integral
#' \dontshow{\}}
#'
#' @export
surface_integral_z <- function(
    gfun,
    xlim, ylim,
    nx = 160, ny = 160,
    scalar_phi = NULL,
    vector_F  = NULL,
    orientation = c("up", "down"),
    plot = TRUE,
    title = "Surface integral over z = g(x,y)"
) {
  # --- basic checks (avoid cryptic errors) ---
  if (!is.function(gfun)) {
    stop("'gfun' must be a function of two arguments (x, y).", call. = FALSE)
  }
  if (!is.numeric(xlim) || length(xlim) != 2L || !all(is.finite(xlim))) {
    stop("'xlim' must be a numeric vector c(x_min, x_max) with finite entries.", call. = FALSE)
  }
  if (!is.numeric(ylim) || length(ylim) != 2L || !all(is.finite(ylim))) {
    stop("'ylim' must be a numeric vector c(y_min, y_max) with finite entries.", call. = FALSE)
  }
  if (xlim[2] <= xlim[1]) {
    stop("'xlim' must satisfy xlim[2] > xlim[1].", call. = FALSE)
  }
  if (ylim[2] <= ylim[1]) {
    stop("'ylim' must satisfy ylim[2] > ylim[1].", call. = FALSE)
  }
  if (!is.numeric(nx) || length(nx) != 1L || !is.finite(nx) || nx < 2 ||
      abs(nx - round(nx)) > .Machine$double.eps^0.5) {
    stop("'nx' must be an integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(ny) || length(ny) != 1L || !is.finite(ny) || ny < 2 ||
      abs(ny - round(ny)) > .Machine$double.eps^0.5) {
    stop("'ny' must be an integer >= 2.", call. = FALSE)
  }

  orientation <- match.arg(orientation)

  if (isTRUE(plot) && !requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for plotting. Install it with install.packages('plotly').",
         call. = FALSE)
  }

  nx <- as.integer(nx)
  ny <- as.integer(ny)

  # --- grid in (x,y) ---
  x <- seq(xlim[1], xlim[2], length.out = nx)
  y <- seq(ylim[1], ylim[2], length.out = ny)
  X <- matrix(rep(x, each = ny), nrow = ny, ncol = nx)
  Y <- matrix(rep(y, times = nx), nrow = ny, ncol = nx)

  # Evaluate surface z = g(x,y)
  gv <- function(xx, yy) vapply(seq_along(xx), function(i) gfun(xx[i], yy[i]), numeric(1))
  Z  <- matrix(gv(as.vector(X), as.vector(Y)), nrow = ny, ncol = nx)

  # --- finite differences for gx, gy (central inside, one-sided on edges) ---
  d_dx <- function(M, dx) {
    out <- M * 0
    nc <- ncol(M)
    if (nc >= 3L) {
      out[, 2:(nc - 1L)] <- (M[, 3:nc] - M[, 1:(nc - 2L)]) / (2 * dx)
    }
    out[, 1L]   <- (M[, 2L] - M[, 1L]) / dx
    out[, nc]   <- (M[, nc] - M[, nc - 1L]) / dx
    out
  }
  d_dy <- function(M, dy) {
    out <- M * 0
    nr <- nrow(M)
    if (nr >= 3L) {
      out[2:(nr - 1L), ] <- (M[3:nr, ] - M[1:(nr - 2L), ]) / (2 * dy)
    }
    out[1L, ]   <- (M[2L, ] - M[1L, ]) / dy
    out[nr, ]   <- (M[nr, ] - M[nr - 1L, ]) / dy
    out
  }

  dx <- (xlim[2] - xlim[1]) / (nx - 1L)
  dy <- (ylim[2] - ylim[1]) / (ny - 1L)
  gx <- d_dx(Z, dx)
  gy <- d_dy(Z, dy)

  # Composite trapezoid weights on the (x,y) grid
  wx <- rep(1, nx); wx[c(1L, nx)] <- 0.5
  wy <- rep(1, ny); wy[c(1L, ny)] <- 0.5
  W  <- outer(wy, wx)  # ny x nx

  # --- 1) Scalar surface integral: Integral_S phi dS ---
  #     Integral_D phi(x,y,g) * sqrt(1 + gx^2 + gy^2) dx dy
  area_density_integral <- NA_real_
  scalar_integrand <- NULL
  if (is.function(scalar_phi)) {
    phiv <- function(xx, yy, zz) {
      vapply(seq_along(xx),
             function(i) scalar_phi(xx[i], yy[i], zz[i]),
             numeric(1))
    }
    PHI <- matrix(
      phiv(as.vector(X), as.vector(Y), as.vector(Z)),
      nrow = ny, ncol = nx
    )
    J <- sqrt(1 + gx^2 + gy^2)
    scalar_integrand <- PHI * J
    area_density_integral <- sum(W * scalar_integrand) * dx * dy
  }

  # --- 2) Flux integral: Integral_S F dot n dS ---
  #     For upward orientation: n non-unit is (-gx, -gy, 1)
  #     For downward orientation: (gx, gy, -1)
  flux_integral  <- NA_real_
  flux_integrand <- NULL
  if (is.function(vector_F)) {
    Fv <- function(xx, yy, zz) {
      val <- vector_F(xx, yy, zz)
      if (!is.numeric(val) || length(val) != 3L) {
        stop("vector_F must return numeric length-3 c(Fx, Fy, Fz).", call. = FALSE)
      }
      val
    }

    Fx <- matrix(0, ny, nx)
    Fy <- matrix(0, ny, nx)
    Fz <- matrix(0, ny, nx)
    for (j in seq_len(nx)) {
      for (i in seq_len(ny)) {
        vv <- Fv(X[i, j], Y[i, j], Z[i, j])
        Fx[i, j] <- vv[1L]
        Fy[i, j] <- vv[2L]
        Fz[i, j] <- vv[3L]
      }
    }

    if (orientation == "up") {
      Nx <- -gx; Ny <- -gy; Nz <-  1
    } else {
      Nx <-  gx; Ny <-  gy; Nz <- -1
    }

    flux_integrand <- Fx * Nx + Fy * Ny + Fz * Nz
    flux_integral  <- sum(W * flux_integrand) * dx * dy
  }

  # --- Optional plot: color by scalar integrand if available, else by flux ---
  plt <- NULL
  if (isTRUE(plot)) {
    col_mat <- if (!is.null(scalar_integrand)) scalar_integrand else flux_integrand
    if (is.null(col_mat)) col_mat <- Z * 0

    colorbar_title <- if (!is.null(scalar_integrand)) {
      "phi * sqrt(1 + gx^2 + gy^2)"
    } else if (!is.null(flux_integrand)) {
      "F dot N"
    } else {
      "z"
    }

    main_title <- title
    if (!is.na(area_density_integral)) {
      main_title <- paste0(
        main_title,
        sprintf(" | scalar integral approx %.6g", area_density_integral)
      )
    }
    if (!is.na(flux_integral)) {
      main_title <- paste0(
        main_title,
        sprintf(" | flux integral approx %.6g", flux_integral)
      )
    }

    plt <- plotly::plot_ly() |>
      plotly::add_surface(
        x = X, y = Y, z = Z, surfacecolor = col_mat,
        colorbar = list(title = colorbar_title),
        name = "Surface"
      ) |>
      plotly::layout(
        title = main_title,
        scene = list(
          xaxis = list(title = "x"),
          yaxis = list(title = "y"),
          zaxis = list(title = "z")
        )
      )
  }

  list(
    area_density_integral = area_density_integral,
    flux_integral         = flux_integral,
    plot                  = plt,
    grids = list(
      X = X, Y = Y, Z = Z,
      gx = gx, gy = gy,
      W = W, dx = dx, dy = dy
    ),
    fields = list(
      scalar_integrand = scalar_integrand,
      flux_integrand   = flux_integrand
    )
  )
}
