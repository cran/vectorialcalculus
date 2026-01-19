#' Plot a parametric surface and estimate its area
#'
#' This function plots a smooth parametric surface defined by
#' functions x(u,v), y(u,v), and z(u,v). It also estimates the
#' surface area using a trapezoidal approximation based on the
#' magnitudes of partial-derivative cross products.
#'
#' The parametric domain is given by ranges for u and v, and
#' the surface is evaluated on a regular grid with sizes
#' specified by \code{nu} and \code{nv}. Finite differences are used to
#' approximate partial derivatives with respect to u and v.
#'
#' @param xfun Function of two arguments (u, v) returning the
#'   x-coordinate of the surface.
#' @param yfun Function of two arguments (u, v) returning the
#'   y-coordinate of the surface.
#' @param zfun Function of two arguments (u, v) returning the
#'   z-coordinate of the surface.
#' @param urange Numeric vector of length 2 giving the interval
#'   for the parameter u, \code{c(u_min, u_max)} with \code{u_max > u_min}.
#' @param vrange Numeric vector of length 2 giving the interval
#'   for the parameter v, \code{c(v_min, v_max)} with \code{v_max > v_min}.
#' @param nu Integer, number of grid points along u (recommended:
#'   at least 20 for a reasonable surface).
#' @param nv Integer, number of grid points along v (recommended:
#'   at least 20).
#' @param h_u Numeric step size for finite differences in u.
#'   If \code{NULL}, a default based on the grid spacing is used.
#' @param h_v Numeric step size for finite differences in v.
#'   If \code{NULL}, a default based on the grid spacing is used.
#' @param title_prefix Character string used in the plot title.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{plot}: a \pkg{plotly} surface object showing the
#'     parametric surface.
#'   \item \code{area}: numeric estimate of the surface area.
#'   \item \code{grid}: list with elements \code{U}, \code{V}, \code{X},
#'     \code{Y}, and \code{Z}, representing the parameter values and
#'     evaluated surface.
#' }
#'
#' @examples
#' \dontshow{if (interactive()) \{}
#' # Example: torus-like parametric surface
#' xfun <- function(u, v) (2 + cos(v)) * cos(u)
#' yfun <- function(u, v) (2 + cos(v)) * sin(u)
#' zfun <- function(u, v) sin(v)
#' result <- surface_parametric_area(
#'   xfun, yfun, zfun,
#'   urange = c(0, 2*pi),
#'   vrange = c(0, 2*pi),
#'   nu = 80, nv = 80
#' )
#' result$area
#' \dontshow{\}}
#'
#' @export
surface_parametric_area <- function(
    xfun, yfun, zfun,
    urange = c(0, 2*pi),
    vrange = c(0, 2*pi),
    nu = 160, nv = 160,
    h_u = NULL, h_v = NULL,
    title_prefix = "r(u,v)"
) {
  # --- basic checks ---
  if (!is.function(xfun) || !is.function(yfun) || !is.function(zfun)) {
    stop("'xfun', 'yfun' and 'zfun' must be functions of two arguments (u, v).",
         call. = FALSE)
  }
  if (!is.numeric(urange) || length(urange) != 2L || any(!is.finite(urange))) {
    stop("'urange' must be numeric c(u_min, u_max) with finite entries.", call. = FALSE)
  }
  if (!is.numeric(vrange) || length(vrange) != 2L || any(!is.finite(vrange))) {
    stop("'vrange' must be numeric c(v_min, v_max) with finite entries.", call. = FALSE)
  }
  if (urange[2] <= urange[1]) {
    stop("'urange' must satisfy urange[2] > urange[1].", call. = FALSE)
  }
  if (vrange[2] <= vrange[1]) {
    stop("'vrange' must satisfy vrange[2] > vrange[1].", call. = FALSE)
  }

  if (!is.numeric(nu) || length(nu) != 1L || !is.finite(nu) ||
      nu < 2 || abs(nu - round(nu)) > .Machine$double.eps^0.5) {
    stop("'nu' must be an integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(nv) || length(nv) != 1L || !is.finite(nv) ||
      nv < 2 || abs(nv - round(nv)) > .Machine$double.eps^0.5) {
    stop("'nv' must be an integer >= 2.", call. = FALSE)
  }

  nu <- as.integer(nu)
  nv <- as.integer(nv)

  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required. Install it with install.packages('plotly').",
         call. = FALSE)
  }

  # --- parameter grid ---
  u <- seq(urange[1], urange[2], length.out = nu)
  v <- seq(vrange[1], vrange[2], length.out = nv)
  U <- matrix(rep(u, each = nv), nrow = nv, ncol = nu)
  V <- matrix(rep(v, times = nu), nrow = nv, ncol = nu)

  # vectorized coordinate evaluation (safe over vectors)
  xv <- function(uu, vv) {
    vapply(seq_along(uu), function(i) xfun(uu[i], vv[i]), numeric(1))
  }
  yv <- function(uu, vv) {
    vapply(seq_along(uu), function(i) yfun(uu[i], vv[i]), numeric(1))
  }
  zv <- function(uu, vv) {
    vapply(seq_along(uu), function(i) zfun(uu[i], vv[i]), numeric(1))
  }

  # Evaluate and force matrix shape (nv x nu)
  X <- matrix(xv(as.vector(U), as.vector(V)), nrow = nv, ncol = nu)
  Y <- matrix(yv(as.vector(U), as.vector(V)), nrow = nv, ncol = nu)
  Z <- matrix(zv(as.vector(U), as.vector(V)), nrow = nv, ncol = nu)

  # --- finite-difference steps in parameter space ---
  if (is.null(h_u)) {
    h_u <- max(1e-6, (urange[2] - urange[1]) / max(1L, nu - 1L))
  }
  if (is.null(h_v)) {
    h_v <- max(1e-6, (vrange[2] - vrange[1]) / max(1L, nv - 1L))
  }

  d_du <- function(M, delta) {
    if (is.null(dim(M))) stop("Expected a matrix in d_du; got a vector.", call. = FALSE)
    out <- M * 0
    nc <- ncol(M)
    if (nc >= 3L) {
      out[, 2:(nc - 1L)] <- (M[, 3:nc] - M[, 1:(nc - 2L)]) / (2 * delta)
    }
    out[, 1L] <- (M[, 2L] - M[, 1L]) / delta
    out[, nc] <- (M[, nc] - M[, nc - 1L]) / delta
    out
  }
  d_dv <- function(M, delta) {
    if (is.null(dim(M))) stop("Expected a matrix in d_dv; got a vector.", call. = FALSE)
    out <- M * 0
    nr <- nrow(M)
    if (nr >= 3L) {
      out[2:(nr - 1L), ] <- (M[3:nr, ] - M[1:(nr - 2L), ]) / (2 * delta)
    }
    out[1L, ] <- (M[2L, ] - M[1L, ]) / delta
    out[nr, ] <- (M[nr, ] - M[nr - 1L, ]) / delta
    out
  }

  # --- partial derivatives r_u and r_v ---
  Xu <- d_du(X, h_u); Yu <- d_du(Y, h_u); Zu <- d_du(Z, h_u)
  Xv <- d_dv(X, h_v); Yv <- d_dv(Y, h_v); Zv <- d_dv(Z, h_v)

  # cross product r_u x r_v and Jacobian magnitude
  Cx <- Yu * Zv - Zu * Yv
  Cy <- Zu * Xv - Xu * Zv
  Cz <- Xu * Yv - Yu * Xv
  mag <- sqrt(Cx^2 + Cy^2 + Cz^2)

  # --- composite trapezoid weights in (u,v) ---
  w_u <- rep(1, nu); w_u[c(1L, nu)] <- 0.5
  w_v <- rep(1, nv); w_v[c(1L, nv)] <- 0.5
  W   <- outer(w_v, w_u)

  du <- (urange[2] - urange[1]) / max(1L, nu - 1L)
  dv <- (vrange[2] - vrange[1]) / max(1L, nv - 1L)

  area_est <- sum(W * mag) * du * dv

  # --- plot ---
  main_title <- paste0(
    title_prefix, " | surface area (approx): ",
    formatC(area_est, digits = 8, format = "g")
  )

  p <- plotly::plot_ly() |>
    plotly::add_surface(
      x = X, y = Y, z = Z,
      showscale = FALSE,
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

  message(sprintf("Estimated surface area (approx): %.10g", area_est))

  list(
    plot = p,
    area = area_est,
    grid = list(U = U, V = V, X = X, Y = Y, Z = Z)
  )
}
