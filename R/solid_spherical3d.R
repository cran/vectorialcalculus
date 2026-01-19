#' Solid in spherical coordinates with Plotly visualization and volume
#'
#' @description
#' Draws a three-dimensional solid described in spherical coordinates by:
#' \itemize{
#'   \item a radial variable \code{r} between \code{R1(theta, phi)} and
#'         \code{R2(theta, phi)},
#'   \item an azimuthal angle \code{theta} in the interval
#'         \code{[theta_range[1], theta_range[2]]} (in radians),
#'   \item a polar angle \code{phi} in the interval
#'         \code{[phi_range[1], phi_range[2]]} (in radians).
#' }
#'
#' The function uses the standard convention for spherical coordinates:
#' \code{theta} is the azimuth (angle in the xy-plane) and \code{phi} is the
#' polar angle measured from the positive z-axis.
#'
#' Optionally, the volume of the solid is computed using the spherical volume
#' element. The exact integral has the form:
#' \itemize{
#'   \item inner integral: from \code{r = R1(theta, phi)} to \code{r = R2(theta, phi)}
#'         of \code{r^2 * sin(phi) dr},
#'   \item outer integrals: over \code{phi} in \code{[phi_min, phi_max]} and
#'         \code{theta} in \code{[theta_min, theta_max]}.
#' }
#' Equivalently, for each pair \code{(theta, phi)} one integrates
#' \code{(R2(theta, phi)^3 - R1(theta, phi)^3) / 3 * sin(phi)} over the angular
#' rectangle.
#'
#' @param R1,R2 Functions \code{function(theta, phi)} giving the inner and outer
#'        radius, respectively, as numeric scalars.
#' @param theta_range Numeric vector of length 2, \code{c(theta_min, theta_max)},
#'        giving the azimuth interval in radians.
#' @param phi_range   Numeric vector of length 2, \code{c(phi_min, phi_max)},
#'        giving the polar angle interval in radians.
#' @param n_theta,n_phi Mesh resolution for the two boundary surfaces. Each
#'        surface is sampled on an \code{n_phi x n_theta} grid.
#' @param plot Logical. If \code{TRUE}, the solid boundaries are drawn with
#'        \pkg{plotly}.
#' @param show_surfaces Logical vector of length 2 indicating which spherical
#'        shells to show in the plot, in the order \code{c(r = R1, r = R2)}.
#' @param colorscales Colorscales for the two surfaces. You can pass:
#'        \itemize{
#'          \item a single Plotly colorscale (string, single color, or vector
#'                of colors) applied to both surfaces, or
#'          \item a list of length 2, with one colorscale per surface.
#'        }
#'        Flat colors such as \code{"#2a9d8f"} or \code{"rgba(0,0,0,0.6)"} are
#'        also accepted.
#' @param opacities Numeric vector of length 1 or 2 giving the opacity of the
#'        two surfaces.
#' @param show_surface_grid Logical. If \code{TRUE}, draws grid lines on the
#'        surfaces.
#' @param surface_grid_color,surface_grid_width Color and width for the surface
#'        grid lines.
#' @param scene Plotly 3D scene settings. By default the aspect mode is
#'        \code{"data"} and each axis has a title.
#' @param bg Background colors, typically a list of the form
#'        \code{list(paper = "white", plot = "white")}.
#' @param compute_volume Logical. If \code{TRUE}, the volume of the solid is
#'        approximated numerically using the spherical volume integral.
#' @param vol_method Character string indicating the integration method for
#'        the volume: \code{"adaptive"} uses nested \code{stats::integrate},
#'        while \code{"grid"} uses the trapezoidal rule on a regular grid.
#' @param n_th_vol,n_ph_vol Integer resolutions for the grid method in the
#'        azimuth and polar directions, respectively (ignored when
#'        \code{vol_method = "adaptive"}).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{theta_seq}, \code{phi_seq}: the parameter sequences used for
#'         plotting,
#'   \item \code{R1_surf}, \code{R2_surf}: lists containing the matrices
#'         \code{X}, \code{Y}, \code{Z} for the two boundary surfaces
#'         (or \code{NULL} if the corresponding surface is not shown),
#'   \item \code{fig}: a \pkg{plotly} figure when \code{plot = TRUE},
#'         otherwise \code{NULL},
#'   \item \code{volume}: \code{NULL} if \code{compute_volume = FALSE},
#'         or a list with the numeric volume estimate, the method used and
#'         additional metadata when \code{compute_volume = TRUE}.
#' }
#'
#' @examples
#' \dontshow{if (interactive()) \{}
#' # Example 1: Spherical shell: a <= r <= b, independent of angles
#' R1 <- function(th, ph) 0.8
#' R2 <- function(th, ph) 1.2
#' out <- solid_spherical3d(
#'   R1, R2,
#'   theta_range = c(0, 2*pi),
#'   phi_range   = c(0, pi),
#'   plot = TRUE,
#'   colorscales = list("Blues", "Reds"),
#'   opacities   = c(0.25, 0.35),
#'   compute_volume = TRUE
#' )
#' out$volume$estimate       # approximately 4/3 * pi * (1.2^3 - 0.8^3)
#'
#' # Example 2: Spherical cap: 0 <= r <= 1, phi in [0, pi/3]
#' R1 <- function(th, ph) 0
#' R2 <- function(th, ph) 1
#' out2 <- solid_spherical3d(
#'   R1, R2,
#'   theta_range = c(0, 2*pi),
#'   phi_range   = c(0, pi/3),
#'   plot = TRUE,
#'   compute_volume = TRUE
#' )
#' out2$volume$estimate      # analytic value matches the standard spherical cap formula
#' \dontshow{\}}
#'
#' @export
solid_spherical3d <- function(
    R1, R2,
    theta_range = c(0, 2*pi),
    phi_range   = c(0, pi),
    n_theta = 160, n_phi = 120,
    plot = TRUE,
    show_surfaces = c(TRUE, TRUE),                 # show R1 (inner), R2 (outer)
    colorscales = list("Blues", "Reds"),
    opacities   = c(0.30, 0.35),
    show_surface_grid  = TRUE,
    surface_grid_color = "rgba(60,80,200,0.25)",
    surface_grid_width = 1,
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z")
    ),
    bg = list(paper = "white", plot = "white"),
    compute_volume = FALSE,
    vol_method = c("adaptive", "grid"),
    n_th_vol = 600, n_ph_vol = 400
) {
  # ---- validations
  if (!is.function(R1) || !is.function(R2))
    stop("'R1' and 'R2' must be functions (theta, phi) -> scalar radius.", call. = FALSE)
  if (!is.numeric(theta_range) || length(theta_range) != 2L || any(!is.finite(theta_range)))
    stop("'theta_range' must be c(min, max) numeric.", call. = FALSE)
  if (!is.numeric(phi_range) || length(phi_range) != 2L || any(!is.finite(phi_range)))
    stop("'phi_range' must be c(min, max) numeric.", call. = FALSE)
  if (theta_range[2] <= theta_range[1]) stop("theta_range must have max > min.", call. = FALSE)
  if (phi_range[2]   <= phi_range[1])   stop("phi_range must have max > min.",   call. = FALSE)

  vol_method <- match.arg(vol_method)

  # Coerce mesh sizes to integers and validate
  n_theta  <- as.integer(n_theta)
  n_phi    <- as.integer(n_phi)
  n_th_vol <- as.integer(n_th_vol)
  n_ph_vol <- as.integer(n_ph_vol)

  if (!is.finite(n_theta) || n_theta < 2L ||
      !is.finite(n_phi)   || n_phi   < 2L) {
    stop("'n_theta' and 'n_phi' must be integer values >= 2.", call. = FALSE)
  }

  if (vol_method == "grid" &&
      (!is.finite(n_th_vol) || n_th_vol < 2L ||
       !is.finite(n_ph_vol) || n_ph_vol < 2L)) {
    stop("'n_th_vol' and 'n_ph_vol' must be integer values >= 2 for grid volume.",
         call. = FALSE)
  }

  # ---- helpers: colorscale handling (flat color or plotly scale)
  is_color <- function(x) {
    if (!is.character(x) || length(x) != 1) return(FALSE)
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
    sprintf("rgba(%g,%g,%g,%g)", 255*rgb[1], 255*rgb[2], 255*rgb[3], a_)
  }
  as_colorscale <- function(x, alpha = NULL) {
    if (is.null(x)) return(NULL)
    # x can be plotly named scale, single color, or vector of colors
    if (is.list(x) && length(x) >= 2 && is.numeric(x[[1]][[1]])) return(x)
    if (is.character(x) && length(x) == 1) {
      if (is_color(x)) {
        ccol <- as_rgba(x, alpha)
        return(list(list(0, ccol), list(1, ccol)))
      } else {
        return(x)
      }
    }
    if (is.character(x) && length(x) > 1) {
      cols <- vapply(x, as_rgba, character(1), alpha = alpha)
      pos  <- seq(0, 1, length.out = length(cols))
      return(lapply(seq_along(cols), function(i) list(pos[i], cols[i])))
    }
    stop("Unrecognized colorscale format.", call. = FALSE)
  }
  to_two <- function(x) {
    if (length(x) == 1L) rep(x, 2L) else x
  }

  # Normalize colorscales input: single scale or list of length 2
  if (length(colorscales) == 1L && !is.list(colorscales)) {
    colorscales <- list(colorscales, colorscales)
  }
  if (length(colorscales) != 2L) {
    stop("'colorscales' must be a single scale or a list (length 2).",
         call. = FALSE)
  }

  cs1 <- as_colorscale(colorscales[[1]])
  cs2 <- as_colorscale(colorscales[[2]])
  opacities <- to_two(opacities)

  # ---- parameter grids
  th_seq <- seq(theta_range[1], theta_range[2], length.out = n_theta)
  ph_seq <- seq(phi_range[1],   phi_range[2],   length.out = n_phi)

  # ---- surface builder (R(theta, phi) -> (x,y,z))
  build_surface <- function(Rfun) {
    # matrices: phi varies by row, theta by column
    TH <- matrix(rep(th_seq, each = n_phi), nrow = n_phi)  # n_phi x n_theta
    PH <- matrix(rep(ph_seq, times = n_theta), nrow = n_phi)
    # evaluate R safely (elementwise)
    Rvals <- matrix(
      mapply(function(th, ph) as.numeric(Rfun(th, ph)), as.numeric(TH), as.numeric(PH)),
      nrow = n_phi, ncol = n_theta
    )
    # spherical -> Cartesian
    X <- Rvals * sin(PH) * cos(TH)
    Y <- Rvals * sin(PH) * sin(TH)
    Z <- Rvals * cos(PH)
    list(X = X, Y = Y, Z = Z)
  }

  R1_surf <- if (isTRUE(show_surfaces[1])) build_surface(R1) else NULL
  R2_surf <- if (isTRUE(show_surfaces[2])) build_surface(R2) else NULL

  # ---- plot
  fig <- NULL
  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("Plotly is required for plotting.", call. = FALSE)
    } else {
      contours_arg <- if (isTRUE(show_surface_grid)) list(
        x = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
        y = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
        z = list(show = FALSE)
      ) else NULL

      plt <- plotly::plot_ly()
      if (!is.null(R1_surf)) {
        plt <- plt |>
          plotly::add_surface(
            x = R1_surf$X, y = R1_surf$Y, z = R1_surf$Z,
            colorscale = cs1, showscale = FALSE,
            opacity = opacities[1], contours = contours_arg
          )
      }
      if (!is.null(R2_surf)) {
        plt <- plt |>
          plotly::add_surface(
            x = R2_surf$X, y = R2_surf$Y, z = R2_surf$Z,
            colorscale = cs2, showscale = FALSE,
            opacity = opacities[2], contours = contours_arg
          )
      }
      plt <- plt |>
        plotly::layout(
          title = "Spherical solid",
          scene = scene, paper_bgcolor = bg$paper, plot_bgcolor = bg$plot
        )
      fig <- plt
      print(fig)
    }
  }

  # ---- volume
  volume <- NULL
  if (isTRUE(compute_volume)) {
    th_min <- theta_range[1]; th_max <- theta_range[2]
    ph_min <- phi_range[1];   ph_max <- phi_range[2]

    # scalar-safe inner integrand
    shell_term_scalar <- function(th, ph) {
      r1 <- as.numeric(R1(th, ph)); r2 <- as.numeric(R2(th, ph))
      # guard swapped/negative radii
      if (!is.finite(r1) || !is.finite(r2)) return(0)
      if (r2 < r1) {
        tmp <- r1
        r1  <- r2
        r2  <- tmp
      }
      if (r2 <= 0) return(0)
      r1 <- max(0, r1)  # clamp inner radius at 0
      ((r2^3 - r1^3) / 3) * sin(ph)
    }

    if (vol_method == "adaptive") {
      # integrate over phi
      inner_phi <- function(th) {
        fphi <- function(ph) vapply(ph, function(p_) shell_term_scalar(th, p_), numeric(1))
        stats::integrate(function(ph) fphi(ph), lower = ph_min, upper = ph_max, rel.tol = 1e-6)$value
      }
      # integrate over theta
      fth <- function(th) vapply(th, function(t_) inner_phi(t_), numeric(1))
      V <- stats::integrate(function(th) fth(th), lower = th_min, upper = th_max, rel.tol = 1e-6)$value
      volume <- list(estimate = V, method = "adaptive")

    } else {
      # grid trapezoid rule in (theta, phi)
      th_g <- seq(th_min, th_max, length.out = n_th_vol)
      ph_g <- seq(ph_min, ph_max, length.out = n_ph_vol)
      dth  <- (th_max - th_min) / (n_th_vol - 1)
      dph  <- (ph_max - ph_min) / (n_ph_vol - 1)

      # compute shell term on grid
      TERM <- matrix(NA_real_, nrow = n_ph_vol, ncol = n_th_vol)
      for (j in seq_along(th_g)) {
        thj <- th_g[j]
        TERM[, j] <- vapply(ph_g, function(ph_) shell_term_scalar(thj, ph_), numeric(1))
      }
      # trapezoid in phi for each theta
      Iphi <- dph * (colSums(TERM) - 0.5 * TERM[1, ] - 0.5 * TERM[n_ph_vol, ])
      # trapezoid in theta
      V <- dth * (sum(Iphi) - 0.5 * Iphi[1] - 0.5 * Iphi[length(Iphi)])
      volume <- list(estimate = V, method = "grid", n_th = n_th_vol, n_ph = n_ph_vol)
    }
  }

  list(
    theta_seq = th_seq,
    phi_seq   = ph_seq,
    R1_surf   = R1_surf,
    R2_surf   = R2_surf,
    fig       = fig,
    volume    = volume
  )
}
