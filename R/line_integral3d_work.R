#' 3D line integral with work visualization
#'
#' Compute a numerical approximation of the line integral of a vector field
#' along a parametric space curve and optionally draw a three dimensional
#' visualization of the curve together with arrows of the field.
#'
#' The parameter t runs from a to b. The curve r(t) must return a numeric
#' vector of length three. The field field(x, y, z) may optionally also depend
#' on t through a fourth argument.
#'
#' @param field Function that represents the vector field. It must be a function
#'   of three or four numeric arguments. In the three argument form the
#'   arguments are x, y and z. In the four argument form the arguments are
#'   x, y, z and t. The function must return a numeric vector of length three.
#' @param r Function of one numeric argument t that returns a numeric vector
#'   c(x, y, z) of length three with the coordinates of the curve.
#' @param a Numeric scalar with the lower value of the parameter interval.
#' @param b Numeric scalar with the upper value of the parameter interval.
#'   It must satisfy \code{b > a}.
#' @param plot Logical value. If \code{TRUE}, a plotly object with the field
#'   arrows and the curve is created.
#' @param n_curve Integer number of sampled points on the curve.
#' @param n_field Integer number of grid points per axis for the field
#'   arrows. The total number of arrows is \code{n_field^3}.
#' @param field_ranges Optional list with named components \code{x}, \code{y}
#'   and \code{z}, each a numeric vector of length two giving the range used
#'   to build the field grid. If \code{NULL}, the ranges are taken from the
#'   bounding box of the curve and expanded by \code{pad}.
#' @param pad Numeric fraction used to expand the automatic field ranges.
#' @param arrows Character string that selects which arrow style is drawn.
#'   One of \code{"both"}, \code{"line"}, \code{"cone"} or \code{"none"}.
#' @param arrow_scale Numeric factor that controls the length of the arrows
#'   as a fraction of the size of the domain.
#' @param normalize_bias Numeric value used to regularize the scaling of
#'   the field magnitude when computing arrows.
#' @param scene List of plotly scene options passed to \code{plotly::layout()}.
#' @param bg List with two character elements named \code{paper} and
#'   \code{plot} that control the background colours in the plotly layout.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{value}: numeric value of the line integral (work).
#'   \item \code{fig}: plotly object when \code{plot = TRUE}, otherwise
#'     \code{NULL}.
#'   \item \code{field_box}: list with numeric ranges used for the field
#'     grid, with components \code{x}, \code{y} and \code{z}.
#' }
#'
#' @examples
#' \donttest{
#' field <- function(x, y, z) c(-y, x, 0.2*z)
#' r <- function(t) c(cos(t), sin(t), 0.25*t)
#' out <- line_integral3d_work(
#'   field = field, r = r, a = 0, b = 2*pi,
#'   plot = FALSE, n_curve = 200, n_field = 5
#' )
#' out$value
#' }
#'
#' @export
line_integral3d_work <- function(
    field, r, a, b,
    plot = TRUE,
    n_curve = 600,
    n_field = 7,
    field_ranges = NULL,
    pad = 0.15,
    arrows = c("both", "line", "cone", "none"),
    arrow_scale = 0.10,
    normalize_bias = 1,
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z")
    ),
    bg = list(paper = "white", plot = "white")
) {
  arrows <- match.arg(arrows)

  if (!is.function(field)) {
    stop("'field' must be function(x,y,z) or function(x,y,z,t).", call. = FALSE)
  }
  if (!is.function(r)) {
    stop("'r' must be function(t) -> c(x,y,z).", call. = FALSE)
  }
  if (!is.numeric(a) || !is.numeric(b) || length(a) != 1L || length(b) != 1L ||
      !is.finite(a) || !is.finite(b) || b <= a) {
    stop("'a' and 'b' must be finite scalars with b > a.", call. = FALSE)
  }
  if (!is.numeric(n_curve) || length(n_curve) != 1L || !is.finite(n_curve) || n_curve < 10) {
    stop("'n_curve' must be an integer >= 10.", call. = FALSE)
  }
  if (!is.numeric(n_field) || length(n_field) != 1L || !is.finite(n_field) || n_field < 2) {
    stop("'n_field' must be an integer >= 2.", call. = FALSE)
  }

  # helper to call field with or without time
  field_eval <- function(x, y, z, t) {
    if (length(formals(field)) >= 4L) field(x, y, z, t) else field(x, y, z)
  }

  # curve sampling
  t_curve <- seq(a, b, length.out = n_curve)
  Rxyz <- t(vapply(t_curve, r, numeric(3)))
  cx <- Rxyz[, 1]
  cy <- Rxyz[, 2]
  cz <- Rxyz[, 3]

  # numerical r'(t) (central differences along sampled grid)
  dt <- (b - a) / (n_curve - 1L)
  Rp <- rbind(
    (Rxyz[2, ] - Rxyz[1, ]) / dt,
    (Rxyz[3:n_curve, ] - Rxyz[1:(n_curve - 2L), ]) / (2 * dt),
    (Rxyz[n_curve, ] - Rxyz[n_curve - 1L, ]) / dt
  )

  # integrand field(r(t)) . r'(t)
  Fr <- t(mapply(function(x, y, z, t) field_eval(x, y, z, t), cx, cy, cz, t_curve))
  dot <- rowSums(Fr * Rp)

  # composite Simpson rule (ensure even number of subintervals)
  n_seg <- n_curve - 1L
  if (n_seg %% 2L == 1L) {
    t_curve <- t_curve[-length(t_curve)]
    dot <- dot[-length(dot)]
    n_seg <- n_seg - 1L
  }
  hS <- (t_curve[length(t_curve)] - t_curve[1L]) / n_seg
  work <- hS * (
    dot[1L] + dot[length(dot)] +
      4 * sum(dot[seq(2L, length(dot) - 1L, by = 2L)]) +
      2 * sum(dot[seq(3L, length(dot) - 2L, by = 2L)])
  ) / 3

  # arrow domain (auto from curve bbox unless provided)
  if (is.null(field_ranges)) {
    rx <- range(cx)
    ry <- range(cy)
    rz <- range(cz)
    pad_range <- function(rg) {
      span <- diff(rg)
      if (!is.finite(span) || span <= 0) span <- 1
      c(rg[1] - pad * span, rg[2] + pad * span)
    }
    xlim <- pad_range(rx)
    ylim <- pad_range(ry)
    zlim <- pad_range(rz)
  } else {
    if (!is.list(field_ranges) || is.null(field_ranges$x) ||
        is.null(field_ranges$y) || is.null(field_ranges$z)) {
      stop("'field_ranges' must be a list with components x, y, z.", call. = FALSE)
    }
    xlim <- field_ranges$x
    ylim <- field_ranges$y
    zlim <- field_ranges$z
  }

  # grid for arrows
  xs <- seq(xlim[1], xlim[2], length.out = n_field)
  ys <- seq(ylim[1], ylim[2], length.out = n_field)
  zs <- seq(zlim[1], zlim[2], length.out = n_field)
  G <- as.matrix(expand.grid(xs, ys, zs))
  colnames(G) <- c("x", "y", "z")

  V <- t(mapply(function(x, y, z) field_eval(x, y, z, 0), G[, "x"], G[, "y"], G[, "z"]))
  vmag <- sqrt(rowSums(V * V))
  Vsat <- V / sqrt(vmag^2 + normalize_bias)
  lsat <- sqrt(rowSums(Vsat * Vsat))
  nz <- lsat > 0
  diru <- Vsat
  diru[nz, ] <- Vsat[nz, , drop = FALSE] / lsat[nz]
  diru[!nz, ] <- 0
  span <- max(diff(xlim), diff(ylim), diff(zlim))
  if (!is.finite(span) || span <= 0) span <- 1
  L <- arrow_scale * span * lsat
  tips <- G + diru * L

  fig <- NULL
  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("Install 'plotly' to plot.", call. = FALSE)
    } else {
      plt <- plotly::plot_ly()

      # field body (lines)
      if (arrows %in% c("both", "line", "cone")) {
        xsL <- as.numeric(rbind(G[, "x"], tips[, "x"], NA))
        ysL <- as.numeric(rbind(G[, "y"], tips[, "y"], NA))
        zsL <- as.numeric(rbind(G[, "z"], tips[, "z"], NA))
        plt <- plotly::add_trace(
          plt, x = xsL, y = ysL, z = zsL,
          type = "scatter3d", mode = "lines",
          line = list(color = "rgba(60,60,60,0.55)", width = 2),
          showlegend = FALSE, hoverinfo = "none"
        )
      }

      # field heads (cones if available)
      if (arrows %in% c("both", "cone")) {
        if ("add_cone" %in% getNamespaceExports("plotly")) {
          plt <- plotly::add_trace(
            plt,
            type = "cone",
            x = tips[, "x"], y = tips[, "y"], z = tips[, "z"],
            u = diru[, 1], v = diru[, 2], w = diru[, 3],
            anchor = "tip",
            sizemode = "absolute",
            sizeref = max(1e-8, arrow_scale * span * 0.35),
            colorscale = list(list(0, "#444"), list(1, "#444")),
            showscale = FALSE, opacity = 0.9
          )
        }
      }

      # curve coloured by instantaneous work density field . r'
      plt <- plotly::add_trace(
        plt, x = cx, y = cy, z = cz,
        type = "scatter3d", mode = "lines+markers",
        line = list(width = 6, color = dot, colorscale = "RdBu", reversescale = TRUE),
        marker = list(size = 2, color = dot, colorscale = "RdBu", reversescale = TRUE),
        showlegend = FALSE
      )

      # start marker
      plt <- plotly::add_markers(
        plt,
        x = cx[1L], y = cy[1L], z = cz[1L],
        marker = list(size = 4, color = "forestgreen"),
        showlegend = FALSE, hoverinfo = "none"
      )
      # end marker
      plt <- plotly::add_markers(
        plt,
        x = cx[length(cx)], y = cy[length(cy)], z = cz[length(cz)],
        marker = list(size = 4, color = "firebrick"),
        showlegend = FALSE, hoverinfo = "none"
      )

      ttl <- sprintf("3D line integral (work) = %.6g", work)
      fig <- plotly::layout(
        plt, title = ttl, scene = scene,
        paper_bgcolor = bg$paper, plot_bgcolor = bg$plot
      )
      print(fig)
    }
  }

  list(
    value = work,
    fig = fig,
    field_box = list(x = xlim, y = ylim, z = zlim)
  )
}

