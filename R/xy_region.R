#' Planar region between two curves y = H1(x) and y = H2(x)
#'
#' @description
#' Constructs a numerical representation of the planar region bounded by
#' two functions of one variable. The region consists of all points whose
#' horizontal coordinate lies between \code{a} and \code{b}, and whose
#' vertical coordinate lies between the values returned by \code{H1(x)}
#' and \code{H2(x)}. Optionally, the region can be displayed either as a
#' two-dimensional filled subset of the plane or as a thin surface in
#' three dimensions using \pkg{plotly}.
#'
#' @details
#' The function samples the interval \code{[a, b]} at \code{n_curve} points
#' to represent the boundary curves. The interval \code{[a, b]} is also
#' subdivided into \code{D} vertical strips. For each strip, the values
#' \code{H1(x)} and \code{H2(x)} are evaluated to define the vertical bounds
#' of the region.
#'
#' Depending on the arguments, the function can:
#' \itemize{
#'   \item build a data grid suitable for numerical integration or
#'         visualization,
#'   \item draw a two-dimensional depiction of the region, possibly filled
#'         with a selected color,
#'   \item generate a simple three-dimensional visualization where the
#'         region is drawn as a thin plate at a chosen height.
#' }
#'
#' Additional options allow drawing grid lines, showing the boundary curves,
#' controlling colors and transparency, and adjusting the aspect ratio.
#'
#' @param H1,H2 Functions of one variable returning the lower and upper
#'   vertical bounds at each value of \code{x}.
#' @param a,b Numeric endpoints of the interval for \code{x}. It is assumed
#'   that \code{a < b}.
#' @param D Integer giving the number of subdivisions in the horizontal
#'   direction (number of vertical strips).
#' @param plot Logical; if \code{TRUE}, produces a visualization using
#'   \pkg{plotly}.
#' @param n_curve Integer giving the number of points for sampling the
#'   boundary curves.
#' @param fill Logical; if \code{TRUE}, fills the two-dimensional region.
#' @param fillcolor Character string defining the fill color in 2D mode.
#' @param boundary_line List with \pkg{plotly} style options for drawing the
#'   boundary curves.
#' @param partition_line List with style parameters for vertical partition
#'   lines.
#' @param show_end_edges Logical; if \code{TRUE}, draws boundary segments at
#'   the endpoints \code{x = a} and \code{x = b}.
#' @param axis_equal Logical; if \code{TRUE}, enforces equal scaling on both
#'   axes in two dimensions.
#' @param as_3d Logical; if \code{TRUE}, draws the region as a thin
#'   three-dimensional plate.
#' @param plane_z Numeric height at which to draw the region when
#'   \code{as_3d = TRUE}.
#' @param n_u Integer number of internal subdivisions used for
#'   discretization of the region in the cross-section (between lower and
#'   upper boundary).
#' @param surface_colorscale Character string specifying a \pkg{plotly}
#'   colorscale for the three-dimensional mode.
#' @param surface_opacity Numeric value between 0 and 1 controlling the
#'   transparency of the surface in 3D mode.
#' @param show_surface_grid Logical; if \code{TRUE}, overlays grid lines on
#'   the plotted surface in 3D mode.
#' @param surface_grid_color Character string giving the color of the grid
#'   lines in 3D mode.
#' @param surface_grid_width Numeric width of the grid lines in 3D mode.
#' @param scene Optional list of \pkg{plotly} scene parameters for
#'   three-dimensional rendering.
#' @param bg Optional list defining the background colors of the figure,
#'   typically with components \code{paper} and \code{plot}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{x}: the sample points along the horizontal axis,
#'   \item \code{y1}, \code{y2}: the sampled boundary values,
#'   \item \code{y_low}, \code{y_high}: the lower and upper envelopes
#'         \code{pmin(H1,H2)} and \code{pmax(H1,H2)},
#'   \item \code{x_part}: the partition points used in the horizontal
#'         direction,
#'   \item \code{fig}: a \pkg{plotly} object for visualization if
#'         \code{plot = TRUE} and \pkg{plotly} is available; otherwise
#'         \code{NULL}.
#' }
#'
#' @examples
#' H1 <- function(x) 0
#' H2 <- function(x) 1 - x
#' xy_region(H1, H2, a = 0, b = 1, D = 20, plot = FALSE)
#'
#' @export
xy_region <- function(
    H1, H2,
    a, b,
    D,
    plot = TRUE,
    n_curve = 800,
    fill = FALSE,
    fillcolor = "rgba(49,130,189,0.25)",
    boundary_line  = list(color = "blue", width = 2),
    partition_line = list(color = "blue", width = 1),
    show_end_edges = TRUE,
    axis_equal = TRUE,
    as_3d = FALSE,
    plane_z = 0,
    n_u = 30,
    surface_colorscale = "Blues",
    surface_opacity    = 0.30,
    show_surface_grid  = TRUE,
    surface_grid_color = "rgba(60,80,200,0.25)",
    surface_grid_width = 1,
    scene = list(
      aspectmode = "data",
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z")
    ),
    bg = list(paper = "white", plot = "white")
) {
  if (!is.function(H1) || !is.function(H2)) {
    stop("'H1' and 'H2' must be functions of one variable.", call. = FALSE)
  }
  if (!is.numeric(a) || !is.numeric(b) || length(a) != 1L || length(b) != 1L || b <= a) {
    stop("'a' and 'b' must be numeric scalars with b > a.", call. = FALSE)
  }

  D <- as.integer(D)
  if (!is.finite(D) || D < 1L) {
    stop("'D' must be an integer >= 1 (number of subdivisions).", call. = FALSE)
  }

  n_curve <- as.integer(n_curve)
  if (!is.finite(n_curve) || n_curve < 2L) {
    stop("'n_curve' must be an integer >= 2.", call. = FALSE)
  }

  n_u <- as.integer(n_u)
  if (!is.finite(n_u) || n_u < 2L) {
    stop("'n_u' must be an integer >= 2.", call. = FALSE)
  }

  # Number of vertical strips
  n <- D

  # Curve mesh
  xs  <- seq(a, b, length.out = n_curve)
  y1  <- H1(xs)
  y2  <- H2(xs)
  ylo <- pmin(y1, y2)
  yhi <- pmax(y1, y2)

  # Partition positions (includes a and b)
  x_part <- seq(a, b, length.out = n + 1L)

  fig <- NULL
  if (isTRUE(plot)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("To plot you need 'plotly' installed.", call. = FALSE)
    } else {
      if (!isTRUE(as_3d)) {
        # ---------- 2D MODE ----------
        fig <- plotly::plot_ly()

        if (isTRUE(fill)) {
          fig <- fig |>
            plotly::add_trace(
              x = c(xs, rev(xs)),
              y = c(ylo, rev(yhi)),
              type = "scatter",
              mode = "lines",
              fill = "toself",
              fillcolor = fillcolor,
              line = list(color = "rgba(0,0,0,0)"),
              hoverinfo = "none",
              showlegend = FALSE
            )
        }

        fig <- fig |>
          plotly::add_trace(
            x = xs, y = y1,
            type = "scatter", mode = "lines",
            line = boundary_line, hoverinfo = "none", showlegend = FALSE
          ) |>
          plotly::add_trace(
            x = xs, y = y2,
            type = "scatter", mode = "lines",
            line = boundary_line, hoverinfo = "none", showlegend = FALSE
          )

        k_seq <- if (isTRUE(show_end_edges)) 0:n else 1:(n - 1L)
        if (length(k_seq) > 0L) {
          for (k in k_seq) {
            xk <- a + (b - a) * k / n
            yk1 <- H1(xk); yk2 <- H2(xk)
            yklo <- min(yk1, yk2); ykhi <- max(yk1, yk2)
            fig <- fig |>
              plotly::add_trace(
                x = c(xk, xk),
                y = c(yklo, ykhi),
                type = "scatter", mode = "lines",
                line = partition_line, hoverinfo = "none", showlegend = FALSE
              )
          }
        }

        layout_args <- list(
          xaxis = list(title = "x"),
          yaxis = list(title = "y")
        )
        if (isTRUE(axis_equal)) {
          layout_args$yaxis$scaleanchor <- "x"
          layout_args$yaxis$scaleratio  <- 1
        }

        fig <- fig |>
          plotly::layout(
            title = "Region in the xy-plane",
            xaxis = layout_args$xaxis,
            yaxis = layout_args$yaxis
          )

      } else {
        # ---------- 3D MODE ----------
        fig <- plotly::plot_ly()

        # Fill as a surface S(u,j) with u in [0,1]: y(u) = (1-u)*ylo + u*yhi
        if (isTRUE(fill)) {
          U <- seq(0, 1, length.out = n_u)
          Xmat <- matrix(rep(xs, each = n_u), nrow = n_u)
          Ymat <- outer(U, seq_along(xs), function(u, j) (1 - u) * ylo[j] + u * yhi[j])
          Zmat <- matrix(plane_z, nrow = n_u, ncol = length(xs))

          contours_arg <- if (isTRUE(show_surface_grid)) list(
            x = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
            y = list(show = TRUE, color = surface_grid_color, width = surface_grid_width),
            z = list(show = FALSE)
          ) else NULL

          fig <- fig |>
            plotly::add_surface(
              x = Xmat, y = Ymat, z = Zmat,
              colorscale = surface_colorscale,
              showscale  = FALSE,
              opacity    = surface_opacity,
              contours   = contours_arg
            )
        }

        # Boundary curves over z = plane_z
        fig <- fig |>
          plotly::add_trace(
            x = xs, y = y1, z = rep(plane_z, length(xs)),
            type = "scatter3d", mode = "lines",
            line = boundary_line, hoverinfo = "none", showlegend = FALSE
          ) |>
          plotly::add_trace(
            x = xs, y = y2, z = rep(plane_z, length(xs)),
            type = "scatter3d", mode = "lines",
            line = boundary_line, hoverinfo = "none", showlegend = FALSE
          )

        # Vertical partitions over z = plane_z
        k_seq <- if (isTRUE(show_end_edges)) 0:n else 1:(n - 1L)
        if (length(k_seq) > 0L) {
          for (k in k_seq) {
            xk <- a + (b - a) * k / n
            yk1 <- H1(xk); yk2 <- H2(xk)
            yklo <- min(yk1, yk2); ykhi <- max(yk1, yk2)
            fig <- fig |>
              plotly::add_trace(
                x = c(xk, xk),
                y = c(yklo, ykhi),
                z = c(plane_z, plane_z),
                type = "scatter3d", mode = "lines",
                line = partition_line, hoverinfo = "none", showlegend = FALSE
              )
          }
        }

        fig <- fig |>
          plotly::layout(
            title = sprintf("Region on plane z = %g", plane_z),
            scene = scene,
            paper_bgcolor = bg$paper,
            plot_bgcolor  = bg$plot
          )
      }

      print(fig)
    }
  }

  list(
    x      = xs,
    y1     = y1,
    y2     = y2,
    y_low  = ylo,
    y_high = yhi,
    x_part = x_part,
    fig    = fig
  )
}
