#' Plot a 3D parametric curve with plotly
#'
#' Creates an interactive 3D plot of a parametric curve given
#' a tibble with columns \code{t}, \code{x}, \code{y}, \code{z}
#' (typically produced by [curve_sample3d()]).
#'
#' This function requires the \pkg{plotly} package to be installed.
#'
#' @param data Tibble with columns \code{t}, \code{x}, \code{y}, \code{z}.
#' @param mode Character string. Plotly trace mode, e.g. \code{"lines"} or
#'   \code{"lines+markers"}.
#' @param line List with line styling options, such as
#'   \code{list(color = "blue", width = 3, dash = "solid")}.
#' @param marker Optional list with marker styling options, or \code{NULL}
#'   to omit markers.
#' @param title Optional character string for the plot title.
#' @param scene List specifying 3D axis titles and options, passed to
#'   \code{plotly::layout()}, typically including \code{xaxis}, \code{yaxis},
#'   and \code{zaxis}.
#' @param bg List with background colors, e.g.
#'   \code{list(paper = "white", plot = "white")}.
#'
#' @return
#' A \pkg{plotly} object, which is printed for interactive visualization.
#'
#' @seealso [curve_sample3d()], [arc_length3d()]
#'
#' @examples
#' data <- curve_sample3d(
#'   function(t) 2 * cos(t),
#'   function(t) 3 * sin(t),
#'   function(t) t / 5,
#'   0, 2 * pi, 100
#' )
#' \donttest{
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   plot_curve3d(data, line = list(color = "red", width = 4))
#' }
#' }
#'
#' @export
plot_curve3d <- function(
    data,
    mode  = "lines",
    line  = list(color = "blue", width = 3, dash = "solid"),
    marker = NULL,
    title  = NULL,
    scene  = list(
      xaxis = list(title = "x(t)"),
      yaxis = list(title = "y(t)"),
      zaxis = list(title = "z(t)")
    ),
    bg = list(paper = "white", plot = "white")
) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("You need to install 'plotly' to use plot_curve3d().")
  }

  plt <- plotly::plot_ly(
    data = data, x = ~x, y = ~y, z = ~z,
    type = "scatter3d", mode = mode, line = line
  )

  if (!is.null(marker)) {
    plt <- plt |>
      plotly::style(marker = marker)
  }

  plotly::layout(
    plt,
    title = title,
    scene = scene,
    paper_bgcolor = bg$paper,
    plot_bgcolor  = bg$plot
  )
}
