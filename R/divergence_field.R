#' Numerical divergence of a vector field
#'
#' Computes the divergence of a vector field at a given point using
#' central finite differences. The vector field \code{field} must take a
#' numeric vector \code{x} and return a numeric vector of the same length.
#'
#' Optionally, if the dimension is 2 or 3 and \code{plot = TRUE},
#' a simple visualization is produced using \pkg{plotly}.
#'
#' @param field Function of the form \code{field(x)} returning a numeric vector
#'   of the same length as \code{x}.
#' @param x0 Numeric vector giving the evaluation point.
#' @param h Step size for finite differences. Can be:
#'   \itemize{
#'     \item \code{NULL}: an automatic step is selected as
#'           \code{1e-4 * (1 + abs(x0))};
#'     \item A numeric scalar: same step for all components;
#'     \item A numeric vector of the same length as \code{x0}.
#'   }
#' @param plot Logical; if \code{TRUE} and the dimension is 2 or 3,
#'   a basic visualization is drawn with \pkg{plotly}.
#'
#' @return A numeric scalar: the divergence evaluated at \code{x0}.
#'
#' @examples
#' field <- function(x) c(x[1] + x[2], x[2] - x[1])
#' divergence_field(field, c(0.5, -0.2))
#'
#' @export
divergence_field <- function(field, x0, h = NULL, plot = FALSE) {
  if (!is.function(field)) {
    stop("'field' must be a function(x).", call. = FALSE)
  }
  if (!is.numeric(x0)) {
    stop("'x0' must be numeric.", call. = FALSE)
  }
  x0 <- as.numeric(x0)
  n  <- length(x0)

  Fx0 <- field(x0)
  if (!is.numeric(Fx0) || length(Fx0) != n) {
    stop("'field(x)' must return a numeric vector of length n.", call. = FALSE)
  }

  if (is.null(h)) {
    h <- 1e-4 * (1 + abs(x0))
  } else if (length(h) == 1L) {
    h <- rep(as.numeric(h), n)
  } else if (length(h) != n) {
    stop("'h' must be scalar, a vector of length n, or NULL.", call. = FALSE)
  }

  if (any(!is.finite(h)) || any(h <= 0)) {
    stop("'h' must contain positive finite values.", call. = FALSE)
  }

  div <- 0
  for (i in seq_len(n)) {
    ei <- rep(0, n)
    ei[i] <- h[i]
    Fp <- field(x0 + ei)
    Fm <- field(x0 - ei)
    if (!is.numeric(Fp) || length(Fp) != n ||
        !is.numeric(Fm) || length(Fm) != n) {
      stop("field(x) must return a numeric vector of length n for all x.", call. = FALSE)
    }
    div <- div + (Fp[i] - Fm[i]) / (2 * h[i])
  }

  # Optional visualization for 2D / 3D
  if (isTRUE(plot) && (n %in% c(2L, 3L))) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("Package 'plotly' is required for visualization.", call. = FALSE)
    } else {
      if (n == 2L) {
        # small grid around x0
        xs <- seq(x0[1] - 1, x0[1] + 1, length.out = 5L)
        ys <- seq(x0[2] - 1, x0[2] + 1, length.out = 5L)
        grid <- expand.grid(x = xs, y = ys)
        vecs <- t(apply(grid, 1L, field))

        fig <- plotly::plot_ly()
        fig <- fig |>
          plotly::add_trace(
            x = grid$x, y = grid$y,
            u = vecs[, 1], v = vecs[, 2],
            type = "cone",
            sizemode = "absolute",
            sizeref = 0.5,
            colorscale = if (div >= 0) "Reds" else "Blues",
            showscale = FALSE
          ) |>
          plotly::layout(
            title = sprintf(
              "Divergence in 2D at x0 = (%.2f, %.2f): %.4f",
              x0[1], x0[2], div
            ),
            xaxis = list(title = "x"),
            yaxis = list(title = "y")
          )
        print(fig)
      } else if (n == 3L) {
        # simple starburst around x0
        directions <- rbind(
          c(1, 0, 0), c(-1, 0, 0),
          c(0, 1, 0), c(0, -1, 0),
          c(0, 0, 1), c(0, 0, -1)
        )
        arrows <- sweep(directions, 2L, x0, "+")

        fig <- plotly::plot_ly()
        for (i in seq_len(nrow(directions))) {
          fig <- fig |>
            plotly::add_trace(
              x = c(x0[1], arrows[i, 1]),
              y = c(x0[2], arrows[i, 2]),
              z = c(x0[3], arrows[i, 3]),
              type = "scatter3d",
              mode = "lines",
              line = list(
                color = if (div >= 0) "red" else "blue",
                width = 6
              ),
              showlegend = FALSE,
              hoverinfo = "none"
            )
        }
        fig <- fig |>
          plotly::layout(
            title = sprintf(
              "Divergence in 3D at x0 = (%.2f, %.2f, %.2f): %.4f",
              x0[1], x0[2], x0[3], div
            ),
            scene = list(
              xaxis = list(title = "x"),
              yaxis = list(title = "y"),
              zaxis = list(title = "z")
            )
          )
        print(fig)
      }
    }
  }

  div
}
