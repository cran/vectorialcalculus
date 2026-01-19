#' @title Unified Numerical Double Integration
#'
#' @description
#' Calculates the definite double integral of a function f(x, y) over a
#' general region D, which can be defined as Type I (dy dx) or Type II (dx dy).
#' Uses the Composite Simpson's Rule for numerical approximation.
#'
#' @param f A function in R of two variables, f(x, y), returning a numeric value.
#' @param const_min The constant lower limit of the outer integration (a for Type I, c for Type II).
#' @param const_max The constant upper limit of the outer integration (b for Type I, d for Type II).
#' @param limit1 A function in R of one variable defining the inner integral's lower limit (h1(x) or h1(y)).
#' @param limit2 A function in R of one variable defining the inner integral's upper limit (h2(x) or h2(y)).
#' @param region_type A string specifying the region type: "type1" (dy dx) or "type2" (dx dy). Default is "type1".
#' @param n_outer Number of subintervals for the outer integration. Must be even. Default is 100.
#' @param n_inner Number of subintervals for the inner integration. Must be even. Default is 100.
#' @param plot_domain Logical. If TRUE, generates a ggplot2 plot of the integration domain. Default is TRUE.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{integral_value}: The calculated approximation of the integral.
#'     \item \code{domain_plot}: The ggplot2 object representing the domain (if plot_domain = TRUE).
#'   }
#'
#' @export
integrate_double_xy <- function(
    f,
    const_min, const_max,
    limit1, limit2,
    region_type  = "type1",
    n_outer      = 100,
    n_inner      = 100,
    plot_domain  = TRUE
) {
  # --- 1. Validation and setup ----------------------------------------------
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("The 'purrr' package is required. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning(
      "The 'ggplot2' package is required for plotting but is not installed. 'plot_domain' set to FALSE.",
      call. = FALSE
    )
    plot_domain <- FALSE
  }

  if (!is.numeric(const_min) || !is.numeric(const_max) ||
      length(const_min) != 1L || length(const_max) != 1L ||
      !is.finite(const_min) || !is.finite(const_max)) {
    stop("'const_min' and 'const_max' must be finite numeric scalars.", call. = FALSE)
  }
  if (const_max <= const_min) {
    stop("'const_max' must be greater than 'const_min'.", call. = FALSE)
  }
  if (!is.function(f)) {
    stop("'f' must be a function of two variables.", call. = FALSE)
  }
  if (!is.function(limit1) || !is.function(limit2)) {
    stop("'limit1' and 'limit2' must be functions of one variable.", call. = FALSE)
  }
  if (!region_type %in% c("type1", "type2")) {
    stop("Invalid 'region_type'. Must be 'type1' or 'type2'.", call. = FALSE)
  }
  if (!is.numeric(n_outer) || !is.numeric(n_inner) ||
      length(n_outer) != 1L || length(n_inner) != 1L) {
    stop("'n_outer' and 'n_inner' must be numeric scalars.", call. = FALSE)
  }
  n_outer <- as.integer(n_outer)
  n_inner <- as.integer(n_inner)
  if (n_outer <= 0L || n_inner <= 0L ||
      n_outer %% 2L != 0L || n_inner %% 2L != 0L) {
    stop("'n_outer' and 'n_inner' must be positive even integers for Simpson's Rule.", call. = FALSE)
  }

  # Simpson coefficients
  coef_simp_inner <- c(1, rep(c(4, 2), n_inner / 2L - 1L), 4, 1)
  coef_simp_outer <- c(1, rep(c(4, 2), n_outer / 2L - 1L), 4, 1)

  # --- 2. Inner integration (Simpson) ---------------------------------------
  integral_inner <- function(val_outer, n_inner) {
    inner_min <- limit1(val_outer)
    inner_max <- limit2(val_outer)

    if (!is.numeric(inner_min) || !is.numeric(inner_max) ||
        length(inner_min) != 1L || length(inner_max) != 1L ||
        !is.finite(inner_min) || !is.finite(inner_max)) {
      stop("Inner limits must be finite numeric scalars.", call. = FALSE)
    }

    if (inner_max <= inner_min) {
      return(0)
    }

    h_inner   <- (inner_max - inner_min) / n_inner
    inner_vals <- seq(inner_min, inner_max, length.out = n_inner + 1L)

    if (region_type == "type1") {
      # outer is x, inner is y
      f_vals <- purrr::map_dbl(inner_vals, function(y) f(val_outer, y))
    } else {
      # region_type == "type2": outer is y, inner is x
      f_vals <- purrr::map_dbl(inner_vals, function(x) f(x, val_outer))
    }

    if (any(!is.finite(f_vals))) {
      stop("Function 'f' returned non-finite values inside the domain.", call. = FALSE)
    }

    sum(f_vals * coef_simp_inner) * (h_inner / 3)
  }

  # --- 3. Outer integration (Simpson) ---------------------------------------
  h_outer   <- (const_max - const_min) / n_outer
  outer_vals <- seq(const_min, const_max, length.out = n_outer + 1L)

  integral_outer_results <- purrr::map_dbl(
    outer_vals,
    integral_inner,
    n_inner = n_inner
  )

  if (any(!is.finite(integral_outer_results))) {
    stop("Non-finite values encountered in outer integration.", call. = FALSE)
  }

  integral_final <- sum(integral_outer_results * coef_simp_outer) * (h_outer / 3)

  # --- 4. Domain plotting (ggplot2) -----------------------------------------
  domain_plot <- NULL
  if (isTRUE(plot_domain)) {
    outer_seq <- seq(const_min, const_max, length.out = 100L)
    inner_min_vals <- limit1(outer_seq)
    inner_max_vals <- limit2(outer_seq)

    if (!is.numeric(inner_min_vals) || !is.numeric(inner_max_vals) ||
        length(inner_min_vals) != length(outer_seq) ||
        length(inner_max_vals) != length(outer_seq)) {
      stop(
        "Functions 'limit1' and 'limit2' must be vectorizable over the outer variable.",
        call. = FALSE
      )
    }

    data_plot <- data.frame(
      outer_var = outer_seq,
      inner_min = inner_min_vals,
      inner_max = inner_max_vals
    )

    if (region_type == "type1") {
      # Type I: x in [a, b], y in [h1(x), h2(x)]
      domain_plot <- ggplot2::ggplot(data_plot, ggplot2::aes(x = outer_var)) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = inner_min, ymax = inner_max),
          fill = "skyblue", alpha = 0.5
        ) +
        ggplot2::geom_line(
          ggplot2::aes(y = inner_min),
          color = "blue", linetype = "dashed"
        ) +
        ggplot2::geom_line(
          ggplot2::aes(y = inner_max),
          color = "blue", linetype = "dashed"
        ) +
        ggplot2::labs(
          title = "Domain Plot (Type I: dy dx)",
          x = "x (outer variable)",
          y = "y (inner variable)"
        ) +
        ggplot2::theme_minimal()
    } else {
      # Type II: y in [c, d], x in [h1(y), h2(y)]
      domain_plot <- ggplot2::ggplot(data_plot, ggplot2::aes(y = outer_var)) +
        ggplot2::geom_ribbon(
          ggplot2::aes(xmin = inner_min, xmax = inner_max),
          fill = "coral", alpha = 0.5
        ) +
        ggplot2::geom_line(
          ggplot2::aes(x = inner_min),
          color = "red", linetype = "dashed"
        ) +
        ggplot2::geom_line(
          ggplot2::aes(x = inner_max),
          color = "red", linetype = "dashed"
        ) +
        ggplot2::labs(
          title = "Domain Plot (Type II: dx dy)",
          x = "x (inner variable)",
          y = "y (outer variable)"
        ) +
        ggplot2::theme_minimal()
    }
  }

  # --- 5. Return -------------------------------------------------------------
  list(
    integral_value = integral_final,
    domain_plot    = domain_plot
  )
}

# -------------------------------------------------------------------
# Block of global variables (added here) to avoid NOTE in CRAN checks.

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "x", "y", "outer_var", "inner_min", "inner_max",
    "x_min", "x_max", "y_min", "y_max"
  ))
}
