#' Optimality check with Lagrange multipliers and bordered Hessian
#'
#' @description
#' Evaluates first- and second-order optimality conditions for a smooth
#' constrained optimization problem with equality constraints at a given
#' candidate point. The function checks the Lagrange conditions,
#' builds the bordered Hessian, and classifies the candidate as a
#' minimum, maximum or indeterminate/saddle according to the signs of
#' the leading principal minors.
#'
#' @details
#' Consider a problem of minimizing or maximizing a scalar function
#' \code{f(x)} subject to \code{m} equality constraints collected in
#' \code{g(x) = 0}, where \code{x} is a vector in R^n and
#' \code{g(x)} is a vector in R^m.
#'
#' At the candidate point \code{x}, the function:
#' \itemize{
#'   \item Approximates the gradient of \code{f} and the gradients
#'         of each constraint using second-order central finite
#'         differences.
#'   \item Builds the Jacobian matrix \code{J} of the constraints
#'         (rows are gradients of each constraint).
#'   \item Approximates the Hessian matrix of \code{f} and the Hessian
#'         of each constraint, also by central finite differences.
#'   \item Forms the Hessian of the Lagrangian by combining the Hessian
#'         of \code{f} and the Hessians of the constraints with the
#'         Lagrange multipliers.
#'   \item Builds the bordered Hessian matrix using the Jacobian and
#'         the Hessian of the Lagrangian.
#'   \item Computes leading principal minors of the bordered Hessian
#'         and uses their signs to classify the candidate.
#' }
#'
#' The classification is based on the standard bordered Hessian test:
#' after multiplying each leading principal minor by (-1)^m, if all
#' resulting values are positive the point is classified as a minimum;
#' if their signs alternate (negative, positive, negative, and so on),
#' the point is classified as a maximum. In any other case, the result
#' is reported as indeterminate.
#' All derivatives (gradients and Hessians) are computed numerically
#' with central finite differences of second order. The step sizes can
#' be given explicitly or chosen automatically from the scale of
#' the point \code{x}.
#'
#' @param f Objective function. It must be \code{function(x)} and
#'   return a single numeric value. The argument \code{x} is a
#'   numeric vector of length \code{n}.
#' @param g Equality constraints. Either a single function \code{function(x)}
#'   returning a numeric vector of length \code{m}, or a list of scalar
#'   functions \code{function(x)}, one per constraint.
#' @param x Numeric vector giving the candidate point at which the
#'   optimality conditions are evaluated.
#' @param h Step size for finite differences. It can be:
#'   \itemize{
#'     \item a single numeric value used for all coordinates,
#'     \item a numeric vector of length \code{n} with one step per
#'           coordinate,
#'     \item or \code{NULL}, in which case step sizes are chosen as
#'           \code{1e-4 * (1 + abs(x))} componentwise.
#'   }
#' @param tol Numeric tolerance used to judge feasibility of the
#'   constraints, the effective rank of the Jacobian, near singularity
#'   of matrices and very small principal minors.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{ok_stationarity}: numeric value of the norm of the
#'         stationarity residual. When constraints are present, this
#'         measures how close the gradient of \code{f} is to the
#'         linear combination given by the Jacobian transpose and the
#'         Lagrange multipliers.
#'   \item \code{ok_feasible}: maximum absolute value of the
#'         constraint vector \code{g(x)} at the candidate point.
#'   \item \code{lambda}: numeric vector of length \code{m} with the
#'         Lagrange multipliers.
#'   \item \code{J}: Jacobian matrix of the constraints at \code{x},
#'         with dimension \code{m x n}.
#'   \item \code{H_f}: Hessian matrix of the objective function
#'         at \code{x}, of size \code{n x n}.
#'   \item \code{H_g}: list of Hessian matrices corresponding to each
#'         constraint function, each of size \code{n x n}.
#'   \item \code{H_L}: Hessian matrix of the Lagrangian at \code{x}.
#'   \item \code{B}: bordered Hessian matrix, of size
#'         \code{(m + n) x (m + n)} when constraints are present.
#'   \item \code{minors}: \code{data.frame} with one row per leading
#'         principal minor. It contains the order of the minor, its
#'         signed value and the logarithm of the absolute determinant
#'         used in the computation.
#'   \item \code{clasificacion}: character value equal to
#'         \code{"minimo"}, \code{"maximo"} or \code{"indeterminado"},
#'         according to the bordered Hessian criterion.
#'   \item \code{notas}: character vector with diagnostic messages
#'         about the rank of the Jacobian, near singularity of matrices
#'         or any numerical issues detected.
#' }
#'
#' @examples
#' ## 1) Minimum with one constraint:
#' ##    f(x,y) = x^2 + y^2,  g(x,y) = x + y - 1 = 0  -> (0.5, 0.5)
#' f1 <- function(x) x[1]^2 + x[2]^2
#' g1 <- function(x) x[1] + x[2] - 1
#' lagrange_check(f1, g1, x = c(0.5, 0.5))
#'
#' ## 2) Maximum with one constraint:
#' ##    f(x,y) = -(x^2 + y^2),  g(x,y) = x + y - 1 = 0  -> (0.5, 0.5)
#' f2 <- function(x) -(x[1]^2 + x[2]^2)
#' lagrange_check(f2, g1, x = c(0.5, 0.5))
#'
#' ## 3) Two constraints in R^3 (minimum norm with two planes)
#' f3 <- function(x) sum(x^2)
#' g3 <- list(
#'   function(x) x[1] + x[2] + x[3] - 1,
#'   function(x) x[1] - x[3]
#' )
#' ## Candidate solution: x1 = x3, 2*x1 + x2 = 1  ->  x = (1/3, 1/3, 1/3)
#' lagrange_check(f3, g3, x = c(1, 1, 1) / 3)
#'
#' @export
lagrange_check <- function(f, g, x, h = NULL, tol = 1e-6) {
  # ---------- internal helpers (not exported) -------------------------------
  .steps <- function(x, h) {
    n <- length(x)
    if (is.null(h)) {
      hi <- 1e-4 * (1 + abs(x))
    } else if (length(h) == 1L) {
      hi <- rep(as.numeric(h), n)
    } else {
      hi <- as.numeric(h)
      if (length(hi) != n) {
        stop("'h' must be scalar, length n, or NULL.", call. = FALSE)
      }
    }
    if (any(!is.finite(hi)) || any(hi <= 0)) {
      stop("Finite-difference steps 'h' must be positive and finite.", call. = FALSE)
    }
    hi
  }

  .as_list_g <- function(g, x) {
    if (is.null(g)) return(list())
    if (is.list(g)) return(g)
    if (is.function(g)) {
      out <- g(x)
      if (!is.numeric(out)) stop("'g(x)' must be numeric.", call. = FALSE)
      m <- length(out)
      return(lapply(seq_len(m), function(i) function(xx) g(xx)[i]))
    }
    stop("'g' must be function(x) or a list of function(x).", call. = FALSE)
  }

  .num_grad <- function(fun, x, h) {
    n  <- length(x)
    hi <- .steps(x, h)
    g  <- numeric(n)
    for (i in seq_len(n)) {
      ei <- rep(0, n); ei[i] <- hi[i]
      g[i] <- (fun(x + ei) - fun(x - ei)) / (2 * hi[i])
    }
    g
  }

  .num_hess <- function(fun, x, h) {
    n  <- length(x)
    hi <- .steps(x, h)
    H  <- matrix(0, n, n)
    fx <- fun(x)
    # diagonal
    for (i in seq_len(n)) {
      ei <- rep(0, n); ei[i] <- hi[i]
      fpp <- fun(x + ei); fmm <- fun(x - ei)
      H[i, i] <- (fpp - 2 * fx + fmm) / (hi[i]^2)
    }
    # off-diagonal
    if (n >= 2L) {
      for (i in 1:(n - 1L)) for (j in (i + 1L):n) {
        ei <- rep(0, n); ej <- rep(0, n)
        ei[i] <- hi[i]; ej[j] <- hi[j]
        fpp <- fun(x + ei + ej)
        fpm <- fun(x + ei - ej)
        fmp <- fun(x - ei + ej)
        fmm <- fun(x - ei - ej)
        val <- (fpp - fpm - fmp + fmm) / (4 * hi[i] * hi[j])
        H[i, j] <- H[j, i] <- val
      }
    }
    H
  }

  .det_sign_logabs <- function(M) {
    d <- suppressWarnings(determinant(M, logarithm = TRUE))
    sign   <- as.numeric(d$sign)      # 0 if singular
    logabs <- as.numeric(d$modulus)   # log |det|
    list(sign = sign, logabs = logabs)
  }

  .classify <- function(m, deltas, tol) {
    # Classification from signs of (-1)^m * Delta_k
    if (length(deltas) == 0L) return("indeterminado")
    if (any(abs(deltas) < tol)) return("indeterminado")
    s <- (-1)^m * deltas
    if (all(s > 0)) return("minimo")
    alt <- all(s[seq(1, length(s), 2)] < 0) &&
      all(s[seq(2, length(s), 2)] > 0)
    if (alt) return("maximo")
    "indeterminado"
  }

  # ---------- validations ----------------------------------------------------
  if (!is.function(f)) {
    stop("'f' must be function(x) returning a scalar.", call. = FALSE)
  }
  if (!is.numeric(x) || any(!is.finite(x))) {
    stop("'x' must be a finite numeric vector.", call. = FALSE)
  }
  x <- as.numeric(x)
  n <- length(x)

  if (!is.null(h) &&
      (!is.numeric(h) || any(!is.finite(h)) || any(h <= 0))) {
    stop("'h' must be positive, finite numeric (or NULL).", call. = FALSE)
  }

  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("'tol' must be a positive finite numeric scalar.", call. = FALSE)
  }

  glist <- .as_list_g(g, x)
  m     <- length(glist)

  # ---------- gradient of f, Jacobian J, Hessians ----------------------------
  gvals <- if (m > 0L) {
    vapply(glist, function(gi) gi(x), numeric(1))
  } else {
    numeric(0)
  }

  if (m > 0L && any(!is.finite(gvals))) {
    stop("Constraint functions in 'g' returned non-finite values.", call. = FALSE)
  }

  grad_f <- .num_grad(f, x, h)
  H_f    <- .num_hess(f, x, h)

  if (any(!is.finite(grad_f))) {
    stop("Gradient of 'f' contains non-finite values.", call. = FALSE)
  }
  if (any(!is.finite(H_f))) {
    stop("Hessian of 'f' contains non-finite values.", call. = FALSE)
  }

  if (m > 0L) {
    J   <- matrix(NA_real_, nrow = m, ncol = n)
    H_g <- vector("list", m)
    for (i in seq_len(m)) {
      gi <- glist[[i]]
      J[i, ]   <- .num_grad(gi, x, h)
      H_g[[i]] <- .num_hess(gi, x, h)
    }
    if (any(!is.finite(J))) {
      stop("Jacobian 'J' contains non-finite values.", call. = FALSE)
    }
  } else {
    J   <- matrix(0, 0L, n)
    H_g <- list()
  }

  # ---------- lambda from J J^T lambda = J grad_f ----------------------------
  notas <- character()
  if (m > 0L) {
    sv  <- svd(J)
    rJ  <- sum(sv$d > tol * max(1, sv$d[1]))
    if (rJ < m) {
      notas <- c(notas,
                 sprintf("Warning: rank(J) = %d < m = %d; regularity fails.", rJ, m))
    }

    rhs <- as.numeric(J %*% grad_f)
    JJt <- J %*% t(J)

    # Check near-singularity of JJt via determinant
    dJJ <- .det_sign_logabs(JJt)
    if (dJJ$sign == 0 || exp(dJJ$logabs) < tol) {
      notas <- c(notas,
                 "JJ^T nearly singular; using pseudoinverse (SVD).")
      sv2  <- svd(JJt)
      invd <- ifelse(sv2$d > tol, 1 / sv2$d, 0)
      pinv <- sv2$v %*% (t(sv2$u) * invd)
      lambda <- as.numeric(pinv %*% rhs)
    } else {
      lambda <- as.numeric(solve(JJt, rhs))
    }
  } else {
    lambda <- numeric(0)
  }

  # ---------- Hessian of the Lagrangian -------------------------------------
  H_L <- H_f
  if (m > 0L) {
    for (i in seq_len(m)) {
      H_L <- H_L - lambda[i] * H_g[[i]]
    }
  }

  if (any(!is.finite(H_L))) {
    notas <- c(notas, "Hessian of the Lagrangian contains non-finite values.")
  }

  # ---------- Bordered Hessian B --------------------------------------------
  if (m > 0L) {
    B <- rbind(
      cbind(matrix(0, m, m), J),
      cbind(t(J), H_L)
    )
  } else {
    B <- H_L
    notas <- c(notas,
               "No constraints: using H_L (Hessian of f) for the standard test.")
  }

  # ---------- Leading principal minors Delta_k -------------------------------
  deltas <- numeric(n)
  signs  <- integer(n)
  logs   <- numeric(n)

  for (k in seq_len(n)) {
    idx <- seq_len(m + k)
    subM <- B[idx, idx, drop = FALSE]
    d    <- .det_sign_logabs(subM)
    signs[k] <- d$sign
    logs[k]  <- d$logabs
    deltas[k] <- if (d$sign == 0) 0 else d$sign * exp(d$logabs)
  }

  clas <- .classify(m, deltas, tol)

  ok_stationarity <- if (m > 0L) {
    resid <- grad_f - as.vector(t(J) %*% lambda)
    sqrt(sum(resid * resid))
  } else {
    sqrt(sum(grad_f * grad_f))
  }

  ok_feasible <- if (m > 0L) max(abs(gvals)) else 0

  list(
    ok_stationarity = ok_stationarity,
    ok_feasible     = ok_feasible,
    lambda          = lambda,
    J               = J,
    H_f             = H_f,
    H_g             = H_g,
    H_L             = H_L,
    B               = B,
    minors          = data.frame(
      k         = seq_len(n),
      order     = m + seq_len(n),
      Delta_k   = deltas,
      sign      = signs,
      logAbsDet = logs
    ),
    clasificacion = clas,
    notas         = notas
  )
}
