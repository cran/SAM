#' Unified Entry Point for Sparse Additive Modelling
#'
#' Fit a sparse additive model by dispatching to the appropriate family-specific
#' function (\code{samQL}, \code{samLL}, \code{samEL}, or \code{samHL}).
#'
#' @param X Numeric training matrix with \code{n} rows (samples) and \code{d}
#'   columns (features).
#' @param y Response vector of length \code{n}.
#' @param p The number of basis spline functions. The default value is 3.
#' @param family A string specifying the loss family. One of \code{"gaussian"}
#'   (default), \code{"binomial"}, \code{"poisson"}, or \code{"hinge"}.
#' @param \dots Additional arguments passed to the family-specific function.
#' @return An S3 object of class \code{samQL}, \code{samLL}, \code{samEL}, or
#'   \code{samHL}, depending on the chosen family.
#' @seealso \code{\link{samQL}}, \code{\link{samLL}}, \code{\link{samEL}},
#'   \code{\link{samHL}}
#' @examples
#' n <- 100; d <- 50
#' X <- matrix(runif(n * d), n, d)
#' y <- rnorm(n)
#' fit <- sam(X, y, family = "gaussian")
#' fit
#' @export
sam <- function(X, y, p = 3, family = c("gaussian", "binomial", "poisson", "hinge"), ...) {
  family <- match.arg(family)
  switch(family,
    gaussian = samQL(X, y, p = p, ...),
    binomial = samLL(X, y, p = p, ...),
    poisson  = samEL(X, y, p = p, ...),
    hinge    = samHL(X, y, p = p, ...)
  )
}
