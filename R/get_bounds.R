#' Function to return ACE bounds from P(X = 1 | Z), P(Y = 1 | Z), or P(X,Y | Z).
#'
#' Using the [polymake](https://www.polymake.org) program, a matrix \eqn{A} is created such that \eqn{A \cdot p} will have non-negative
#' elements if the probabilities \eqn{p} specify a probability distribution that does not violate the IV model. \eqn{p} can be probabilities
#' of either P(X = x | Z = z) and P(Y = y | Z = z), or P(Y = y, X = x | Z = z).
#'
#' @param gammas vector of length \eqn{k} giving \eqn{P(Y = 1 | Z = z), z = 0, ..., (k - 1)}
#' @param thetas vector of length \eqn{k} giving \eqn{P(X = 1 | Z = z), z = 0, ..., (k - 1)}
#' @param zetas array of dimensions \eqn{2\cdot 2\cdot k} giving P(Y = y, X = x | Z = z) at position (y+1,x+1,z+1),
#' @param stop logical; should function stop with error if constraints are not met?
#' @param warning logical; if stop is \code{FALSE}, this specifies whether a warning should be returned or not. Default: \code{TRUE}
#' @param from_polymake result of \code{read_polymake_results}. If \code{NULL} (default), \code{x_mono} and \code{y_mono} must be specified,
#'   and previously created matrices will be used.
#' @param x_mono logical; should P(X = 1 | Z = z) be assumed monotonically increasing? (Default: FALSE)
#' @param y_mono logical; should P(Y = 1 | Z = z) be assumed monotonically increasing? (Default: FALSE)
#'
#' @return A list with the following elements:
#'
#' * `lower_bounds`: a vector of values that the lower bound must be greater than
#' * `upper_bounds`: a vector of values that the upper bound must be smaller than
#' * `interval`: a vector with the lower bound (`max(lower_bounds)`) and upper bound (`min(upper_bounds)`).
#' * `width`: simply the width of the interval
#' * `constraints_violated`: a logical value that indicates if any constraints are violated (`TRUE`) or not (`FALSE`)
#' * `assumptions`: logical vector of length two giving `x_mono` and `y_mono`. For bookkeeping.
#'
#' @export
get_bounds <- function(gammas = NULL, thetas = NULL,
                       zetas = NULL,
                       stop = TRUE,
                       warning = TRUE,
                       from_polymake = NULL,
                       x_mono = FALSE,
                       y_mono = FALSE){

  if(is.null(gammas) + is.null(thetas) == 1)
    stop("If one of gammas or thetas is not NULL, the other must also be specified.")

  if(any(is.null(gammas) + is.null(zetas) + is.null(thetas) == c(0,3)))
    stop("Either zetas, or gammas and thetas must be non-null. zetas =", is.null(zetas), "gammas =", is.null(gammas), "thetas =", is.null(thetas))

  if(length(gammas) != length(thetas))
    stop("gammas and thetas must be same length.")

  if(!is.null(from_polymake) && (length(thetas) != (ncol(from_polymake) - 1)/4) & is.null(zetas))
    stop(paste("Length of gammas and thetas must both be (ncol(from_polymake)-1)/4 =", (ncol(from_polymake) - 1)/4))

  if(is.null(from_polymake) & any(c(is.null(x_mono), is.null(y_mono))))
    stop("If from_polymake is not provided, x_mono AND y_mono must be specified.")

  if(is.null(zetas)){
    probs <- c(1-gammas, gammas,
               1-thetas, thetas)

    if(is.null(from_polymake)){
      from_polymake <- matrices_from_polymake %>%
        filter(data_format == "bivariate",
               n_z_levels == length(gammas),
               x_monotone == x_mono,
               y_monotone == y_mono) %>%
        pull(matrix)

      if(length(from_polymake) != 1)
        stop("The number of categories is not supported.")

      from_polymake <- from_polymake[[1]]
    }
  } else {
    if(any(abs(apply(zetas, 3, sum) - 1) > 10^-10))
      stop("zetas not valid. Must sum to one for each level of Z")

    k <- length(zetas)/4

    probs <- c(zetas[x="0",y="0",1:k],
               zetas[x="1",y="0",1:k],
               zetas[x="0",y="1",1:k],
               zetas[x="1",y="1",1:k])

    if(is.null(from_polymake)){
      from_polymake <- matrices_from_polymake %>%
        filter(data_format == "trivariate",
               n_z_levels == k,
               x_monotone == x_mono,
               y_monotone == y_mono) %>%
        pull(matrix)

      if(length(from_polymake) != 1)
        stop("The number of categories is not supported.")

      from_polymake <- from_polymake[[1]]
    }
  }

  ## Check constrains
  constraints_violated <- from_polymake %>%
    filter(alpha == 0) %>%
    select(-alpha) %>%
    as.matrix() %*% probs < 0

  if(sum(constraints_violated) > 0 & stop)
    stop("The probability distributions specified do not satisfy the IV constraints")

  if(sum(constraints_violated) > 0 & warning & !stop){
    warning("The probability distributions specified do not satisfy the IV constraints")
  }

  upper_bounds <- from_polymake %>%
    filter(alpha == -1) %>%
    select(-alpha) %>%
    as.matrix(.) %*% probs

  lower_bounds <- from_polymake %>%
    filter(alpha == 1) %>%
    select(-alpha) %>%
    as.matrix(.) %*% probs

  combine_bounds <- list(lower_bounds = -c(lower_bounds),
                         upper_bounds = c(upper_bounds))

  actual_bounds <- list(interval = c(lower = max(combine_bounds$lower_bounds),
                                     upper = min(combine_bounds$upper_bounds)),
                        width = min(combine_bounds$upper_bounds) -
                          max(combine_bounds$lower_bounds))

  if (actual_bounds$width < 0 & stop){
    stop("Upper bound is smaller than lower bound.")
  }

  if (actual_bounds$width < 0 & warning & !stop){
    warning("Upper bound is smaller than lower bound.")
  }

  output <- c(actual_bounds, combine_bounds,
              list(assumptions = c(x_mono = x_mono, y_mono = y_mono)),
              constraints_violated = sum(constraints_violated) > 0)

  return(output)
}
