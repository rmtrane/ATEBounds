#' Get P(X = 1 | Z_j), P(X) from a specified logistic model
#'
#' This code calculates \eqn{P(X = 1 | Z_j = z)}, \eqn{z = 0,1,2,...,k}, for any number of \eqn{Z_j}'s. Assuming
#' \eqn{P(X = 1 | Z_1 = z_1, ..., Z_n = z_n) = \text{expit}(\beta_0 + \sum_j \beta_j z_j)}, this code numerical integrates
#' out all but one IV for all IVs. Here, it is assumed that all IVs are independent. Values of \eqn{\beta}'s and \eqn{P(Z_j = z), z=0,1,...,k}
#' are user specified.
#'
#' This code is fast if all \eqn{\left\{P(Z_j = z)\right\}_{z=0}^k} are the same for all \eqn{j}, and all coefficients \eqn{\beta_j}'s are the same.
#' In this case, it can handle an arbitrary number of IVs.
#'
#' For other cases, it will brute force the integral to compute exactly when possible. When number of IVs is greater than 15, it will
#' switch to stochastic integration (essentially MCMC).
#'
#' @param b0 numeric; value of intercept \eqn{\beta_0}.
#' @param b a numeric vector of slope coefficients \eqn{\beta_1, ..., \beta_n} in logit probability model. Length of vector specifies number of IVs.
#' @param p.z a list of vectors of probabilities of \eqn{P(Z_j = z), z = 0,1,2,...,k}. Length of list must be same as
#'   length of `b`, and each element of of the list must have length \eqn{l}, number of categories of the IVs.
#' @param nSim number of simulations to use when using stochastic integration
#'
#' @return tibble with one row for each IV, and a total of k+2 columns: one column to id the IVs (`j`), one column for each \eqn{z} giving
#'   \eqn{P(X = 1 | Z_j = z)}, and one column giving \eqn{P(X = 1)}.
#'
#' @importFrom tibble as_tibble
#'
#' @export
marginals_from_logit <- function(b0, b, p.z = rep(list(c(1, 1, 1)/3),3), nSim = 10000) {
  if(is.null(b0) | is.null(b) | is.null(p.z))
    stop("b0, b, and p.z must all be specified")

  if(!is.numeric(b0) | !is.numeric(b) | !all(unlist(lapply(p.z, is.numeric))))
    stop("both b0, b, and p.z must all be numeric")

  if(length(b) < 3)
    stop("b must be of length >= 3")

  if(any(unlist(lapply(p.z, min)) < 0) | any(unlist(lapply(p.z, max)) > 1) | any(unlist(lapply(p.z, sum)) != 1))
    stop("Each element of p.z must be a vector with values between 0 and 1, and sum(p.z[[i]]) must be equal to 1 for all i.")

  if(length(unique(unlist(lapply(p.z, length)))) != 1)
    stop("All elements of p.z must have same length!")

  # Matrix to hold results
  margProb = matrix(0,length(b),length(p.z[[1]]) + 1)

  # Appropriate column names
  colnames(margProb) = c(paste0("P(X = 1 | Z_j = ", 0:(length(p.z[[1]]) - 1),")"), "P(X = 1)")

  # Appropriate row names
  rownames(margProb) = 1:length(b)

  # If all intercepts and probabilities are identical, things are simpler.
  if(length(unique(b)) == 1 && length(unique(p.z)) == 1){
    # In this case, (Z_1,...,Z_n) ~ Multinomial with size k, P(Z_1 = j) = p.z[[1]][j].

    # Need to find all possible combinations of outcomes of IVs except 1. Basically,
    # we distribute the n IVs into the k+1 categories 0,1,...,k.

    # First, all combinations of number of IVs in the different categories.
    all_combos <- expand.grid(rep(list(0:(length(b)-1)), length(p.z[[1]])))

    # Restrict to combinations where the total number of IVs is n-1.
    all_combos <- all_combos[rowSums(all_combos) == length(b) - 1, ]

    # Now, each row of all_combos gives one way of distributing the n IVs in the k+1
    # values 0,1,2,...,k. For example, the row c(2, 3, 0) correspdonds to a scenario
    # where we have five three-leveled IVs where two take the value 0, three take the
    # value 1, and none take the value 0.

    # Since all slopes are equal, \sum b_i Z_i = b_1 \sum Z_i. The matrix all_combos
    # has the number of IVs with outcome 0 in the first column, 1 in the second column,
    # etc. I.e. for each row x, \sum Z_i is the inner product of the row x and the numbers
    # 0,1,2,..., k (k is the number of categories, i.e. number of columns of all_combos, i.e.
    # length of each row vector.)
    sums <- apply(all_combos, 1, function(x) sum(x * 0:(length(x) - 1)))

    # Using dmultinom, we find probability of each outcome.
    probs <- apply(all_combos, 1, function(x) dmultinom(x = x, size = sum(x), prob = p.z[[1]]))

    # Get linear combination of Z's and coefficients, in this case simply b * \sum Z_j
    temp <- b[1] * sums

    # Then for each IV...
    for(i in 1:length(b)){
      # ... we find the marginal probabilities P(X = 1 | Z_i = 0), P(X = 1 | Z_i = 1), ..., P(X = 1 | Z_i = k).
      margProb[i,-ncol(margProb)] <- sapply(0:(length(p.z[[1]])-1),
                                            # This functions finds P(X = 1 | Z_i = z), and is applied to the vector c(0,1,...,k).
                                            function(z){
                                              # \sum_{Z_j: j\neq i} P(X = 1 | Z_1 = z_1, ..., Z_n = Z_n) * P(Z_1 = z_1, ..., Z_{j-1} = z_{j-1}, Z_{j+1} = z_{j+1}, ..., Z_n = z_n)
                                              sum (1/(1 + exp(-(b0 + b[i]*z + temp))) * probs )
                                            })
      margProb[i,ncol(margProb)] <- sum(margProb[i, -ncol(margProb)] * p.z[[1]])
    }
  } else {
    if(length(p.z) > 15) {
      message("You'll run out of memory. Switching to MCMC approx...")
      for(i in 1:length(b)) {
        # Sample "observed" values of other Z's.
        Z_minus_i <- sapply(seq_along(b)[-i], FUN = function(y) sample(x = 0:(length(p.z[[1]])-1), size = nSim, replace = TRUE, prob = p.z[[y]]))

        # Linear combination of coefficients and "observed" values of Z's.
        temp = as.numeric( Z_minus_i %*% b[-i])


        margProb[i,-ncol(margProb)] <- sapply(1:length(p.z[[1]]) - 1, function(x) mean(1/(1 + exp(-(b0 + b[i]*x + temp)))))
        margProb[i, ncol(margProb)] <- sum(margProb[i,-ncol(margProb)] * p.z[[i]])
      }
    } else  {
      # All possible combinations of observed values of Z_1, Z_2, ..., Z_{i-1}, Z_{i+1}, ..., Z_k.
      Z_minus_i <- expand.grid(lapply(p.z[-1], function(x) seq_along(x)-1))

      # For each IV...
      for(i in 1:length(b)) {
        # ... calculate probability of the other IVs being any given combination
        probs <- sapply(1:ncol(Z_minus_i),
                        # P(combination)
                        FUN = function(x) p.z[-i][[x]][Z_minus_i[, x] + 1])

        # ... and take product to get joint probability.
        probs <- apply(probs, 1, prod)

        # Linear combination of coefficients and "observed" values.
        temp <- as.numeric( as.matrix(Z_minus_i) %*% b[-i])

        # Marginalize as earlier
        margProb[i,-ncol(margProb)] <- sapply(1:length(p.z[[1]]) - 1, function(x) sum (1/(1 + exp(-(b0 + b[i]*x + temp))) * probs))
        margProb[i, ncol(margProb)] <- sum(margProb[i,-ncol(margProb)] * p.z[[i]])
      }
    }
  }
  return(as_tibble(margProb, rownames = "j"))
}
