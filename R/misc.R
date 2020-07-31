#' Simulate dependent categorical variables
#'
#' @export
rbivNorm <- function(n, mu = 0, var1 = 1, var2 = 1, rho = 0, probs = c(0.25, 0.5, 0.25)) {
  tmp <- MASS::mvrnorm(n = n,
                       mu = rep(mu, 2),
                       Sigma = matrix(data = c(var1, rho*sqrt(var1)*sqrt(var2),
                                               rho*sqrt(var1)*sqrt(var2), var2),
                                      nrow = 2, byrow = FALSE))

  tmp <- matrix(tmp, nrow = n)
  out <- apply(tmp, 2,
               function(x){
                 as.numeric(cut(x, breaks = qnorm(unique(c(0, cumsum(probs), 1))))) - 1
               })

  return(out)
}

