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

#' Calculate ATE from simulated data
#'
#' Calculates the ATE based on the results from ATEBounds::simulate_data.
#'
#' @param from_simulate_data the object returned by ATEBounds::simulate_data
#'
#' @return a numeric value giving the ATE = P(Y = 1 | X = 1, U) - P(Y = 1 | X = 0, U). To make sure this
#' is accurate, use large n in your ATEBounds::simulate_data call.
#'
#' @export
ATE_from_simulated_data <- function(from_simulate_data){
  intercept <- filter(from_simulate_data$coefficients, effect == "Yintercept")$coef
  x_beta <- filter(from_simulate_data$coefficients, effect == "X_on_Y")$coef
  u_beta <- filter(from_simulate_data$coefficients, effect == "U_on_Y")$coef

  pY1X0 <- 1 / (1 + exp(-intercept - u_beta*from_simulate_data$simulated_data$U))
  pY1X1 <- 1 / (1 + exp(-intercept - x_beta - u_beta*from_simulate_data$simulated_data$U))

  return(mean(pY1X1 - pY1X0))
}
