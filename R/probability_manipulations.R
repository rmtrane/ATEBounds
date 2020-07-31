#' Calculate thetas from joint conditionals
#'
#' Calculates P(X = 1 | Z = z) from a table of P(Y = y, X = x | Z = z)
#'
#' @param tabp a 2x2xn_z_levels table given P(Y = y, X = x | Z = z), y = 0,1, x = 0,1, and z = 0,1,..., n_z_levels - 1.
#'
#' @export
thetas_from_tabp <- function(tabp) apply(tabp, c(1, 3), sum)[2, ]

#' Calculate gammas from joint conditionals
#'
#' Calculates \eqn{P(Y = 1 | Z = z)} from a table of \eqn{P(Y = y, X = x | Z = z)}
#'
#' @param tabp a \eqn{2x2xk} table given \eqn{P(Y = y, X = x | Z = z), y = 0,1, x = 0,1}, and \eqn{z = 0,1,..., k - 1}.
#'
#' @export
gammas_from_tabp <- function(tabp) apply(tabp, c(2, 3), sum)[2, ]

#' Calculate gammas/thetas from full data
#'
#' Calculates P(Y = 1 | Z = z) (`gammas`) and P(X = 1 | Z = z) (`thetas`) from data.frame/tibble with columns X,Y,Z
#'
#' @param dat data.frame with columns containing observations of X,Y,Z
#' @param X variable with observations on X
#' @param Y variable with observations on Y
#' @param Z variable with observations on Z
#'
#' @return A list with two vectors: `gammas` of length number of different Z values giving P(Y = 1 | Z = z), and `thetas` of length number of different Z values giving P(X = 1 | Z = z)
#'
#' @export
summary_probs_from_data <- function(dat, X = X, Y = Y, Z = Z){
  tmp <- dat %>%
    select(X = {{X}}, Y = {{Y}}, Z = {{Z}}) %>%
    pivot_longer(cols = c(X, Y),
                 names_to = "variable") %>%
    count(Z, variable, value) %>%
    group_by(Z, variable) %>%
    mutate(p = n/sum(n)) %>%
    filter(value == 1) %>%
    arrange(variable, Z)

  return(list(gammas = filter(tmp, variable == "Y")$p,
              thetas = filter(tmp, variable == "X")$p))
}
