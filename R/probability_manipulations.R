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

#' Calculate gammas/thetas/zetas from full data
#'
#' Calculates either P(Y = 1 | Z = z) (`gammas`) and P(X = 1 | Z = z) (`thetas`), or P(X = x, Y = y | Z = z) (`zetas`) from data.frame/tibble with columns X,Y,Z
#'
#' @param dat data.frame with columns containing observations of X,Y,Z
#' @param X variable with observations on X
#' @param Y variable with observations on Y
#' @param Z variable with observations on Z
#' @param data_format bivariate if you want P(X = 1 | Z = z) and P(Y = 1 | Z = z), trivariate if you want P(X = x, Y = y | Z = z).
#'   If `NULL` (default), trivariate is chosen if both X and Y are
#'
#' @return a list with the following elements:
#'
#'   * if `data_format` is bivariate, a vector `gammas` of length number of different Z values giving P(Y = 1 | Z = z) is returned if `Y` specified, and `thetas` of length number of different Z values giving P(X = 1 | Z = z) if `X` specified.
#'   * if `data_format` is trivariate, a 2x2xn_z_leveles array `zetas` giving P(Y = y, X = x | Z = z) as entry (x+1, y+1, z+1).
#'
#' @export
probs_from_data <- function(dat, X, Y, Z, data_format){

  if(missing(Z))
     stop("You must provide a variable to use as IV")

  if(missing(data_format)){
    data_format <- if_else(!missing(X) & !missing(Y), "trivariate", "bivariate")
  }

  if(!data_format %in% c("bivariate", "trivariate"))
    stop("data_format must be either NULL, 'bivariate' or 'trivariate'")

  if(data_format == "trivariate" & (missing(X) | missing(Y)))
    stop('When specifying data_format = "trivariate", both X and Y must be specified.')

  if(missing(X) & missing(Y)){
    stop("You must provide at least one of X or Y")
  } else {
    if (missing(X))
      X <- NULL
    if (missing(Y))
      Y <- NULL
  }

  out <- list()

  if(data_format == "bivariate"){
    tmp <- dat %>%
      select(X = {{X}}, Y = {{Y}}, Z = {{Z}}) %>%
      pivot_longer(cols = c(X, Y),
                   names_to = "variable") %>%
      count(Z, variable, value) %>%
      group_by(Z, variable) %>%
      mutate(p = n/sum(n)) %>%
      filter(value == 1) %>%
      arrange(variable, Z)

    if("Y" %in% tmp$variable)
      out$gammas <- filter(tmp, variable == "Y")$p

    if("X" %in% tmp$variable)
      out$thetas <- filter(tmp, variable == "X")$p
  }

  if (data_format == "trivariate"){
    tmp <- dat %>%
      select(x = {{X}}, y = {{Y}}, z = {{Z}})

    tmp_out <- xtabs(~x+y+z,data = tmp)
    out$zetas <- tmp_out / rep(apply(tmp_out, 3, sum), each = 4)

  }

  return(out)
}
