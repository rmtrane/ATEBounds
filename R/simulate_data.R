#' Simulate data from IV model
#'
#' This function allows the user to simulate data from an IV model
#' with an arbitrary number of IVs with the option to include IVs
#' that are dependent on each other, and IVs with a direct effect
#' on the outcome.
#'
#' Data are simulated from a logistic model as follows:
#' \itemize{
#'   \item A confounder \eqn{U} is simulated based on the distribution specified.
#'   \item A user specified number of independent IVs (\eqn{Z_1, ..., Z_p}) that are either
#'     di- or trichotomous. Probabilities \eqn{P(Z_j = z), z=0,1} (and \eqn{P(Z_j = 2)} if trichotomous
#'     case is considered) are specified by the user.
#'   \item Optionally, a pair of dependent IVs, \eqn{Z_{D1}, Z_{D2}}, simulated as follows:
#'     \itemize{
#'       \item User specifies \eqn{\rho}
#'       \item \eqn{(N_{D1}, N_{D2}) \sim N\left(\begin{pmatrix} 0 \\ 0 \end{pmatrix},\begin{pmatrix} 1 & \rho \\ \rho & 1 \end{pmatrix} \right)},
#'         then \deqn{
#'           Z_{Di} = \left\{\begin{array}{cl} 0 & \text{ if } N_{Di} \le z_{p_1} \\ 1 & \text{ if } z_{p_0} < N_{Di} \le z_{p_0 + p_1} \\ 2 & \text{ otherwise} \end{array}  \right. ,
#'         } where \eqn{z_k} is the \eqn{k'th} percentile of the standard normal (\eqn{P(Z \le z_k) = k} when \eqn{Z \sim N(0,1)}).
#'     }
#'   \item A binary treatment is simulated as \eqn{X | Z_1,...,Z_p, Z_{D1}, Z_{D2}, U \sim \text{Bernoulli}(p_X)},
#'     where \deqn{\text{logit}(p_X) = \gamma_0 + \gamma_1 \cdot Z_1 + ... + \gamma_p \cdot Z_p + \gamma_{D1} \cdot Z_{D1} + \gamma_{D2} \cdot Z_{D2} + \gamma_U \cdot U},
#'     and all parameters \eqn{\gamma_j} are user specified.
#'   \item A binary outcome is simulated as \eqn{Y | X, U, Z_1,...,Z_p, Z_{D1}, Z_{D2}, U \sim \text{Bernoulli}(p_Y)},
#'     where \deqn{\text{logit}(p_Y) = \beta_0 + \beta_X \cdot X + \beta_1 \cdot Z_1 + ... + \beta_p \cdot Z_p + \beta_{D1} \cdot Z_{D1} + \beta_{D2} \cdot Z_{D2} + \beta_U \cdot U},
#'     and all parameters \eqn{\beta_j} are user specified.
#' }
#'
#' @param sample_size default 10000; number of observations to simulate.
#' @param IVs_ps list of numeric vectors. Length of list specifies number of independent IVs to include. Each vector
#'   gives probabilities of \eqn{P(Z = i)}, and must therefore sum to 1.
#' @param indIVs_on_X numeric vector specifying the effects of the independent IVs on X (i.e. \eqn{\gamma_1, ..., \gamma_n}).
#' @param indIVs_on_Y numeric vector specifying the effects of the independent IVs on Y (i.e. \eqn{\beta_1, ..., \beta_n}). Default: 0
#' @param U `distributions3` object specifying the distribution to be used for the unmeasured confounder. Default is `Bernoulli(0.5)`.
#' @param X_intercept numeric value to use for \eqn{\gamma_x}. Default: 0.05.
#' @param Y_intercept numeric value to use for \eqn{\gamma_y}. Default: 0.05.
#' @param U_on_X numeric value to use for \eqn{\alpha_U}. Default: 0.1.
#' @param U_on_Y numeric value to use for \eqn{\beta_U}. Default: 0.1.
#' @param X_on_Y numeric value to use for \eqn{\beta_X}. Default: 0.05.
#' @param rho numeric value specifying \eqn{\rho}. Default: NULL.
#' @param depIVs_on_X numeric vector specifying the effects of the dependent IVs on X (i.e. \eqn{\alpha_{D1}} and \eqn{\alpha_{D2}}).
#' @param depIVs_on_Y numeric vector specifying the effects of the dependent IVs on Y (i.e. \eqn{\beta_{D1}} and \eqn{\beta_{D2}}).
#' @param dep_probs numeric vector. The length specifies the number of categories for the IVs, and the values give \eqn{P(Z_i = z)}.
#'
#' @import stringr
#'
#' @export
simulate_data <- function(sample_size = 10000,
                          IVs_ps = list(c(0.5, 0.5), c(0.5, 0.5)),
                          indIVs_on_X = NULL, indIVs_on_Y = rep(0, length(IVs_ps)),
                          U = distributions3::Bernoulli(0.5),
                          X_intercept = 0.05, U_on_X = 0.1,
                          Y_intercept = 0.05, X_on_Y = 0.05, U_on_Y = 0.1,
                          rho = NULL, depIVs_on_X = NULL, dep_probs = c(0.5, 0.5), depIVs_on_Y = NULL){

  ## Check that probabilities of P(Z_i = z) sum to 1 for all i
  IVs_ps_sums <- purrr::map_dbl(IVs_ps, sum)

  if (any(IVs_ps_sums != 1))
    stop("Each vector of IVs_ps must sum to 1.")

  ## Check that length of IVs_ps is the same as length of indIVs_on_X
  if (length(IVs_ps) != length(indIVs_on_X))
    stop("Make sure the number of coefficients in indIVs_on_X is the same as number of distributions provided in IVs_ps.")

  ## Check that rho and depIVs_on_X are either both NULL, or both specified
  if (is.null(depIVs_on_X) != is.null(rho))
    stop("If rho is specified, so must depIVs_on_X be, and vice versa.")

  ## Check that if depIVs_on_Y is specified, it's the same length as depIVs_on_X
  if (!is.null(depIVs_on_Y) & (length(depIVs_on_X) != length(depIVs_on_Y)))
    stop("If depIVs_on_Y is specified, length must be same as length of depIVs_on_X.")
  ## Default IV effects on Y to 0
  if (is.null(indIVs_on_Y))
    indIVs_on_Y <- rep(0, length(indIVs_on_X))

  if (is.null(depIVs_on_Y))
    depIVs_on_Y <- rep(0, length(depIVs_on_X))

  ## Start tibble with just U. Include effects on X and Y for later use.
  final_data <- tibble::tibble(var = "U",
                               var_on_X = U_on_X,
                               var_on_Y = U_on_Y,
                               simulated = map(var, function(x){
                                 distributions3::random(U, n = sample_size)
                               }))

  ## If independent IVs are included...
  if(length(indIVs_on_X) > 0){
    ## ... create tibble with simulations of these. Include effects on X and Y for later use.
    out_indIVs <- tibble::tibble(var = paste0("Z", seq_along(indIVs_on_X)),
                                 var_on_X = indIVs_on_X,
                                 var_on_Y = indIVs_on_Y,
                                 simulated = map2(var_on_X, var, function(x,y){
                                   ## Simulated IVs by sampling from 1,2 (or 1,2,3) with replacement.
                                   ## Prob of each is given as argument IVs_ps. Subtract one to turn into
                                   ## 0,1 (or 0,1,2) variable.
                                   var_i <- as.numeric(stringr::str_extract_all(y, pattern = "[0-9]+"))
                                   sample(1:length(IVs_ps[[var_i]]), replace = TRUE, size = sample_size,
                                          prob = IVs_ps[[var_i]]) - 1
                                 }))

    ## Add to full data
    final_data <- bind_rows(
      final_data,
      out_indIVs
    )

  }

  ## If dependent IVs are included...
  if(length(depIVs_on_X) > 0){
    tmp <- rbivNorm(n = sample_size,
                    rho = rho, probs = dep_probs)

    colnames(tmp) <- paste0("Z", length(indIVs_on_X) + c(1,2))
    out_depIVs <- tibble::as_tibble(tmp) %>%
      gather(key = var, value = simulated) %>%
      nest(simulated = simulated) %>%
      mutate(var_on_X = depIVs_on_X,
             var_on_Y = depIVs_on_Y,
             simulated = map(simulated, function(x) pull(x, simulated)))

    final_data <- bind_rows(
      final_data,
      out_depIVs
    )
  }

  ## Some manipulations
  tmp_data <- final_data %>%
    ## Multiple vars by effects on X and Y, respectively
    mutate(weighted_for_X = map2(var_on_X, simulated,
                                 function(x,y){
                                   tibble(weights_for_X = x*y)
                                 }),
           weighted_for_Y = map2(var_on_Y, simulated,
                                 function(x,y){
                                   tibble(weights_for_Y = x*y)
                                 })) %>%
    ## Expand simulated and weighted observations
    unnest(cols = c(simulated, weighted_for_X, weighted_for_Y)) %>%
    ## Group by var, and create id column
    group_by(var) %>%
    mutate(i = row_number()) %>%
    ## Group by id, and get contribution from U and Z's on X and Y
    group_by(i) %>%
    mutate(pX = sum(weights_for_X) + X_intercept,
           pY = sum(weights_for_Y) + Y_intercept) %>%
    ungroup() %>%
    select(var, simulated, pX, pY, i)

  ## Simulate X
  out_data <- tmp_data %>%
    pivot_wider(id_cols = c(i, var, pX, pY), names_from = var, values_from = simulated) %>%
    mutate(X = map_dbl(pX, ~rbinom(n = 1, size = 1, prob = plogis(.x))),
           pY = pY + X*X_on_Y,
           Y = map_dbl(pY, ~rbinom(n = 1, size = 1, prob = plogis(.x)))) %>%
    select(-i)

  #message(paste("Number of pX > 1:", sum(out_data$pX > 1)))
  #message(paste("Number of pY > 1:", sum(out_data$pY > 1)))

  coefs <- tibble::tibble(
    coef = c(X_intercept,
             indIVs_on_X,
             depIVs_on_X,
             U_on_X,
             Y_intercept,
             X_on_Y,
             U_on_Y),
    effect = c("Xintercept",
               if(length(indIVs_on_X) > 0){
                 glue::glue("Z{1:length(indIVs_on_X)}_on_X")
               } else {
                 NULL
               },
               if(length(depIVs_on_X) > 0){
                 glue::glue("Z{length(indIVs_on_X) + c(1,2)}_on_X")
               } else {
                 NULL
               },
               "U_on_X",
               "Yintercept",
               "X_on_Y",
               "U_on_Y")
  )

  if(length(depIVs_on_X) > 0){
    coefs <- coefs %>%
      bind_rows(
        tibble::tibble(coef = rho,
                       effect = glue::glue("Z{length(indIVs_on_X)+1}_on_Z{length(indIVs_on_X)+2}"))
      )
  }

  if(any(indIVs_on_Y != 0)){
    coefs <- bind_rows(
      coefs,
      tibble::tibble(coef = indIVs_on_Y,
                     effect = glue::glue("Z{1:length(indIVs_on_Y)}_on_Y"))
    ) %>%
      filter(coef != 0)
  }

  if(any(depIVs_on_Y != 0)){
    coefs <- bind_rows(
      coefs,
      tibble::tibble(coef = depIVs_on_Y,
                     effect = glue::glue("Z{length(depIVs_on_Y) + c(1,2)}_on_Y"))
    ) %>%
      filter(coef != 0)
  }

  return(list(simulated_data = select(out_data, -pX, -pY),
              coefficients = coefs,
              IVs_ps = IVs_ps,
              dep_probs = dep_probs))
}
