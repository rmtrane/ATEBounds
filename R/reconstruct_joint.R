#' P(X,Y|Z) from P(X|Z), P(Y|Z), and Cov(X,Y|Z)
#'
#' For specified values of \eqn{P(X = 1 | Z = z)}, \eqn{P(Y = 1 | Z = z)}, and \eqn{\text{Cov}(X, Y | Z = z)}, \eqn{z = 0,...,k-1},
#' the joint \eqn{P(Y = y, X = x | Z = z)} is reconstructed.
#'
#' @param gammas vector of length `n_z_levels` giving P(Y = 1 | Z = z), z = 0, ..., (`n_z_levels` - 1).
#' @param thetas vector of length `n_z_levels` giving P(X = 1 | Z = z), z = 0, ..., (`n_z_levels` - 1).
#' @param cov vector of length `n_z_levels` giving Cov(X, Y | Z = z), z = 0,..., (`n_z_levels` - 1).
#' @param stop logical; if `TRUE` (default), an error is thrown if the resulting probability distribution
#'   is invalid, or does not satisfy the IV inequalities. Otherwise, a warning is printed, and all
#'   probabilities returned set to NA.
#' @param message logical; only used when `stop = TRUE`. If false, no warning message is printed when the
#'   resulting probability distribution is invalid. (Default: `TRUE`)
#' @param from_polymake matrix created by the polymake program with the constraints a valid joint distribution
#'   should meet. If `NULL` (default), a previously prepared matrix is used if data seems to match one of these.
#' @param mono_x logical; should \eqn{P(X = 1 | Z = z)} be monotonically increasing?
#' @param mono_y logical; should \eqn{P(Y = 1 | Z = z)} be monotonically increasing?
#' @param return_bounds logical; should the resulting bounds be included in the output?
#'
#' @return A \eqn{2x2xk} array with \eqn{P(Y = y, X = x | Z = z)} as entry (y+1, x+1, z+1).
#'
#' @export
joint_from_marginal <- function(thetas, gammas, cov, stop = TRUE, message = TRUE,
                                x_mono = FALSE, y_mono = FALSE, from_polymake = NULL,
                                return_bounds = FALSE){
  if(!(length(thetas) == length(gammas) & length(gammas) == length(cov)))
    stop("All three vectors must be of equal length.")

  tabp <- array(data = NA,
                dim = c(2,2,length(thetas)),
                dimnames = list(x = 0:1, y = 0:1, z = 0:(length(thetas)-1)))

  tabp[x="1",y="1",] <- purrr::map_dbl(seq_along(thetas), ~thetas[.x]*gammas[.x] + cov[.x])
  tabp[x="1",y="0",] <- purrr::map_dbl(seq_along(thetas), ~thetas[.x]*(1-gammas[.x]) - cov[.x])
  tabp[x="0",y="1",] <- purrr::map_dbl(seq_along(thetas), ~(1-thetas[.x])*gammas[.x] - cov[.x])
  tabp[x="0",y="0",] <- purrr::map_dbl(seq_along(thetas), ~(1-thetas[.x])*(1-gammas[.x]) + cov[.x])

  if(any(tabp < 0 | tabp > 1)){
    if (stop)
      stop("Specified values result in probabilities below 0 or above 1.")

    if (message)
      message("Specified values result in probabilities below 0 or above 1.")

    tabp <- array(data = NA, dim = dim(tabp), dimnames = dimnames(tabp))

  }

  if(max(apply(apply(tabp, c(1,2), max), 1, sum)) > 1){
    if (stop)
      stop("Specified values do not satisfy the IV model.")

    if (message)
      message("Specified values do not satisfy the IV model.")

    tabp <- array(data = NA, dim = dim(tabp), dimnames = dimnames(tabp))

    return(tabp)
  }

  bivariate_bounds <- get_bounds(gammas = gammas, thetas = thetas, stop = FALSE, warning = FALSE,
                                 x_mono = x_mono, y_mono = y_mono, from_polymake = from_polymake)

  trivariate_bounds <- get_bounds(zetas = tabp, stop = FALSE, warning = FALSE,
                                  x_mono = x_mono, y_mono = y_mono, from_polymake = from_polymake)

  if(bivariate_bounds$constraints_violated | bivariate_bounds$width < 0){
    if (stop)
      stop("thetas and gammas do not satisfy the constraints, or upper bound less than lower bound.")

    if (message)
      message("thetas and gammas do not satisfy the constraints, or upper bound less than lower bound.")

    tabp <- array(data = NA, dim = dim(tabp), dimnames = dimnames(tabp))
  }

  if(trivariate_bounds$constraints_violated | bivariate_bounds$width < 0){
    if (stop)
      stop("The reconstructed P(Y,X|Z) do not satisfy the constraints, or upper bound less than lower bound.")

    if (message)
      message("The reconstructed P(Y,X|Z) do not satisfy the constraints, or upper bound less than lower bound.")

    tabp <- array(data = NA, dim = dim(tabp), dimnames = dimnames(tabp))
  }

  if(return_bounds)
    return(tibble(joint = list(tabp),
                  bounds = list(trivariate_bounds$interval)))

  return(tabp)
}

#' Valid range of Covs
#'
#' From values of P(X = 1 | Z = z) and P(Y = 1 | Z = z) (for all z = 0, ..., `n_z_levels` - 1),
#' this function returns ranges of valid values for the quantities Cov(X,Y|Z = z) and Cov(X,Y|Z = z1) - Cov(X,Y|Z = z2),
#' such that there probability distribution defined by the input and any set of covariances chosen to fulfill the
#' constraints will be valid, and not violated the IV inequalities.
#'
#' @param gammas vector of length `n_z_levels` giving P(Y = 1 | Z = z), z = 0, ..., (`n_z_levels` - 1).
#' @param thetas vector of length `n_z_levels` giving P(X = 1 | Z = z), z = 0, ..., (`n_z_levels` - 1).
#' @param x_mono logical; should monotonicity of P(X = 1 | Z = z) be assumed?
#'
#' @return A list with four elements:
#'
#' * `thetas` and `gammas` provided
#' * `potential_cov_range` is a tibble with three columns: `min` and `max` give bounds on Cov(X,Y|Z=z) for the
#'   value of `z`
#' * `constraints`is a tibble with four columns: `min` and `max` give bounds on Cov(X,Y|Z = z1) - Cov(X,Y | Z = z2) for the
#'   corresponding values of `z1` and `z2`.
#'
#' @export
potential_covs <- function(thetas, gammas, x_mono = FALSE){

  ##### Limits obtained from 0 <= P(X,Y|Z) <= 1.
  ## Upper limits
  max_11 <- purrr::map_dbl(seq_along(thetas), ~1-thetas[.x]*gammas[.x])
  max_10 <- purrr::map_dbl(seq_along(thetas), ~thetas[.x]*(1-gammas[.x]))
  max_01 <- purrr::map_dbl(seq_along(thetas), ~(1-thetas[.x])*gammas[.x])
  max_00 <- purrr::map_dbl(seq_along(thetas), ~1-(1-thetas[.x])*(1-gammas[.x]))

  ## Lower limits
  min_11 <- purrr::map_dbl(seq_along(thetas), ~0-thetas[.x]*gammas[.x])
  min_10 <- purrr::map_dbl(seq_along(thetas), ~thetas[.x]*(1-gammas[.x])-1)
  min_01 <- purrr::map_dbl(seq_along(thetas), ~(1-thetas[.x])*gammas[.x]-1)
  min_00 <- purrr::map_dbl(seq_along(thetas), ~0-(1-thetas[.x])*(1-gammas[.x]))

  ## Combine
  maxs <- dplyr::bind_cols(`11` = max_11, `10` = max_10, `01` = max_01, `00` = max_00)
  mins <- dplyr::bind_cols(`11` = min_11, `10` = min_10, `01` = min_01, `00` = min_00)

  ## tibble to return
  potential_cov_range <- bind_cols(z = 1:length(thetas) - 1,
                                   min = apply(mins, 1, max), max = apply(maxs, 1, min))

  ## From the IV inequalities, we can find constraints on the differences cov(X, Y | Z = z1) - cov(X , Y | Z = z2).
  constraints <- expand_grid(z1 = seq_along(thetas), z2 = seq_along(thetas)) %>%
    filter(z1 < z2) %>%
    rowwise() %>%
    # mutate(const = list(tibble(max = min(c(thetas[z1]*(1-gammas[z1]) + thetas[z2]*gammas[z2],
    #                                        1 - (1-thetas[z1])*(1-gammas[z1]) - (1-thetas[z2])*gammas[z2])),
    #                            min = max(c(thetas[z1]*(1-gammas[z1]) + thetas[z2]*gammas[z2],
    #                                        1 - (1-thetas[z1])*(1-gammas[z1]) - (1-thetas[z2])*gammas[z2]) - 1)))
    #        )  %>%
    mutate(const = list(tibble(min = max(-(1-thetas[z1])*(1-gammas[z1]) - (1-thetas[z2])*gammas[z2],
                                         thetas[z1]*(1-gammas[z1]) + thetas[z2]*gammas[z2] - 1,
                                         (1-thetas[z2])*(1-gammas[z2]) + (1-thetas[z1])*gammas[z1] - 1,
                                         -thetas[z2]*(1-gammas[z2]) - thetas[z1]*gammas[z1]),
                               max = min(1 - (1-thetas[z1])*(1-gammas[z1]) - (1-thetas[z2])*gammas[z2],
                                         thetas[z1]*(1-gammas[z1]) + thetas[z2]*gammas[z2],
                                         (1-thetas[z2])*(1-gammas[z2]) + (1-thetas[z1])*gammas[z1],
                                         1 - thetas[z2]*(1-gammas[z2]) - thetas[z1]*gammas[z1])))) %>%
    unnest(const) %>%
    ungroup()

  ## If monotonicity is assumed, we get an extra set of constraints
  if(x_mono){
    constraints <- constraints %>%
      rowwise() %>%
      mutate(mono_min = max(c(thetas[z1]*(1-gammas[z1]) - thetas[z2]*(1-gammas[z2]),
                              (1-thetas[z2])*(1-gammas[z2]) - (1-thetas[z1])*(1-gammas[z1]))),
             mono_max = min(c(thetas[z2]*gammas[z2] - thetas[z1]*gammas[z1],
                              (1-thetas[z1])*gammas[z1] - (1-thetas[z2])*gammas[z2])),
             new_min = max(min, mono_min),
             new_max = min(max, mono_max)) %>%
      ungroup() %>%
      select(starts_with("z"), min = new_min, max = new_max)
  }

  constraints <- constraints %>%
    mutate(
      across(
        .cols = starts_with("z", ignore.case = FALSE),
        ~.x - 1
      )
    )


  return(list(thetas = thetas,
              gammas = gammas,
              potential_cov_range = potential_cov_range,
              constraints = constraints))
}

#' Plot valid covs
#'
#' Visualize the area of valid covariances.
#'
#' @param pot_covs output from `potential_covs()`
#'
#' @export
plot_valid_covs <- function(pot_covs, Zs = c(0,1)){
  if(nrow(pot_covs$potential_cov_range) > 2 & length(Zs) != 2)
    stop("Can only plot when Z is binary.")

  future_proof <- pot_covs$constraints %>%
    filter(z1 == z[1]) %>%
    rename(z = z2,
           diff_min = min,
           diff_max = max) %>%
    left_join(pot_covs$potential_cov_range, by = "z") %>%
    rename(z_min = min,
           z_max = max) %>%
    rowwise() %>%
    mutate(z1_min = min(z_min + c(diff_min, diff_max)),
           z1_max = max(z_max + c(diff_min, diff_max)))

  CovZ0_min <- max(filter(pot_covs$potential_cov_range, z == Zs[1])$min,
                   future_proof$z1_min)
  CovZ0_max <- min(filter(pot_covs$potential_cov_range, z == z[1])$max,
                   future_proof$z1_max)

  CovZ0_seq <- seq(CovZ0_min, CovZ0_max, length.out = 100)

  future_proof <- pot_covs$constraints %>%
    filter(z1 == Zs[2]) %>%
    rename(z = z2,
           diff_min = min,
           diff_max = max) %>%
    left_join(pot_covs$potential_cov_range, by = "z") %>%
    mutate(min_min = min + diff_min,
           max_max = max + diff_max)

  check_past <- pot_covs$constraints %>%
    filter(z2 == Zs[2])

  CovZ1_mins <- map_dbl(CovZ0_seq, ~max(.x - check_past$max,
                                        filter(pot_covs$potential_cov_range, z == Zs[2])$min,
                                        future_proof$min_min))
  CovZ1_maxs <- map_dbl(CovZ0_seq, ~min(.x - check_past$min,
                                        filter(pot_covs$potential_cov_range, z == Zs[2])$max,
                                        future_proof$max_max))

  ggplot(data = data.frame(x = CovZ0_seq, ymin = CovZ1_mins, ymax = CovZ1_maxs),
         aes(x = x)) +
    geom_hline(yintercept = filter(pot_covs$potential_cov_range, z == Zs[2]) %>% select(min, max) %>% unlist,
               color = "red", linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = filter(pot_covs$potential_cov_range, z == Zs[1]) %>% select(min, max) %>% unlist,
               color = "red", linetype = "dashed", alpha = 0.5) +
    geom_ribbon(aes(xmin = min(x), xmax = max(x),
                    ymin = ymin, ymax = ymax),
                fill = NA, color = "black",
                outline.type = "full") +
    # geom_line(aes(y = ymin)) +
    # geom_line(aes(y = ymax)) +
    lims(x = c(-1,1), y = c(-1,1)) +
    labs(x = glue::glue("Cov(X,Y|Z={Zs[1]})"), y = glue::glue("Cov(X,Y|Z={Zs[2]})"))
}

#' Sample valid distributions P(X,Y|Z)
#'
#' Returns a `n x 4` tibble with one id column, and one column for each level
#' of Z containing a random draw of a potential value of Cov(X,Y|Z=z). These values are drawn such that
#' the joint conditional distribution P(X,Y|Z) that can be constructed from these and the values of
#' P(X = 1 | Z) and P(Y = 1 | Z) that were provided to `potential_covs` satisfy the IV inequalities.
#'
#' @param pot_covs output from `potential_covs()`
#' @param n number of samples to draw
#' @param max_rejections number of rejections at which point we move on
#' @param x_mono should it be assumed that P(X = 1 | Z = z) is monotonically increasing?
#' @param y_mono should it be assumed that P(Y = 1 | Z = z) is monotonically increasing?
#' @param from_polymake if NULL (default), we look through the pre-constructed matrices to find one that matches our setup. Otherwise, a matrix can be provided.
#' @param return_bounds logical; should the bounds be returned for each sampled joint distribution? (Default: FALSE)
#' @param print_progress logical; should process be printed?
#' @param print_as_progress if not NULL (default), whatever is provided will be printed in beginning. Useful for batch runs.
#'
#' @return `n x 4` tibble with the following columns:
#'
#' * `id`: simply ids the run, i.e. is `1:n`
#' * `covs`: list of vectors where each vector gives the sampled Cov(X,Y|Z = z) for each z.
#' * `joint`: list of `2 x 2 x k` arrays giving the joint conditional distribution reconstructed using
#'   the drawn covariances
#' * `n_rejected`: number of draws that were rejected due to the reconstructed conditional distribution
#'   not satisfying the constraints.
#'
#' @export
sample_joint_probs <- function(pot_covs,
                               n = 1000,
                               max_rejections = Inf,
                               x_mono = FALSE,
                               y_mono = FALSE,
                               from_polymake = NULL,
                               return_bounds = FALSE,
                               print_progress = TRUE,
                               print_as_progress = NULL){

  if(!is.null(print_as_progress))
    cat("j =", print_as_progress, "\n")

  # Simple vector with values for Z
  z <- pot_covs$potential_cov_range$z

  # List with entry for each value of Z which gives constraints on
  # Cov(X,Y|Z=z) in terms of possible values of Cov(X,Y|Z=z1), where
  # z1 > z.
  future_proof <- map(
    z[-length(z)],
    ~pot_covs$constraints %>%
      filter(z1 == .x) %>%
      rename(z = z2,
             diff_min = min,
             diff_max = max) %>%
      left_join(pot_covs$potential_cov_range, by = "z") %>%
      rename(z_min = min,
             z_max = max) %>%
      mutate(z1_min = z_min + diff_min,
             z1_max = z_max + diff_max)
  )

  # List with entry for each value of Z which gives constraints on
  # Cov(X,Y|Z=z) in terms of possible values of Cov(X,Y|Z=z1), where
  # z1 < z.
  check_past <- map(z, ~pot_covs$constraints %>%
                      filter(z2 == .x))

  # Put in placeholder entry for last value of z.
  future_proof[[length(z)]] <- tibble(z1_min = NA, z1_max = NA)

  # Get initial candidates for Cov(X,Y | Z = 0).
  z0 <- runif(n = n,
              min = max(filter(pot_covs$potential_cov_range, z == z[1])$min,
                        future_proof[[1]]$z1_min),
              max = min(filter(pot_covs$potential_cov_range, z == z[1])$max,
                        future_proof[[1]]$z1_max))


  output <- tibble(id = 1:n) %>%
    mutate(tmp = map(id, function(i){

      if(i %% 100 == 0 & print_progress) cat("i = ", i, "\n")

      ## Reset number of rejections.
      n_rejected <- 0
      suggested_covs <- c(z0[i], rep(NA, length(z)-1))

      # As long as not all values of Cov(X,Y|Z=z) has been selected,
      # and we haven't yet hit the max number of rejections...
      while(any(is.na(suggested_covs)) & n_rejected < max_rejections){
        # ... starting at the second value of z
        j <- 2
        # ... and while we haven't gone through all values of z or hit max number of rejections...
        while (j <= length(z) & n_rejected < max_rejections){
          ## ... find the smallest and largest values of Cov(X,Y | Z = z[j])
          cur_min <- max(suggested_covs[!is.na(suggested_covs)] - check_past[[j]]$max,
                         filter(pot_covs$potential_cov_range, z == z[j])$min,
                         future_proof[[j]]$z1_min,
                         na.rm = TRUE)
          cur_max <- min(suggested_covs[!is.na(suggested_covs)] - check_past[[j]]$min,
                         filter(pot_covs$potential_cov_range, z == z[j])$max,
                         future_proof[[j]]$z1_max,
                         na.rm = TRUE)

          ## If smallest value is greater than largest value, we reject the
          ## chosen value of Cov(X,Y | Z = z[j-1]) by taking a step back. This
          ## means we redraw the value of Cov(X,Y | Z = z[j-1]).
          if(cur_min > cur_max | is.na(suggested_covs[j-1])){
            j <- j - 1
            n_rejected <- n_rejected + 1
          }

          suggested_covs[j] <- runif(n = 1,
                                     min = max(suggested_covs[!is.na(suggested_covs)] - check_past[[j]]$max,
                                               filter(pot_covs$potential_cov_range, z == z[j])$min,
                                               future_proof[[j]]$z1_min,
                                               na.rm = TRUE),
                                     max = min(suggested_covs[!is.na(suggested_covs)] - check_past[[j]]$min,
                                               filter(pot_covs$potential_cov_range, z == z[j])$max,
                                               future_proof[[j]]$z1_max,
                                               na.rm = TRUE))

          # Proceed to next j
          j <- j + 1
        }

        # If we got values of Cov(X,Y | Z=z) for all z's before hitting max_rejections,
        # we find the joint distribution.
        if (all(!is.na(suggested_covs))){
          joint <- joint_from_marginal(thetas = pot_covs$thetas, gammas = pot_covs$gammas, cov = suggested_covs,
                                       stop = FALSE, message = FALSE, x_mono = x_mono, y_mono = y_mono, return_bounds = return_bounds,
                                       from_polymake = from_polymake)

          # If the joint is not valid, reset suggested values of Cov(X,Y|Z=z).
          if((is.array(joint) && any(is.na(joint))) | (!is.array(joint) && any(is.na(joint$joint[[1]])))){
            suggested_covs <- rep(NA, length(z))
            n_rejected <- n_rejected + 1
          } else {
            out <- tibble(covs = list(setNames(suggested_covs, paste0("z", 0:(length(z)-1)))),
                          joint = list(joint),
                          n_rejected = n_rejected)
          }
        } else {
          suggested_covs <- rep(NA, length(z))
          n_rejected <- n_rejected + 1
        }
      }

      # If we hit max_rejections, return missing values.
      if (n_rejected >= max_rejections){
        out <- tibble(covs = list(NA), joint = list(NA), n_rejected = n_rejected); message("max_rejections hit for!")
      }
      return(out)
    })
    )

  return(unnest(output, tmp))
}
