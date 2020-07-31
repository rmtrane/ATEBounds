#' Wrapper of bpbounds::bpbounds
#'
#' Function that allows you to get bounds using the bpbounds::bpbounds
#' function from actual data.
#'
#' @param input_data data.frame with at least three columns containing actual data.
#' @param Y column to be used for output
#' @param X column to be used for exposure
#' @param Z columns to be used for IVs. If NULL, all columns with names matching Z[0-9]+ will be used.
#'
#' @export
bpbounds_from_data <- function(input_data,
                               Y = "Y", X = "X", Z = NULL, ...){

  if(is.null(Z)){
    Z <- str_extract_all(colnames(input_data), pattern = "Z[0-9]+") %>% unlist()
  }

  res <- tibble(IV = Z,
                bpbounds = map(IV, function(Z, ...){

                  tabs_formula <- as.formula(glue::glue("~ {X} + {Y} + {Z}"))
                  tabled_data <- xtabs(tabs_formula, data = input_data)
                  return(bpbounds(tabled_data, ...))
                }))
  return(res)
}

#' Tidy up output from bpbounds
#'
#' Return output from bpbounds::bpbounds in a tidy format.
#'
#' @param bpbs Output from bpbounds::bpbounds
#'
#' @export
tidy_bpbounds <- function(bpbs){
  bpbs_summary <- summary(bpbs)

  if(bpbs_summary$inequality){
    output <- bpbs_summary$bounds %>%
      mutate_if(is.factor, as.character) %>%
      mutate(Assumption = "IV Inequality",
             `Assumption Holds` = bpbs_summary$inequality)
  } else {
    output <- tibble(Assumption = "IV Inequality",
                     `Assumption Holds` = bpbs_summary$inequality,
                     `Causal parameter` = NA_character_,
                     `Lower bound` = NA_real_,
                     `Upper bound` = NA_real_)
  }

  if(bpbs_summary$monoinequality){
    output <- bpbs_summary$bounds %>%
      mutate_if(is.factor, as.character) %>%
      mutate(Assumption = "Monotonicity",
             `Assumption Holds` = bpbs_summary$monoinequality) %>%
      bind_rows(output)
  } else {
    output <- bind_rows(output,
                        tibble(Assumption = "Monotonicity",
                               `Assumption Holds` = bpbs_summary$monoinequality,
                               `Causal parameter` = NA_character_,
                               `Lower bound` = NA_real_,
                               `Upper bound` = NA_real_))
  }

  out <- output %>%
    dplyr::select(Assumption, `Assumption Holds`, everything()) %>%
    arrange(Assumption)

  return(as_tibble(out))
}

