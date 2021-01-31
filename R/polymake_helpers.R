#' Function to transform etas and deltas to gammas and thetas
#'
#' This is a helper function that allows us to transform extreme values of \eqn{P(Y = 1 | X = x, U)}, \eqn{P(X = 1 | Z = z, U)} into
#' values of \eqn{P(Y = y | Z = z, U)} and \eqn{P(X = x | Z = z, U)}, or \eqn{P(Y = y, X = x | Z = z, U)}. It is used to programmatically
#' generate vertices needed in the [polymake](https://www.polymake.org) program.
#'
#' @param etas vector of length 2 giving values of P(Y = 1 | X = 0, U), P(Y = 1 | X = 1, U)
#' @param deltas vector of same length as number of categories for IV giving P(X = 1 | Z = z, U), z = 0, ..., k, k = number of IV categories minus 1.
#' @param include_alpha logical; whether to include column for alpha or not.
#' @param data_format one of 'bivariate' or 'trivariate'. If bivariate, output will give extremes in terms of P(Y|Z) and P(X|Z).
#'   If trivariate, output will give extremes in therms of zeta{yxz} = P(Y = y, X = x | Z = z).
#'
#' @export
perform_transformation <- function(etas, deltas, data_format = c("bivariate", "trivariate"), include_alpha = FALSE){
  if(!data_format %in% c("bivariate", "trivariate") | length(data_format) != 1)
    stop("data_format must be one of 'bivariate' or 'trivariate'")

  if(data_format == "bivariate"){
    gammas <- setNames(c((1-etas[1])*(1-deltas) + (1-etas[2])*deltas,
                         etas[1]*(1-deltas) + etas[2]*deltas),
                       c(paste0("gamma", rep(0:1, each = length(deltas)), rep(1:length(deltas), 2))))

    thetas <- setNames(c(1-deltas, deltas),
                       c(paste0("theta", rep(0:1, each = length(deltas)), rep(1:length(deltas), 2))))

    output <- c(gammas, thetas)

  }

  if(data_format == "trivariate"){
    output <- setNames(c((1-etas[1])*(1-deltas),
                         (1-etas[2])*deltas,
                         etas[1]*(1-deltas),
                         etas[2]*deltas),
                       paste0("zeta", rep(c("00", "01", "10", "11"), each = length(deltas)), rep(1:length(deltas), 4)))
  }

  if(include_alpha){
    alpha <- setNames(etas[2] - etas[1], "alpha")
    output <- c(output, alpha)
  }

  return(output)
}

#' Create vertices for polymake program
#'
#' Function to generate vertices that can be fed to the [polymake](https://www.polymake.org) program for a number of different scenarios:
#'
#' * IV with a user specified number of categories
#' * \eqn{P(X = 1 | Z = z)} monotonically increasing
#' * \eqn{P(Y = 1 | Z = z)} monotonically increasing
#' * include column for \eqn{\alpha} (the ACE)
#' * output in terms of bivariate probabilities (\eqn{P(Y = y | Z = z)} and \eqn{P(X = x | Z = z)}), or trivariate probabilities (\eqn{P(Y = y, X = x | Z = z)}).
#'
#' @param n_z_levels number of categories for IV
#' @param x_monotone logical; is P(X = 1 | Z = z) monotonically increasing?
#' @param y_monotone logical; is P(Y = 1 | Z = z) monotonically increasing?
#' @param include_alpha logical; should a column for alpha = P(Y = 1 | X = 1, U) - P(Y = 1 | X = 0, U) be included?
#' @param data_format one of 'bivariate' or 'trivariate'. If bivariate, output will give extremes in terms of P(Y|Z) and P(X|Z).
#'   If trivariate, output will give extremes in therms of P(Y,X|Z).
#'
#' @return matrix of the extreme vertices in the specified scenario.
#'
#' @import purrr dplyr tidyr
#' @importFrom magrittr %>%
#'
#' @export
create_vertices <- function(n_z_levels, data_format = c("bivariate", "trivariate"),
                            x_monotone = FALSE, y_monotone = FALSE, include_alpha = TRUE){

  if(!data_format %in% c("bivariate", "trivariate") | length(data_format) != 1)
    stop("data_format must be one of 'bivariate' or 'trivariate'")

  deltas <- map(0:(n_z_levels-1), ~0:1) %>% setNames(nm = paste0("delta", 0:(n_z_levels-1)))
  etas <- list(eta0 = 0:1, eta1 = 0:1)

  vertices <- cross(c(etas, deltas)) %>%
    map_dfr(~tibble(value = unlist(.), var = names(.)) %>% pivot_wider(names_from = var, values_from = value)) %>%
    rowwise() %>%
    transmute(transformed = list(perform_transformation(etas = c_across(contains("eta")),
                                                        deltas = c_across(contains("delta")),
                                                        data_format = data_format,
                                                        include_alpha = include_alpha))) %>%
    unnest_wider(transformed)

  if(y_monotone){
    for(z in 1:(n_z_levels-1)){
      if(data_format == "bivariate"){
        z1 <- paste0("gamma1", z)
        z2 <- paste0("gamma1", z+1)

        vertices <- vertices[vertices[[z1]] <= vertices[[z2]],]
      }

      if(data_format == "trivariate"){
        z1 <- paste0("zeta", "1", 0:1, z)
        z2 <- paste0("zeta", "1", 0:1, z+1)

        vertices <- vertices[rowSums(vertices[,z1]) <= rowSums(vertices[,z2]),]
      }
    }
  }

  if(x_monotone){
    for(z in 1:(n_z_levels-1)){
      if(data_format == "bivariate"){
        z1 <- paste0("theta1", z)
        z2 <- paste0("theta1", z+1)

        vertices <- vertices[vertices[[z1]] <= vertices[[z2]],]
      }

      if(data_format == "trivariate"){
        z1 <- paste0("zeta", 0:1, "1", z)
        z2 <- paste0("zeta", 0:1, "1", z+1)

        vertices <- vertices[rowSums(vertices[,z1]) <= rowSums(vertices[,z2]),]
      }
    }
  }

  return(vertices)
}

#' Create Polymake Script
#'
#' Function that writes script to run through polymake.
#'
#' @param vertices output from `create_vertices`
#' @param output_file_name name to be used for polymake output
#' @param polymake_script_name name to be used for polymake script
#' @param output_folder folder where you want polymake script and output to be saved
#' @param overwrite logical; if true, existing files will be overwritten
#'
#' @return Nothing is returned, but a file is written that can be fed to polymake.
#'
#' @export
create_script <- function(vertices, output_file_name,
                          polymake_script_name = "polymake_script",
                          overwrite = FALSE,
                          output_folder = "."){

  error <- 0

  if(!dir.exists(output_folder)){
    warning(paste("The folder", output_folder, "does not exist."))
    error <- 1
  }


  if(file.exists(paste(output_folder, polymake_script_name, sep = "/")) & !overwrite){
    warning(paste("The file", paste(output_folder, polymake_script_name, sep = "/"), "already exists. If you wish to overwrite it, specify 'overwrite = TRUE'"))
    error <- 1
  }

  if(file.exists(paste(output_folder, output_file_name, sep = "/")) & !overwrite){
    warning(paste("The file", paste(output_folder, output_file_name, sep = "/"), "already exists. If you wish to overwrite it, specify 'overwrite = TRUE'"))
    error <- 1
  }

  if(error < 1){
    vertices_for_script <- vertices %>%
      rowwise() %>%
      transmute(for_output = paste0("[", paste(c_across(everything()), collapse = ","), "]")) %>%
      unique() %>%
      pull(for_output) %>%
      paste0(., c(rep(",", length(.)-1), "],"))

    for_script <- c("use application 'polytope';",
                    'declare $tmp;',
                    paste0('$tmp = new Polytope<Rational>(VERTICES=>[', vertices_for_script[1]),
                    vertices_for_script[-1],
                    'LINEALITY_SPACE=>[]);',
                    paste0("open(my $f1, '> ", paste(output_folder, output_file_name, sep = "/"), "'); print $f1 $tmp->FACETS; close($f1);"))

    readr::write_lines(x = for_script,
                       file = paste(output_folder,
                                    polymake_script_name, sep = "/"))
  }

  return(list(output_file = paste(output_folder, output_file_name, sep = "/"),
              polymake_script = paste(output_folder, polymake_script_name, sep = "/"),
              error = error))
}

#' Read output from polymake
#'
#' This function reads in the output generated from running the polymake script created by `create_script`.
#'
#' @param output_file path to output file from polymake script. `output_file` output from `create_script`
#' @param data_format one of 'bivaraite' or 'trivariate'. If bivariate, output will give extremes in terms of P(Y|Z) and P(X|Z).
#'   If trivariate, output will give extremes in terms of P(Y,X|Z).
#'
#' @return tibble with matrix of constraints for bounds. I.e. \eqn{A \cdot p} must be non-negative for \eqn{p} to satisfy the IV model.
#'
#' @export
read_polymake_results <- function(output_file, data_format = c("bivariate", "trivariate")){
  if(!data_format %in% c("bivariate", "trivariate") | length(data_format) != 1)
    stop("data_format must be one of 'bivariate' or 'trivariate'")

  # Read first row, and infer n_z_levels
  first_row <- readr::read_delim(output_file, delim = " ", col_names = FALSE, n_max = 1)
  n_z_levels <- (ncol(first_row) - 1)/4


  if(data_format == "bivariate")
    columns <- c(paste0("gamma", rep(0:1, each = n_z_levels), rep(1:n_z_levels, 2)),
                 paste0("theta", rep(0:1, each = n_z_levels), rep(1:n_z_levels, 2)),
                 "alpha")

  if(data_format == "trivariate")
    columns <- c(paste0("zeta", rep(c("00", "01", "10", "11"), each = n_z_levels), rep(1:n_z_levels, 4)),
                 "alpha")

  # Read file
  from_polymake <- readr::read_delim(file = output_file, delim = " ",
                                     col_names = columns)

  return(from_polymake)
}
