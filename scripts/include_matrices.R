library(tidyverse)
matrices_from_polymake <- expand_grid(data_format = c("bivariate", "trivariate"),
            n_z_levels = 2:4,
            x_monotone = c(TRUE, FALSE),
            y_monotone = c(TRUE, FALSE)) %>%
  rowwise() %>%
  mutate(matrix = list(read_polymake_results(paste0("~/Documents/UW-Madison/BoundsFromMultipleIVs/", data_format, "_bound_matrices/",
                                                    "n_z_levels-", n_z_levels,"-x_monotone-", x_monotone, "-y_monotone-", y_monotone),
                                             data_format = data_format))) %>%
  ungroup()

usethis::use_data(matrices_from_polymake, internal = TRUE, overwrite = TRUE)