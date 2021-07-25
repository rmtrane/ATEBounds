library(tidyverse)

## Old code
# matrices_from_polymake <- expand_grid(data_format = c("bivariate", "trivariate"),
#             n_z_levels = 2:4,
#             x_monotone = c(TRUE, FALSE),
#             y_monotone = c(TRUE, FALSE)) %>%
#   rowwise() %>%
#   mutate(matrix = list(
#     read_polymake_results(here::here(paste0(data_format, "_bound_matrices/",
#                                             "n_z_levels-", n_z_levels,"-x_monotone-", x_monotone, "-y_monotone-", y_monotone)),
#                           data_format = data_format)
#   )
#   ) %>%
#   ungroup()


matrices_from_polymake <- tibble(data_format = c("bivariate", "trivariate"),
                                 files = map(data_format, ~list.files(here::here(paste0(.x, "_bound_matrices")),
                                                                      full.names = TRUE))) %>%
  unnest_longer(files) %>%
  mutate(
    matrix = map2(files, data_format, read_polymake_results),
    files = basename(files),
    params = map(files, ~str_split(.x, pattern = "n_z_levels-|-x_monotone-|-y_monotone-", simplify = TRUE) %>% .[,-1] %>% setNames(., c("n_z_levels", "x_monotone", "y_monotone")))
  ) %>%
  select(-files) %>%
  relocate(matrix, .after = params)  %>%
  unnest_wider(params) %>%
  mutate(
    n_z_levels = as.numeric(n_z_levels),
    x_monotone = as.logical(x_monotone),
    y_monotone = as.logical(y_monotone)
  ) %>%
  arrange(
    data_format, n_z_levels, !x_monotone, !y_monotone
  )

usethis::use_data(matrices_from_polymake, internal = TRUE, overwrite = TRUE)