library(tidyverse)
library(ACEBounds)
library(furrr)

plan(multicore, workers = 3)

all_multiple_IV_sims <- list.files(here::here("data/multiple_IV_sims"), pattern = "[0-9]+.Rds", full.names = TRUE)

bounds_and_ATE <- tibble(data_file = all_multiple_IV_sims) %>%
  mutate(
    subset = map(
      data_file,
      function(x) {
        tmp <- read_rds(x)

        out <- tmp %>%
          mutate(
            ATE = map_dbl(sim_data, ATE_from_simulated_data),
            sums_and_bounds = map(
              sim_data,
              function(x) {
                sums_and_bounds <- x$simulated_data %>%
                  pivot_longer(starts_with("Z"), names_to = "IV", values_to = "z") %>%
                  nest_by(IV) %>%
                  ungroup() %>%
                  mutate(
                    sum_stats = future_map(data, ~probs_from_data(.x, X, Y, z, data_format = "bivariate")),
                    get_bounds_res = future_map(sum_stats, ~get_bounds(gammas = .x$gammas, thetas = .x$thetas, stop = FALSE)),
                    bounds = future_map(get_bounds_res, "interval")
                  ) %>%
                  select(-data)

                coefs <- x$coefficients %>%
                  filter(str_sub(effect, end = 1) == "Z") %>%
                  mutate(effect = str_remove(effect, "_on_X")) %>%
                  rename(IV = effect, indIVs_on_X = coef)

                left_join(
                  sums_and_bounds,
                  coefs,
                  by = "IV"
                )
              }
            )
          ) %>%
          select(-sim_data)

        return(out)
      }
    )
  )

write_rds(bounds_and_ATE, here::here("data/multiple_IV_sims/bounds_and_ATE.Rds"))


