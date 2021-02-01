library(tidyverse)
library(ACEBounds)
library(distributions3)
library(furrr)

set.seed(8072659)

sample_size <- 1e6

plan(multicore, workers = 3)

all_indIVs_on_X <- expand_grid(p = c(10, 50),
                               many_weak = c(TRUE, FALSE),
                               scenario = c("MR", "power")) %>%
  rowwise() %>%
  mutate(
    indIVs_on_X = case_when(many_weak & scenario == "MR" ~ list(c(seq(0.0, 0.01, length.out = p-1), 0.2)),
                            many_weak & scenario == "power" ~ list(c(seq(1, 1.2, length.out = p-1), 4)),
                            !many_weak & scenario == "MR" ~ list(seq(0, 0.2, length.out = p)),
                            !many_weak & scenario == "power" ~ list(seq(1, 4, length.out = p-1)),
                            TRUE ~ list(NA))
  ) %>%
  ungroup()

all_combinations <- expand_grid(
  p = c(10, 50),
  many_weak = c(TRUE, FALSE),
  scenario = c("MR", "power"),
  X_on_Y = c(0.25, 0.5, 1, 1.5, 2),
  U_on_XY = c(0.1, 0.5)
) %>%
  left_join(all_indIVs_on_X)

for (i in 1:nrow(all_combinations)) {

  cat(i, "of", nrow(all_combinations), "\n")

  tmp <- all_combinations[i, ] %>%
    mutate(
      sim_data = future_pmap(
        .l = list(indIVs_on_X, X_on_Y, U_on_XY),
        function(x, y, z) {
          simulate_data(
            sample_size = sample_size,
            IVs_ps = rep(list(c(0.25, 0.5, 0.25)), times = length(x)),
            X_intercept = -sum(x),
            Y_intercept = -y / 2,
            indIVs_on_X = x,
            indIVs_on_Y = 0,
            U = distributions3::Normal(),
            U_on_X = z,
            U_on_Y = z,
            X_on_Y = y
          )
        },
        .options = furrr_options(seed = TRUE)
      )
    )

  write_rds(
    tmp,
    here::here("data/multiple_IV_sims", paste0(str_pad(i, width = 2, side = "left", pad = 0), ".Rds"))
  )
}