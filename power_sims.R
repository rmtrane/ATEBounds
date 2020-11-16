library(tidyverse)
library(ACEBounds)
library(distributions3)
library(furrr)

plan(multiprocess)

many_sims <- expand_grid(indIVs_on_X = 1:10/15*2,
                         X_on_Y = c(0.5, 1:10)) %>%
  mutate(
    sim_data = future_map2(
      indIVs_on_X, X_on_Y,
      ~simulate_data(
        sample_size = 10e6,
        IVs_ps = list(c(0.25, 0.5, 0.25)),
        X_intercept = -1,
        indIVs_on_X = .x,
        indIVs_on_Y = 0,
        U = distributions3::Normal(),
        X_on_Y = .y
      ),
      .options = furrr_options(seed = TRUE)
    )
  )

many_sims_w_bounds <- many_sims %>%
  mutate(sum_stats = list(probs_from_data(sim_data$simulated_data, X, Y, Z1, data_format = "bivariate")),
         get_bounds_res = list(get_bounds(gammas = sum_stats$gammas, thetas = sum_stats$thetas, stop = FALSE)),
         bounds = list(get_bounds_res$interval))

# many_sims_w_bounds %>%
#   unnest_wider(bounds) %>%
#   arrange(desc(lower)) %>%
#   ggplot(aes(x = indIVs_on_X)) +
#     geom_errorbar(aes(ymin = lower, ymax = upper)) +
#     facet_grid(~X_on_Y)
