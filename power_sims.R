library(tidyverse)
library(ACEBounds)
library(distributions3)
library(furrr)

plan(multicore)

for (i in 1:10){
many_sims <- expand_grid(i = i,
                         indIVs_on_X = 1:10/15*2,
                         X_on_Y = c(0.5, 1:10)) %>%
  mutate(
    sim_data = future_map2(
      indIVs_on_X, X_on_Y,
      ~simulate_data(
        sample_size = 10e5,
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

#many_sims_w_bounds <- many_sims %>%
#  mutate(sum_stats = list(probs_from_data(sim_data$simulated_data, X, Y, Z1, data_format = "bivariate")),
#         get_bounds_res = list(get_bounds(gammas = sum_stats$gammas, thetas = sum_stats$thetas, stop = FALSE)),
#         bounds = list(get_bounds_res$interval))

write_rds(many_sims, paste0("sims_for_power_", stringr::str_pad(i, width = 2, side = "left", pad = "0"),".Rds"))
}
# many_sims_w_bounds %>%
#   unnest_wider(bounds) %>%
#   arrange(desc(lower)) %>%
#   ggplot(aes(x = indIVs_on_X)) +
#     geom_errorbar(aes(ymin = lower, ymax = upper)) +
#     facet_grid(~X_on_Y)



all_sims <- tibble(sims_for_power = list.files(pattern = "sims_for_power")) %>%
  rowwise() %>%
  mutate(sims = list(read_rds(sims_for_power))) %>%
  ungroup() %>%
  unnest(sims)


all_sims_w_bounds <- all_sims %>%
  mutate(simulated_data = map(sim_data, "simulated_data"),
         sum_stats = map(simulated_data, ~probs_from_data(.x, X, Y, Z1, data_format = "bivariate")),
         get_bounds_res = map(sum_stats, ~get_bounds(gammas = .x$gammas, thetas = .x$thetas, stop = FALSE)),
         bounds = map(get_bounds_res, "interval"))


all_sims_w_bounds %>%
  unnest_wider(bounds) %>%
  mutate(zero = lower < 0 & upper > 0) %>%
  ggplot(aes(x = i, color = zero)) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    geom_errorbar(aes(ymin = lower, ymax = upper)) +
    facet_grid(X_on_Y ~ indIVs_on_X)


all_sims_w_bounds %>%
  unnest_wider(bounds) %>%
  group_by(X_on_Y, indIVs_on_X) %>%
  summarize(
    s_lower = sd(lower),
    s_upper = sd(upper),
    lower_range = max(lower) - min(lower),
    upper_range = max(upper) - min(upper)
  )
