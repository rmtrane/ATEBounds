library(tidyverse)

sim_zetas <- function(n_z_levels = 3){
  tmp_ps <- map(1:n_z_levels, ~runif(4))

  ps <- unlist(map(tmp_ps, ~.x / sum(.x)))

  ps_table <- array(ps,
                    dim = c(2,2,n_z_levels),
                    dimnames = list(x = 0:1, y = 0:1, z = 0:(n_z_levels - 1)))

  return(ps_table)
}

library(ACEBounds)
library(furrr)

plan(multiprocess)

many_tri_bounds <- expand_grid(i = 1:100000,
                               k = 3:4) %>%
  #sample_n(100) %>%
  mutate(ps_table = future_map2(i, k, ~sim_zetas(n_z_levels = .y)),
         bounds = future_map(ps_table, ~get_bounds(zetas = .x,
                                                   warning = FALSE,
                                                   stop = FALSE),
                             .progress = TRUE)) %>%
  rowwise() %>%
  mutate(violations = bounds$constraints_violated,
         intervals = list(bounds$interval),
         thetas = list(thetas_from_tabp(ps_table))) %>%
  ungroup() %>%
  unnest_wider(col = intervals) %>%
  mutate(strength = future_map_dbl(thetas, ~max(abs(outer(.x, .x, FUN = "-"))),
                                   .progress = TRUE),
         width = upper - lower)

write_rds(many_tri_bounds,
          path = here::here("data/many_tri_bounds.Rds"))

if(FALSE){
  many_tri_bounds <- read_rds(here::here("data/many_tri_bounds.Rds"))

  violation_summaries <- many_tri_bounds %>%
    count(k, violations, upper < lower)

  write_csv(violation_summaries,
            here::here("data/many_tri_bounds_violations.csv"))

  width_vs_strength <- many_tri_bounds %>%
    filter(!violations,
           upper >= lower) %>%
    ggplot(aes(x = strength, y = width)) +
    geom_abline(slope = -1, intercept = 1) +
    geom_point(size = 0.3, alpha = 0.1) +
    facet_grid(~k,
               labeller = function(...) label_both(..., sep = " = ")) +
    scale_x_continuous(limits = c(0, 1),
                       expand = expansion(mult = 0, add = 0)) +
    scale_y_continuous(limits = c(0, 1),
                       expand = expansion(mult = 0, add = 0)) +
    theme_bw()

  ggsave(width_vs_strength,
         filename = here::here("figures/trivariate_widths_vs_strengths.png"))
}
