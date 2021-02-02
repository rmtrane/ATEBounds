## ----message = FALSE, warning = FALSE-----------------------------------------
library(tidyverse)

bounds_and_ATE <- read_rds(here::here("data/power_bounds_and_ATE_combined.Rds"))

bounds_and_ATE_unnested <- bounds_and_ATE %>% 
  unnest_wider(sum_stats) %>%
  unnest_wider(bounds)

head(bounds_and_ATE_unnested)

## -----------------------------------------------------------------------------
all_combinations_of_coefs <- bounds_and_ATE %>%
  arrange(X_on_Y, indIVs_on_X, U_on_XY) %>%
  select(`$\\gamma_1$` = indIVs_on_X,
         `$\\beta_1$` = X_on_Y,
         `$\\gamma_U$` = U_on_XY) %>%
  pivot_longer(cols = everything()) %>%
  unique() %>%
  group_by(Coefficient = name) %>%
  summarize(Values = paste(value, collapse = ", ")) %>%
  pivot_wider(names_from = Coefficient, values_from = Values)

pander::pander(all_combinations_of_coefs)

## -----------------------------------------------------------------------------
kableExtra::kable(all_combinations_of_coefs, format = "latex", escape = FALSE, booktabs = TRUE) %>%
  kableExtra::column_spec(2, width = "1.5in") %>%
  write_file(file = here::here("tables", "sim_coefficients_table.tex"))

## ----width = 4, height = 4.7, dpi = 300---------------------------------------
(coefs_vs_strength <- bounds_and_ATE_unnested %>%
  mutate(strength = map_dbl(thetas, ~.x[3] - .x[1]),
         U_on_XY = as.character(U_on_XY)) %>%
  ggplot(aes(y = indIVs_on_X, x = strength,
             color = U_on_XY)) +
    geom_line() +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 6),
                       expand = expansion(mult = 0, add = c(0.1, 0.2))) +
    scale_color_manual(values = c("black", "red", "blue", "purple")) +
    labs(color = bquote(gamma[U]),
         x = "Strength of IV",
         y = bquote(gamma[Z])) +
    guides(color = guide_legend(nrow = 2)) +
    coord_fixed(ratio = 1/6) +
    theme_bw() +
    theme(legend.position = "top"))

ggsave(here::here("figures/MR_coefs_vs_strength.png"),
       coefs_vs_strength, dpi = 300,
       width = 4, height = 4.7)

## -----------------------------------------------------------------------------
(pretty_plot <- bounds_and_ATE_unnested %>%
  mutate(zero = if_else(lower < 0 & upper > 0,
                        "Overlaps Zero", "Does Not Overlap Zero"),
         X_on_Y = paste0("beta[X] ==", X_on_Y),
         U_on_XY = paste0("gamma[U] ==", U_on_XY)) %>%
  arrange(desc(lower)) %>%
  ggplot(aes(y = indIVs_on_X, color = zero)) +
    geom_smooth(se = FALSE, aes(x = upper, group = "loess"), color = "black",
                method = "loess",
                formula = y ~ x,
                size = 0.2) +
    geom_smooth(se = FALSE, aes(x = lower, group = "loess"), color = "black",
                method = "loess",
                formula = y ~ x,
                size = 0.2) +
    geom_vline(xintercept = 0) +
    geom_vline(aes(xintercept = ATE, color = "ATE")) +
    geom_errorbar(aes(xmin = lower, xmax = upper)) +
    facet_grid(U_on_XY ~ X_on_Y,
               labeller = label_parsed) +
    lims(y = c(0, 6),
         x = c(-1, 1)) +
    scale_color_manual(
      values = c("Overlaps Zero" = "black", "Does Not Overlap Zero" = "red", "ATE" = "blue")
    ) +
    labs(
      y = bquote(gamma[1]),
      x = "ATE",
      color = ""
    ) +
    theme_bw() +
    theme(legend.position = "top"))

ggsave(here::here("figures/power.png"),
       pretty_plot, dpi = 300,
       height = 6, width = 8, units = "in")

## -----------------------------------------------------------------------------
loess_model <- bounds_and_ATE_unnested %>%
  loess(data = .,
        lower ~ indIVs_on_X + X_on_Y + U_on_XY,
        span = 0.5)

## -----------------------------------------------------------------------------
roots <- bounds_and_ATE_unnested %>%
  select(X_on_Y, U_on_XY, ATE) %>%
  mutate(ATE = round(ATE, digits = 2)) %>% # to make sure ATE don't vary due to differences from randomness
  unique() %>%
  rowwise() %>%
  mutate(uniroot_res = list(uniroot(f = function(x) predict(loess_model,
                                                            newdata = data.frame(indIVs_on_X = x,
                                                                                 X_on_Y = X_on_Y, U_on_XY = U_on_XY)),
                                    interval = c(0.2, 6))),
         indIVs_on_X = uniroot_res$root) %>%
  ungroup()

roots %>% select(-uniroot_res) %>% rename("$\\gamma_1$" = indIVs_on_X) %>% pander::pander()

## -----------------------------------------------------------------------------
(loess_power <- roots %>%
  ggplot(aes(x = ATE, y = indIVs_on_X, group = U_on_XY, color = as.character(U_on_XY))) +
    geom_point(size = 0.5) +
    geom_line() +
    lims(
      y = c(0, 6),
      x = c(0, 1)
    ) +
    labs(
      color = bquote(gamma[U]),
      y = bquote(gamma[Z])
    ) +
    scale_color_manual(values = c("black", "red", "blue", "purple")) +
    guides(
      color = guide_legend(title.position = "top")
    ) +
    coord_fixed(ratio = 1/12) +
    theme_bw() +
    theme(legend.position = "top"))

ggsave(here::here("figures/loess_power.png"),
       loess_power, dpi = 300,
       height = 5, width = 8)

