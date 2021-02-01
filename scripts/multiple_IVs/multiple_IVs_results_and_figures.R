library(tidyverse)
library(ACEBounds)
library(furrr)

plan(multisession, workers = 3)

ATE_from_simulated_data <- function(from_simulate_data) {
  intercept <- filter(from_simulate_data$coefficients, effect == "Yintercept")$coef
  x_beta <- filter(from_simulate_data$coefficients, effect == "X_on_Y")$coef
  u_beta <- filter(from_simulate_data$coefficients, effect == "U_on_Y")$coef

  pY1X0 <- 1 / (1 + exp(-intercept - u_beta * from_simulate_data$simulated_data$U))
  pY1X1 <- 1 / (1 + exp(-intercept - x_beta - u_beta * from_simulate_data$simulated_data$U))

  return(mean(pY1X1 - pY1X0))
}


if(!file.exists(here::here("data/multiple_IV_sims/bounds_and_ATE.Rds"))){
  message("data/multiple_IV_sims/bounds_and_ATE.Rds doesn't exist. Run the script scripts/multiple_IVs/multiple_IVs_prep.R to create it.")
} else {

  bounds_and_ATE <- read_rds(here::here("data/multiple_IV_sims/bounds_and_ATE.Rds"))

  bounds_and_ATE_unnested <- bounds_and_ATE %>%
    unnest(subset) %>%
    select(-data_file, -indIVs_on_X) %>%
    unnest(sums_and_bounds)

  dilution_effect <- bounds_and_ATE_unnested %>%
    unnest_wider(bounds) %>%
    mutate(strength = map_dbl(sum_stats, ~max(abs(outer(.x$thetas, .x$thetas, "-"))))) %>%
    group_by(p, scenario, X_on_Y, U_on_XY, many_weak) %>%
    filter(indIVs_on_X == max(indIVs_on_X),
           X_on_Y == 1) %>%
    mutate(
      Scenario = case_when(
        !many_weak & scenario == "MR" ~ "Scenario 1",
        !many_weak & scenario == "power" ~ "Scenario 2",
        many_weak & scenario == "MR" ~ "Scenario 3",
        many_weak & scenario == "power" ~ "Scenario 4",
        TRUE ~ NA_character_
      ),
      indIVs_on_X = case_when(round(indIVs_on_X, digits = 3) == 0.162 ~ 0.2,
                              round(indIVs_on_X, digits = 2) %in% c(3.93, 3.67) ~ 4,
                              TRUE ~ round(indIVs_on_X, digits = 1)),
      p = paste("p", p, sep = " = ")
    ) %>%
    ungroup() %>%
    select(p = p, U_on_XY, indIVs_on_X, strength, Scenario) %>%
    pivot_wider(names_from = p, values_from = strength) %>%
    arrange(Scenario)

  write_csv(dilution_effect, file = here::here("tables/dilution_effect.csv"))

  (strength_vs_coefs_power <- bounds_and_ATE_unnested %>%
      filter(scenario == "power") %>%
      unnest_wider(bounds) %>%
      mutate(strength = map_dbl(sum_stats, ~max(abs(outer(.x$thetas, .x$thetas, "-")))),
             U_on_XY = paste0("gamma[U]==", U_on_XY),
             p = factor(paste0("p==", p),
                        levels = paste0("p==", c(10,50))),
             many_weak = if_else(many_weak, "Many Weak", "Uniform")) %>%
      ggplot(aes(x = strength, y = indIVs_on_X, color = many_weak)) +
      geom_point() +
      facet_grid(p ~ U_on_XY,
                 labeller = label_parsed) +
      labs(
        y = bquote(gamma[j]),
        x = "Strength",
        color = "Simulation Scenario"
      ) +
      lims(
        x = c(0, 1),
        y = c(0, 4)
      ) +
      scale_color_manual(values = c("red", "black")) +
      theme_bw() +
      theme(legend.position = "top",
            axis.text.x = element_text(angle = 45, hjust = 1))
  )

  ggsave(here::here("figures/strength_vs_coef_multiple_IVs_power.png"),
         strength_vs_coefs_power,
         width = 4, height = 4, dpi = 300)

  (strength_vs_coefs_MR <- bounds_and_ATE_unnested %>%
      filter(scenario == "MR", X_on_Y == 2) %>%
      unnest_wider(bounds) %>%
      mutate(strength = map_dbl(sum_stats, ~max(abs(outer(.x$thetas, .x$thetas, "-")))),
             U_on_XY = paste0("gamma[U]==", U_on_XY),
             p = factor(paste0("p==", p),
                        levels = paste0("p==", c(10,50))),
             many_weak = if_else(many_weak, "Many Weak", "Uniform")) %>%
      ggplot(aes(x = strength, y = indIVs_on_X, color = many_weak)) +
      geom_point() +
      facet_grid(p ~ U_on_XY,
                 labeller = label_parsed) +
      labs(
        y = bquote(gamma[j]),
        x = "Strength",
        color = "Simulation Scenario"
      ) +
      scale_color_manual(values = c("red", "black")) +
      theme_bw() +
      theme(legend.position = "top",
            axis.text.x = element_text(angle = 45, hjust = 1))
  )

  ggsave(here::here("figures/strength_vs_coef_multiple_IVs_MR.png"),
         strength_vs_coefs_MR,
         width = 4, height = 4, dpi = 300)


  bounds_plot_power <- bounds_and_ATE %>%
    unnest(subset) %>%
    select(-data_file, -indIVs_on_X) %>%
    unnest(sums_and_bounds) %>%
    unnest_wider(bounds) %>%
    filter(scenario == "power",
           !many_weak) %>%
    mutate(U_on_XY = paste0("gamma[U]==", U_on_XY),
           p = factor(paste0("p==", p),
                      levels = paste0("p==", c(10,50))),
           X_on_Y = paste0("beta[X]==", X_on_Y)) %>%
    ggplot(aes(y = indIVs_on_X)) +
    geom_errorbar(aes(xmin = lower, xmax = upper),
                  width = 0.05) +
    geom_vline(xintercept = 0) +
    geom_vline(aes(xintercept = ATE,
                   color = "ATE")) +
    labs(
      y = bquote(gamma[j]),
      x = "ATE",
      color = ""
    ) +
    scale_x_continuous(limits = c(-1, 1),
                       breaks = c(-1, 0, 1)) +
    scale_color_manual(values = c("ATE" = "blue")) +
    facet_grid(U_on_XY ~ X_on_Y + p,
               labeller = label_parsed) +
    theme_bw() +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1))

  bounds_plot_power_many_weak <- bounds_and_ATE %>%
    unnest(subset) %>%
    select(-data_file, -indIVs_on_X) %>%
    unnest(sums_and_bounds) %>%
    unnest_wider(bounds) %>%
    filter(scenario == "power",
           many_weak) %>%
    mutate(U_on_XY = paste0("gamma[U]==", U_on_XY),
           p = factor(paste0("p==", p),
                      levels = paste0("p==", c(10,50))),
           X_on_Y = paste0("beta[X]==", X_on_Y)) %>%
    ggplot(aes(y = indIVs_on_X)) +
    geom_errorbar(aes(xmin = lower, xmax = upper),
                  width = 0.05) +
    geom_vline(xintercept = 0) +
    geom_vline(aes(xintercept = ATE,
                   color = "ATE")) +
    labs(
      y = bquote(gamma[j]),
      x = "ATE",
      color = ""
    ) +
    scale_x_continuous(limits = c(-1,1),
                       breaks = c(-1, 0, 1)) +
    scale_color_manual(values = c("ATE" = "blue")) +
    facet_grid(U_on_XY ~ X_on_Y + p,
               labeller = label_parsed) +
    theme_bw() +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(plot = bounds_plot_power,
         filename = here::here("figures/bounds_from_multiple_IV_sims_power.png"),
         dpi = 300, height = 7.5, width = 10)

  ggsave(plot = bounds_plot_power_many_weak,
         filename = here::here("figures/bounds_from_multiple_IV_sims_power_many_weak.png"),
         dpi = 300, height = 7.5, width = 10)


  bounds_plot_MR <- bounds_and_ATE_unnested %>%
    unnest_wider(bounds) %>%
    filter(scenario == "MR",
           !many_weak) %>%
    mutate(U_on_XY = paste0("gamma[U]==", U_on_XY),
           p = factor(paste0("p==", p),
                      levels = paste0("p==", c(10,50))),
           X_on_Y = paste0("beta[X]==", X_on_Y)) %>%
    ggplot(aes(y = indIVs_on_X)) +
    geom_errorbar(aes(xmin = lower, xmax = upper),
                  width = 0.005) +
    geom_vline(xintercept = 0) +
    geom_vline(aes(xintercept = ATE,
                   color = "ATE")) +
    labs(
      y = bquote(gamma[j]),
      x = "ATE",
      color = ""
    ) +
    scale_x_continuous(limits = c(-1,1),
                       breaks = c(-1, 0, 1)) +
    scale_color_manual(values = c("ATE" = "blue")) +
    facet_grid(U_on_XY ~ X_on_Y + p,
               labeller = label_parsed) +
    theme_bw() +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1))

  bounds_plot_MR_many_weak <- bounds_and_ATE_unnested %>%
    unnest_wider(bounds) %>%
    filter(scenario == "MR",
           many_weak) %>%
    mutate(U_on_XY = paste0("gamma[U]==", U_on_XY),
           p = factor(paste0("p==", p),
                      levels = paste0("p==", c(10,50))),
           X_on_Y = paste0("beta[X]==", X_on_Y)) %>%
    ggplot(aes(y = indIVs_on_X)) +
    geom_errorbar(aes(xmin = lower, xmax = upper),
                  width = 0.005) +
    geom_vline(xintercept = 0) +
    geom_vline(aes(xintercept = ATE,
                   color = "ATE")) +
    labs(
      y = bquote(gamma[j]),
      x = "ATE",
      color = ""
    ) +
    scale_x_continuous(limits = c(-1,1),
                       breaks = c(-1, 0, 1)) +
    scale_color_manual(values = c("ATE" = "blue")) +
    facet_grid(U_on_XY ~ X_on_Y + p,
               labeller = label_parsed) +
    theme_bw() +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(plot = bounds_plot_MR,
         filename = here::here("figures/bounds_from_multiple_IV_sims_MR.png"),
         dpi = 300, height = 7.5, width = 10)

  ggsave(plot = bounds_plot_MR_many_weak,
         filename = here::here("figures/bounds_from_multiple_IV_sims_MR_many_weak.png"),
         dpi = 300, height = 7.5, width = 10)

  ## Figures with subset of panels for mains
  bounds_plot_power_subset <- bounds_and_ATE %>%
    unnest(subset) %>%
    select(-data_file, -indIVs_on_X) %>%
    unnest(sums_and_bounds) %>%
    unnest_wider(bounds) %>%
    filter(scenario == "power",
           X_on_Y %in% c(0.5, 1.5),
           !many_weak) %>%
    mutate(U_on_XY = paste0("gamma[U]==", U_on_XY),
           p = factor(paste0("p==", p),
                      levels = paste0("p==", c(10,50))),
           X_on_Y = paste0("beta[X]==", X_on_Y)) %>%
    ggplot(aes(y = indIVs_on_X)) +
    geom_errorbar(aes(xmin = lower, xmax = upper),
                  width = 0.05) +
    geom_vline(xintercept = 0) +
    geom_vline(aes(xintercept = ATE,
                   color = "ATE")) +
    labs(
      y = bquote(gamma[j]),
      x = "ATE",
      color = ""
    ) +
    scale_x_continuous(limits = c(-1,1),
                       breaks = c(-1, 0, 1)) +
    scale_color_manual(values = c("ATE" = "blue")) +
    facet_grid(U_on_XY ~ X_on_Y + p,
               labeller = label_parsed) +
    theme_bw() +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1))

  bounds_plot_power_many_weak_subset <- bounds_and_ATE %>%
    unnest(subset) %>%
    select(-data_file, -indIVs_on_X) %>%
    unnest(sums_and_bounds) %>%
    unnest_wider(bounds) %>%
    filter(scenario == "power",
           many_weak,
           X_on_Y %in% c(0.5, 1.5)) %>%
    mutate(U_on_XY = paste0("gamma[U]==", U_on_XY),
           p = factor(paste0("p==", p),
                      levels = paste0("p==", c(10,50))),
           X_on_Y = paste0("beta[X]==", X_on_Y)) %>%
    ggplot(aes(y = indIVs_on_X)) +
    geom_errorbar(aes(xmin = lower, xmax = upper),
                  width = 0.05) +
    geom_vline(xintercept = 0) +
    geom_vline(aes(xintercept = ATE,
                   color = "ATE")) +
    labs(
      y = bquote(gamma[j]),
      x = "ATE",
      color = ""
    ) +
    scale_x_continuous(limits = c(-1,1),
                       breaks = c(-1, 0, 1)) +
    scale_color_manual(values = c("ATE" = "blue")) +
    facet_grid(U_on_XY ~ X_on_Y + p,
               labeller = label_parsed) +
    theme_bw() +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(plot = bounds_plot_power_subset,
         filename = here::here("figures/bounds_from_multiple_IV_sims_power_subset.png"),
         dpi = 300, height = 5, width = 4)

  ggsave(plot = bounds_plot_power_many_weak_subset,
         filename = here::here("figures/bounds_from_multiple_IV_sims_power_many_weak_subset.png"),
         dpi = 300, height = 5, width = 4)


  bounds_plot_MR_subset <- bounds_and_ATE_unnested %>%
    unnest_wider(bounds) %>%
    filter(scenario == "MR",
           X_on_Y %in% c(0.5, 1.5),
           !many_weak) %>%
    mutate(U_on_XY = paste0("gamma[U]==", U_on_XY),
           p = factor(paste0("p==", p),
                      levels = paste0("p==", c(10,50))),
           X_on_Y = paste0("beta[X]==", X_on_Y)) %>%
    ggplot(aes(y = indIVs_on_X)) +
    geom_errorbar(aes(xmin = lower, xmax = upper),
                  width = 0.005) +
    geom_vline(xintercept = 0) +
    geom_vline(aes(xintercept = ATE,
                   color = "ATE")) +
    labs(
      y = bquote(gamma[j]),
      x = "ATE",
      color = ""
    ) +
    scale_x_continuous(limits = c(-1,1),
                       breaks = c(-1, 0, 1)) +
    scale_color_manual(values = c("ATE" = "blue")) +
    facet_grid(U_on_XY ~ X_on_Y + p,
               labeller = label_parsed) +
    theme_bw() +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1))

  bounds_plot_MR_many_weak_subset <- bounds_and_ATE_unnested %>%
    unnest_wider(bounds) %>%
    filter(scenario == "MR",
           X_on_Y %in% c(0.5, 1.5),
           many_weak) %>%
    mutate(U_on_XY = paste0("gamma[U]==", U_on_XY),
           p = factor(paste0("p==", p),
                      levels = paste0("p==", c(10,50))),
           X_on_Y = paste0("beta[X]==", X_on_Y)) %>%
    ggplot(aes(y = indIVs_on_X)) +
    geom_errorbar(aes(xmin = lower, xmax = upper),
                  width = 0.005) +
    geom_vline(xintercept = 0) +
    geom_vline(aes(xintercept = ATE,
                   color = "ATE")) +
    labs(
      y = bquote(gamma[j]),
      x = "ATE",
      color = ""
    ) +
    scale_x_continuous(limits = c(-1,1),
                       breaks = c(-1, 0, 1)) +
    scale_color_manual(values = c("ATE" = "blue")) +
    facet_grid(U_on_XY ~ X_on_Y + p,
               labeller = label_parsed) +
    theme_bw() +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(plot = bounds_plot_MR_subset,
         filename = here::here("figures/bounds_from_multiple_IV_sims_MR_subset.png"),
         dpi = 300, height = 5, width = 4)

  ggsave(plot = bounds_plot_MR_many_weak_subset,
         filename = here::here("figures/bounds_from_multiple_IV_sims_MR_many_weak_subset.png"),
         dpi = 300, height = 5, width = 4)

  ## Create tables
  # summaries_not_many_weak <- bounds_and_ATE %>%
  #   unnest(subset) %>%
  #   filter(!many_weak) %>%
  #   select(-data_file, -indIVs_on_X) %>%
  #   unnest(sums_and_bounds) %>%
  #   unnest_wider(bounds) %>%
  #   mutate(width = upper - lower,
  #          across(where(is.numeric), round, digits = 3)) %>%
  #   group_by(p, scenario, X_on_Y, ATE, U_on_XY) %>%
  #   summarize(
  #     N = n(),
  #     across(c(lower, upper, width),
  #            list(average = mean, min = min, max = max)),
  #     .groups = "drop"
  #   ) %>%
  #   mutate(
  #     lower_sum = paste0(sprintf(lower_average, fmt = "%.3f"), " (", sprintf(lower_min, fmt = "%.3f"), " to ", sprintf(lower_max, fmt = "%.3f"), ")"),
  #     upper_sum = paste0(sprintf(upper_average, fmt = "%.3f"), " (", sprintf(upper_min, fmt = "%.3f"), " to ", sprintf(upper_max, fmt = "%.3f"), ")"),
  #     width_sum = paste0(sprintf(width_average, fmt = "%.3f"), " (", sprintf(width_min, fmt = "%.3f"), " to ", sprintf(width_max, fmt = "%.3f"), ")")
  #   ) %>%
  #   select(-c(lower_average:width_max))
  #
  # summaries_not_many_weak %>%
  #   select(-lower_sum, -upper_sum, -N, -ATE) %>%
  #   pivot_wider(names_from = X_on_Y,
  #               values_from = width_sum) %>%
  #   mutate(p = if_else(p == lag(p, default = 0), "", as.character(p)),
  #          scenario = if_else(scenario == lag(scenario, default = ""), "", scenario))
  #
  # summaries_not_many_weak %>%
  #   select(-width_sum, -upper_sum, -N, -ATE) %>%
  #   pivot_wider(names_from = X_on_Y,
  #               values_from = lower_sum) %>%
  #   mutate(p = if_else(p == lag(p, default = 0), "", as.character(p)),
  #          scenario = if_else(scenario == lag(scenario, default = ""), "", scenario))
  #
  #
  # summaries_many_weak <- bounds_and_ATE %>%
  #   unnest(subset) %>%
  #   filter(many_weak) %>%
  #   select(-data_file, -indIVs_on_X) %>%
  #   unnest(sums_and_bounds) %>%
  #   unnest_wider(bounds) %>%
  #   mutate(width = upper - lower,
  #          across(where(is.numeric), round, digits = 3)) %>%
  #   group_by(p, scenario, X_on_Y, ATE, U_on_XY) %>%
  #   summarize(
  #     N = n(),
  #     across(c(lower, upper, width),
  #            list(average = mean, min = min, max = max)),
  #     .groups = "drop"
  #   ) %>%
  #   mutate(
  #     lower_sum = paste0(sprintf(lower_average, fmt = "%.3f"), " (", sprintf(lower_min, fmt = "%.3f"), " to ", sprintf(lower_max, fmt = "%.3f"), ")"),
  #     upper_sum = paste0(sprintf(upper_average, fmt = "%.3f"), " (", sprintf(upper_min, fmt = "%.3f"), " to ", sprintf(upper_max, fmt = "%.3f"), ")"),
  #     width_sum = paste0(sprintf(width_average, fmt = "%.3f"), " (", sprintf(width_min, fmt = "%.3f"), " to ", sprintf(width_max, fmt = "%.3f"), ")")
  #   ) %>%
  #   select(-c(lower_average:width_max))
}