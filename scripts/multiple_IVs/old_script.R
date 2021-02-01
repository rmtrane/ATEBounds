library(tidyverse)
library(ACEBounds)
library(distributions3)
library(furrr)

set.seed(8072659)

sample_size <- 1e6
SIM_AND_SAVE <- FALSE
PLOT <- TRUE
OVERWRITE <- FALSE

ATE_from_simulated_data <- function(from_simulate_data) {
  intercept <- filter(from_simulate_data$coefficients, effect == "Yintercept")$coef
  x_beta <- filter(from_simulate_data$coefficients, effect == "X_on_Y")$coef
  u_beta <- filter(from_simulate_data$coefficients, effect == "U_on_Y")$coef

  pY1X0 <- 1 / (1 + exp(-intercept - u_beta * from_simulate_data$simulated_data$U))
  pY1X1 <- 1 / (1 + exp(-intercept - x_beta - u_beta * from_simulate_data$simulated_data$U))

  return(mean(pY1X1 - pY1X0))
}

if (SIM_AND_SAVE) {

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

    # if(all_combinations[i,]$p == 50){
    #   plan(multicore, workers = 1)
    # } else {
    #   plan(multicore, 3)
    # }

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
}


if (PLOT) {
  if (!file.exists(here::here("data/multiple_IV_sims/bounds_and_ATE.Rds")) | OVERWRITE) {
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
  } else {
    bounds_and_ATE <- read_rds(here::here("data/multiple_IV_sims/bounds_and_ATE.Rds"))
  }

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
  summaries_not_many_weak <- bounds_and_ATE %>%
    unnest(subset) %>%
    filter(!many_weak) %>%
    select(-data_file, -indIVs_on_X) %>%
    unnest(sums_and_bounds) %>%
    unnest_wider(bounds) %>%
    mutate(width = upper - lower,
           across(where(is.numeric), round, digits = 3)) %>%
    group_by(p, scenario, X_on_Y, ATE, U_on_XY) %>%
    summarize(
      N = n(),
      across(c(lower, upper, width),
             list(average = mean, min = min, max = max)),
      .groups = "drop"
    ) %>%
    mutate(
      lower_sum = paste0(sprintf(lower_average, fmt = "%.3f"), " (", sprintf(lower_min, fmt = "%.3f"), " to ", sprintf(lower_max, fmt = "%.3f"), ")"),
      upper_sum = paste0(sprintf(upper_average, fmt = "%.3f"), " (", sprintf(upper_min, fmt = "%.3f"), " to ", sprintf(upper_max, fmt = "%.3f"), ")"),
      width_sum = paste0(sprintf(width_average, fmt = "%.3f"), " (", sprintf(width_min, fmt = "%.3f"), " to ", sprintf(width_max, fmt = "%.3f"), ")")
    ) %>%
    select(-c(lower_average:width_max))

  summaries_not_many_weak %>%
    select(-lower_sum, -upper_sum, -N, -ATE) %>%
    pivot_wider(names_from = X_on_Y,
                values_from = width_sum) %>%
    mutate(p = if_else(p == lag(p, default = 0), "", as.character(p)),
           scenario = if_else(scenario == lag(scenario, default = ""), "", scenario))

  summaries_not_many_weak %>%
    select(-width_sum, -upper_sum, -N, -ATE) %>%
    pivot_wider(names_from = X_on_Y,
                values_from = lower_sum) %>%
    mutate(p = if_else(p == lag(p, default = 0), "", as.character(p)),
           scenario = if_else(scenario == lag(scenario, default = ""), "", scenario))



  summaries_many_weak <- bounds_and_ATE %>%
    unnest(subset) %>%
    filter(many_weak) %>%
    select(-data_file, -indIVs_on_X) %>%
    unnest(sums_and_bounds) %>%
    unnest_wider(bounds) %>%
    mutate(width = upper - lower,
           across(where(is.numeric), round, digits = 3)) %>%
    group_by(p, scenario, X_on_Y, ATE, U_on_XY) %>%
    summarize(
      N = n(),
      across(c(lower, upper, width),
             list(average = mean, min = min, max = max)),
      .groups = "drop"
    ) %>%
    mutate(
      lower_sum = paste0(sprintf(lower_average, fmt = "%.3f"), " (", sprintf(lower_min, fmt = "%.3f"), " to ", sprintf(lower_max, fmt = "%.3f"), ")"),
      upper_sum = paste0(sprintf(upper_average, fmt = "%.3f"), " (", sprintf(upper_min, fmt = "%.3f"), " to ", sprintf(upper_max, fmt = "%.3f"), ")"),
      width_sum = paste0(sprintf(width_average, fmt = "%.3f"), " (", sprintf(width_min, fmt = "%.3f"), " to ", sprintf(width_max, fmt = "%.3f"), ")")
    ) %>%
    select(-c(lower_average:width_max))
}
