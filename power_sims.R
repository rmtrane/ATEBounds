library(tidyverse)
library(ACEBounds)
library(distributions3)
library(furrr)


SIM_AND_SAVE <- TRUE
PLOT <- FALSE # TRUE

ATE_from_simulated_data <- function(from_simulate_data){
  intercept <- filter(from_simulate_data$coefficients, effect == "Yintercept")$coef
  x_beta <- filter(from_simulate_data$coefficients, effect == "X_on_Y")$coef
  u_beta <- filter(from_simulate_data$coefficients, effect == "U_on_Y")$coef

  pY1X0 <- 1 / (1 + exp(-intercept - u_beta*from_simulate_data$simulated_data$U))
  pY1X1 <- 1 / (1 + exp(-intercept - x_beta - u_beta*from_simulate_data$simulated_data$U))

  return(mean(pY1X1 - pY1X0))
}


if(SIM_AND_SAVE){
  plan(multicore, workers = 3)

  Z_on_X <- c(1, 1.5, 1.75, 1.9, 2, 2.1, 2.25, 2.5, 3, 3.5, 4)

  for (i in seq_along(Z_on_X)){
    many_sims <- expand_grid(indIVs_on_X = Z_on_X[i],
                             X_on_Y = c(1, 1.5, 2),
                             U_on_XY = c(0.1, 0.5)) %>%
      mutate(
        sim_data = future_pmap(
          .l = list(indIVs_on_X, X_on_Y, U_on_XY),
          function(x,y,z) simulate_data(
            sample_size = 1e7,
            IVs_ps = list(c(0.25, 0.5, 0.25)),
            X_intercept = -x,
            Y_intercept = -y/2,
            indIVs_on_X = x,
            indIVs_on_Y = 0,
            U = distributions3::Normal(),
            U_on_X = z,
            U_on_Y = z,
            X_on_Y = y
          ),
          .options = furrr_options(seed = TRUE),
          .progress = TRUE
        )
      )

    write_rds(
      many_sims,
      here::here("data", paste0("power_sims_", str_pad(i, width = 2, side = "left", pad = 0), ".Rds"))
    )
  }
}

if(PLOT){
  if(!file.exists(here::here("data/power_bounds_and_ATE.Rds"))){

  all_power_sims <- list.files(here::here("data"), pattern = "power_sims", full.names = TRUE)

  bounds_and_ATE <- tibble(data_file = all_power_sims) %>%
    mutate(
      subset = map(
        data_file,
        function(x){
          tmp <- read_rds(x)

          out <- tmp %>%
            mutate(
              sum_stats = future_map(sim_data, ~c(probs_from_data(.x$simulated_data, X, Y, Z1, data_format = "bivariate"),
                                                 ATE = ATE_from_simulated_data(.x))),
              get_bounds_res = future_map(sum_stats, ~get_bounds(gammas = .x$gammas, thetas = .x$thetas, stop = FALSE), .progress = TRUE),
              bounds = future_map(get_bounds_res, "interval")
            ) %>%
            select(-sim_data)

          return(out)
        })
    )

  write_rds(bounds_and_ATE, here::here("data/power_bounds_and_ATE.Rds"))

  } else {
    bounds_and_ATE <- read_rds(here::here("data/power_bounds_and_ATE.Rds"))
  }

  ## Compare indIVs_on_X and resulting strength
  coefs_vs_strength <- bounds_and_ATE %>%
    unnest(subset) %>%
    unnest_wider(sum_stats) %>%
    unnest_wider(bounds) %>%
    mutate(strength = map_dbl(thetas, ~.x[3] - .x[1]),
           U_on_XY = as.character(U_on_XY)) %>%
    ggplot(aes(x = indIVs_on_X, y = strength,
               color = U_on_XY)) +
      geom_point() +
      scale_color_manual(values = c("black", "red")) +
      labs(color = bquote(beta[U]),
           y = "Strength of IV",
           x = bquote(beta[Z])) +
      lims(
        y = c(0, 1)
      ) +
      theme_bw() +
      theme(legend.position = "top")

  ggsave(here::here("figures/MR_coefs_vs_strength.png"),
         coefs_vs_strength, dpi = 300,
         width = 4, height = 4)

  ## Bounds Plot
  pretty_plot <- bounds_and_ATE %>%
    unnest(subset) %>%
    unnest_wider(sum_stats) %>%
    unnest_wider(bounds) %>%
    mutate(zero = if_else(lower < 0 & upper > 0,
                          "Overlaps Zero", "Does Not Overlap Zero")) %>%
    arrange(desc(lower)) %>%
    rename(
      "Coefficient of X on Y" = X_on_Y,
      "Coefficient of U on both X and Y" = U_on_XY,
    ) %>%
    ggplot(aes(x = indIVs_on_X, color = zero)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_hline(aes(yintercept = ATE, color = "ATE")) +
      geom_errorbar(aes(ymin = lower, ymax = upper)) +
      facet_grid(`Coefficient of U on both X and Y` ~ `Coefficient of X on Y`,
                 labeller = label_both) +
      lims(y = c(-1, 1)) +
      scale_color_manual(
        values = c("Overlaps Zero" = "black", "Does Not Overlap Zero" = "red", "ATE" = "blue")
      ) +
      labs(
        x = "Coefficient for effect of Z on X",
        y = "ATE",
        color = ""
      ) +
      theme_bw() +
      theme(legend.position = "top")

  ggsave(here::here("figures/power.png"),
         plot = pretty_plot, dpi = 300,
         height = 6, width = 8, units = "in")


  power_curves <- bounds_and_ATE %>%
    unnest(subset) %>%
    unnest_wider(sum_stats) %>%
    unnest_wider(bounds) %>%
    filter(lower > 0) %>%
    group_by(U_on_XY, X_on_Y) %>%
    filter(indIVs_on_X == min(indIVs_on_X)) %>%
    mutate(U_on_XY = as.character(U_on_XY)) %>%
    ggplot(aes(x = ATE, y = indIVs_on_X, color = U_on_XY)) +
      geom_point() + geom_line() +
      lims(y = c(1.5, 3),
           x = c(0.2, 0.5)) +
      labs(
        x = "Average Treatment Effect",
        y = bquote(beta[Z]),
        color = bquote(beta[U])
      ) +
      scale_color_manual(values = c("black", "red")) +
      theme_bw() +
      theme(legend.position = "top")

  ggsave(here::here("figures/power_curves.png"),
         plot = power_curves, dpi = 300,
         height = 4, width = 4)


  loess_strength <- bounds_and_ATE %>%
    unnest(subset) %>%
    unnest_wider(sum_stats) %>%
    unnest_wider(bounds) %>%
    mutate(strength = map_dbl(thetas, ~.x[3] - .x[1])) %>%
    loess(data = .,
          strength ~ indIVs_on_X + U_on_XY)

  uniroot(f = function(x) predict(loess_strength, newdata = data.frame(indIVs_on_X = x, U_on_XY = 0.1)) - 0.5,
          interval = c(1, 1.5))

  ## Model lower bounds by ATE
  loess_model <- bounds_and_ATE %>%
    unnest(subset) %>%
    unnest_wider(sum_stats) %>%
    unnest_wider(bounds) %>%
    loess(data = .,
          lower ~ indIVs_on_X + X_on_Y + U_on_XY)

  ## Show on plot
  loess_plot <- bounds_and_ATE %>%
    unnest(subset) %>%
    unnest_wider(sum_stats) %>%
    unnest_wider(bounds) %>%
    ggplot(aes(x = indIVs_on_X, y = lower)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point() +
      geom_smooth() +
      facet_grid(U_on_XY ~ X_on_Y,
                 labeller = label_both) +
      lims(y = c(-1, 1)) +
      labs(
        x = "Coefficient for effect of Z on X",
        y = "ATE",
        color = ""
      ) +
      theme_bw() +
      theme(legend.position = "top")

  ggsave(here::here("figures/power_loess_lower_bound.png"),
         loess_plot, dpi = 300,
         height = 6, width = 8, units = "in")

  ## Find roots for loess curves
  roots <- tmp %>%
    select(X_on_Y, U_on_XY, ATE) %>%
    mutate(ATE = round(ATE, digits = 7)) %>%
    unique() %>%
    rowwise() %>%
    mutate(uniroot_res = list(uniroot(f = function(x) predict(tmp_loess, newdata = data.frame(X_on_Y = X_on_Y, U_on_XY = U_on_XY, indIVs_on_X = x)),
                                      interval = c(1.6, 2.4))),
           indIVs_on_X = uniroot_res$root) %>%
    ungroup()


  ## Power curves including for loess model
  power_with_loess_ests <- bounds_and_ATE %>%
    unnest(subset) %>%
    unnest_wider(sum_stats) %>%
    unnest_wider(bounds) %>%
    filter(lower > 0) %>%
    group_by(U_on_XY, X_on_Y) %>%
    filter(indIVs_on_X == min(indIVs_on_X)) %>%
    ggplot(aes(x = ATE, y = indIVs_on_X, group = U_on_XY, color = "From simulated data")) +
      geom_point() + geom_line() +
      geom_point(data = roots, aes(color = "from loess")) + geom_line(data = roots, aes(color = "from loess")) +
      lims(y = c(1.5, 3)) +
      scale_color_manual("Smallest coefficient of Z on X with lower bound > 0",
                         values = c("black", "red")) +
      facet_grid(~U_on_XY,
                 labeller = label_both) +
      guides(
        color = guide_legend(title.position = "top")
      ) +
      theme_bw() +
      theme(legend.position = "top")

  ggsave(here::here("figures/power_with_loess_ests.png"),
         power_with_loess_ests, dpi = 300,
         height = 5, width = 8)
}



# all_sims0 <- tibble(sims_for_power = list.files(pattern = "sims_for_power")) %>%
#   rowwise() %>%
#   mutate(sims = list(read_rds(sims_for_power))) %>%
#   ungroup() %>%
#   unnest(sims)
#
# plot_DAG(all_sims$sim_data[[1]])
#
# all_sims_w_bounds0 <- many_sims %>% #all_sims %>%
#   mutate(simulated_data = map(sim_data, "simulated_data"),
#          sum_stats = map(simulated_data, ~probs_from_data(.x, X, Y, Z1, data_format = "bivariate")),
#          get_bounds_res = map(sum_stats, ~get_bounds(gammas = .x$gammas, thetas = .x$thetas, stop = FALSE)),
#          bounds = map(get_bounds_res, "interval"))
#
#
# all_sims_w_bounds0 %>%
#   unnest_wider(bounds) %>%
#   #filter(i == 1) %>%
#   mutate(zero = lower < 0 & upper > 0) %>%
#   ggplot(aes(x = indIVs_on_X, color = zero)) +
#     geom_hline(yintercept = 0, linetype = "dashed") +
#     geom_errorbar(aes(ymin = lower, ymax = upper)) +
#     facet_wrap(~ X_on_Y, labeller = label_both,
#                nrow = 3) +
#     scale_color_manual(
#       values = c("TRUE" = "black", "FALSE" = "red")
#     )
#
#
# all_sims_w_bounds %>%
#   unnest_wider(bounds) %>%
#   group_by(X_on_Y, indIVs_on_X) %>%
#   summarize(
#     s_lower = sd(lower),
#     s_upper = sd(upper),
#     lower_range = max(lower) - min(lower),
#     upper_range = max(upper) - min(upper)
#   )
