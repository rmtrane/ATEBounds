params <-
list(n_cores = 2L)

## ----message = FALSE----------------------------------------------------------
library(tidyverse)
library(ACEBounds)

## ----echo = FALSE-------------------------------------------------------------
ACEBounds:::matrices_from_polymake %>% 
  filter(!x_monotone, !y_monotone,
         data_format == "bivariate", 
         n_z_levels == 3) %>% 
  pull(matrix) %>% .[[1]] %>% 
  filter(alpha == 0) %>%
  select(where(~sum(abs(.x)) > 0)) %>%
  filter(rowSums(abs(.)) > 1) %>% 
  pander::pander()

## -----------------------------------------------------------------------------
simulate_gammas_from_thetas <- function(thetas){

  gammas <- vector(length = length(thetas))


  gammas[1] <- runif(1)

  gammas[2] <- runif(1,
                     min = max(0,
                               gammas[1] - thetas[1] - thetas[2], # row 3
                               gammas[1] + thetas[1] + thetas[2] - 2), # row 5
                     max = min(1,
                               gammas[1] - thetas[1] - thetas[2] + 2, # row 11
                               gammas[1] + thetas[1] + thetas[2]) # row 16
                     )

  gammas[3] <- runif(1,
                     min = max(0,
                               gammas[1] - thetas[1] - thetas[3], # row 4
                               gammas[1] + thetas[1] + thetas[3] - 2, # row 6
                               gammas[2] - thetas[2] - thetas[3], # row 12
                               gammas[2] + thetas[2] + thetas[3] - 2), # row 9
                     max = min(1,
                               gammas[2] - thetas[2] - thetas[3] + 2, # row 10
                               gammas[2] + thetas[2] + thetas[3], # row 13
                               gammas[1] - thetas[1] - thetas[3] + 2, # row 17
                               gammas[1] + thetas[1] + thetas[3])) # row 15
  return(gammas)
}

## ----many_sample_joints-------------------------------------------------------
if(file.exists(here::here("vignettes_data/many_sample_joints.rds"))){
  many_sample_joints <- read_rds(here::here("vignettes_data/many_sample_joints.rds"))
} else {
  
  library(furrr)
  
  plan(multisession, workers = params$n_cores)
 
  set.seed(7226637)
  sim_probs <- tibble(j = 1:100) %>%
    mutate(thetas = map(j, ~runif(n = 3, min = 0, max = 1)), 
           gammas = map(thetas, simulate_gammas_from_thetas),
           pot_covs = map2(thetas, gammas, potential_covs))
  
  many_sample_joints <- sim_probs %>%
    mutate(
      samp_joints = future_map2(pot_covs, j,
                                ~sample_joint_probs(.x, return_bounds = TRUE, n = 1000, max_rejections = 500, 
                                                    print_progress = TRUE, print_as_progress = .y),
                               .options = furrr_options(seed = TRUE))
    )
  
  write_rds(many_sample_joints,
            here::here("vignettes_data/many_sample_joints.rds"))
}

head(many_sample_joints)

## ----trivariate_bounds--------------------------------------------------------
trivariate_bounds <- many_sample_joints %>% 
  unnest(samp_joints) %>% 
  unnest(joint) %>% 
  unnest_wider(bounds) %>% 
  select(j, joint, trivariate_lower = lower, trivariate_upper = upper, n_rejected)

head(trivariate_bounds)

## ----bivariate_bounds---------------------------------------------------------
bivariate_bounds <- many_sample_joints %>% 
  rowwise() %>% 
  mutate(bounds = list(get_bounds(thetas = thetas, gammas = gammas, stop = FALSE, warning = FALSE)),
         bivariate_constraints_violated = bounds$constraints_violated,
         interval = list(bounds$interval)) %>% 
  ungroup() %>% 
  unnest_wider(interval) %>% 
  rename(bivariate_lower = lower,
         bivariate_upper = upper) %>% 
  select(j, contains("bivariate"), thetas, gammas)

head(bivariate_bounds)

## ----for_plot-----------------------------------------------------------------
for_plot <- trivariate_bounds %>% 
  left_join(bivariate_bounds %>% 
              mutate(center = (bivariate_upper + bivariate_lower) / 2) %>% 
              arrange(center) %>% 
              mutate(k = row_number())) %>% 
  arrange(j, trivariate_lower) %>% 
  group_by(j) %>% 
  mutate(id = row_number()) %>% 
  ungroup() %>% 
  mutate(contains_zero = if_else(trivariate_lower < 0 & trivariate_upper > 0, "Overlaps Zero", "Does not Overlap Zero"),
         facet_row = ceiling(k/10), 
         facet_col = k - (facet_row-1)*10)

head(for_plot)

## ----plot, fig.width = 10, fig.height = 8, dpi = 150--------------------------
(plot <- for_plot %>% 
  ggplot(aes(y = id/max(id))) + 
    geom_vline(aes(xintercept = bivariate_lower, color = "Two-sample Bounds")) + 
    geom_vline(aes(xintercept = bivariate_upper, color = "Two-sample Bounds")) + 
    geom_errorbar(aes(xmin = trivariate_lower, 
                      xmax = trivariate_upper,
                      color = contains_zero)) + 
    geom_text(data = for_plot %>% 
                select(j, k, contains("facet"), bivariate_upper) %>% 
                unique(),
              aes(label = j, y = 0.25, x = if_else(bivariate_upper < 0.6,
                                                   0.85, -0.85))) +
    xlim(c(-1,1)) +
    facet_grid(facet_row ~ facet_col) +
    coord_fixed() +
    scale_color_manual(values = c("black", "red4", "grey50")) + 
    labs(y = "", x = "ATE",
         title = "Bounds on ATE",
         color = "",
         linetype = ""
    ) + 
    theme_bw() +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank()))

## -----------------------------------------------------------------------------
ggsave(
  plot, 
  filename = here::here("figures/trivariate_bounds_plot.png"),
  width = 10, height = 8, dpi = 300
)

## -----------------------------------------------------------------------------
chosen_ones <- c(53, 6, 66, 
                 38, 7, 88,
                 73, 44, 60)

for_subset_plot <- for_plot %>% 
  filter(j %in% chosen_ones) %>% 
  mutate(J = as.numeric(factor(j, levels = chosen_ones)),
         row_i = case_when(J %in% 1:3 ~ "A",
                           J %in% 4:6 ~ "B",
                           J %in% 7:9 ~ "C",
                           TRUE ~ NA_character_),
         col_j = (J-1) %% 3 + 1)

subset_plot_summaries <- for_subset_plot %>% 
  group_by(row_i, col_j, thetas, gammas, J) %>% 
  summarize(bivariate_lower = unique(bivariate_lower),
            bivariate_upper = unique(bivariate_upper),
            p_no_zero = mean(contains_zero == "Does not Overlap Zero")) %>% 
  rowwise() %>% 
  mutate(thetas = list(setNames(thetas, paste0("P(X = 1 | Z = ", 0:2, ")"))),
         gammas = list(setNames(gammas, paste0("P(Y = 1 | Z = ", 0:2, ")")))) %>% 
  unnest_wider(thetas) %>% 
  unnest_wider(gammas)

subset_plot <- ggplot(for_subset_plot,
                      aes(y = id/max(id))) + 
  geom_rect(data = data.frame(x_min = rep(0.2, 3),
                              x_max = rep(0.3, 3),
                              y_min = rep(0.1, 3),
                              y_max = rep(0.2, 3),
                              contains_zero = c("Two-Sample Bounds", "Overlaps Zero", "Does not Overlap Zero")),
            inherit.aes = FALSE,
            aes(xmin = x_min,
                xmax = x_max,
                ymin = y_min,
                ymax = y_max,
                fill = contains_zero),
            alpha = 0) +
  geom_vline(aes(xintercept = bivariate_lower, color = "Two-Sample Bounds"),
             linetype = "dashed") + 
  geom_vline(aes(xintercept = bivariate_upper, color = "Two-Sample Bounds"),
             linetype = "dashed") + 
  geom_errorbar(aes(xmin = trivariate_lower, 
                    xmax = trivariate_upper,
                    color = contains_zero)) + 
  geom_vline(xintercept = 0) +
  xlim(c(-1,1)) +
  facet_grid(row_i ~ col_j) +
  coord_fixed(ratio = 2) +
  scale_color_manual(values = c("black", "grey50", "red"),
                     breaks = c("Two-Sample Bounds", "Overlaps Zero", "Does not Overlap Zero")) +
  scale_fill_manual(values = c("black", "grey50", "red"),
                    breaks = c("Two-Sample Bounds", "Overlaps Zero", "Does not Overlap Zero")) + 
  labs(x = "ATE", 
       y = "",
       color = "",
       fill = "", 
       linetype = ""
  ) + 
  guides(
    color = "none",
    fill = guide_legend(override.aes = list(alpha = 1))
  ) +
  theme_bw() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "top")

subset_plot +
  labs(title = "Bounds on ATE")

## -----------------------------------------------------------------------------
write_rds(x = subset_plot_summaries,
          file = here::here("vignettes_data", "subset_plot_summaries.Rds"))

ggsave(
  plot = subset_plot,
  filename = here::here("figures", "trivariate_bounds_subset_plot.png"),
  width = 8, dpi = 300, height = 8
)

## ----fig.width = 10, fig.height = 10, dpi = 300-------------------------------
simulate_gammas_from_thetas_mono <- function(thetas){

  gammas <- vector(length = length(thetas))

  ## constraints: 
  # matrices_from_polymake %>%
  #   filter(n_z_levels == 3, x_monotone, !y_monotone, data_format == "bivariate") %>%
  #   pull(matrix) %>% .[[1]] %>%
  #   filter(alpha == 0) %>%
  #   select(where(~sum(abs(.x)) > 0)) %>%
  #   filter(rowSums(abs(.)) > 1)
  
  gammas[1] <- runif(1)

  gammas[2] <- runif(1,
                     min = max(0,
                               gammas[1] - thetas[2] + thetas[1]), # row 8
                     max = min(1,
                               gammas[1] + thetas[2] - thetas[1])) # row 4

  gammas[3] <- runif(1,
                     min = max(0,
                               gammas[2] - gammas[1], # row 5
                               gammas[2] + thetas[2] - thetas[3]), # row 3
                     max = min(1,
                               gammas[2] + thetas[3] - thetas[2], # row 2
                               1 + gammas[2] - gammas[1])) # row 7
  return(gammas)
}

if(file.exists(here::here("vignettes_data/many_sample_joints_mono.rds"))){
  many_sample_joints_mono <- read_rds(here::here("vignettes_data/many_sample_joints_mono.rds"))
} else {
  
  library(furrr)

  if(!interactive())
    plan(multisession, workers = params$n_cores)
 
  set.seed(2884193)
  sim_probs_mono <- tibble(j = 1:100) %>%
    mutate(thetas = map(j, ~sort(runif(n = 3, min = 0, max = 1))), 
           gammas = map(thetas, simulate_gammas_from_thetas_mono),
           pot_covs = map2(thetas, gammas, potential_covs, x_mono = TRUE))
  
  many_sample_joints_mono <- sim_probs_mono %>% 
    mutate(
      samp_joints = future_map(pot_covs, 
                               sample_joint_probs, #x, 
                               return_bounds = TRUE, n = 1000, max_rejections = 100,
                               x_mono = TRUE,
                               .options = furrr_options(seed = TRUE))
    )
  
  write_rds(many_sample_joints_mono,
            here::here("vignettes_data/many_sample_joints_mono.rds"))
}


trivariate_bounds_mono <- many_sample_joints_mono %>% 
  unnest(samp_joints) %>% 
  unnest(joint) %>% 
  unnest_wider(bounds) %>% 
  select(j, joint, trivariate_lower = lower, trivariate_upper = upper, n_rejected)

bivariate_bounds_mono <- many_sample_joints_mono %>% 
  rowwise() %>% 
  mutate(bounds = list(get_bounds(thetas = thetas, gammas = gammas, x_mono = TRUE, stop = FALSE, warning = FALSE)),
         bivariate_constraints_violated = bounds$constraints_violated,
         interval = list(bounds$interval)) %>% 
  ungroup() %>% 
  unnest_wider(interval) %>% 
  rename(bivariate_lower = lower,
         bivariate_upper = upper) %>% 
  print() %>% 
  select(j, contains("bivariate"))

for_plot_mono <- trivariate_bounds_mono %>% 
  left_join(bivariate_bounds_mono) %>% 
  arrange(j, trivariate_lower) %>% 
  group_by(j) %>% 
  mutate(id = row_number()) %>% 
  ungroup() %>% 
  mutate(contains_zero = if_else(trivariate_lower < 0 & trivariate_upper > 0, "Overlaps Zero", "Does not Overalp Zero"))

plot_mono <- for_plot_mono %>% 
  mutate(facet_row = ceiling(j/10), 
         facet_col = j - (facet_row-1)*10) %>% 
  ggplot(aes(x = id/max(id))) + 
    geom_hline(aes(yintercept = bivariate_lower, color = "Two-Sample Bounds")) + 
    geom_hline(aes(yintercept = bivariate_upper, color = "Two-Sample Bounds")) + 
    geom_errorbar(aes(ymin = trivariate_lower, ymax = trivariate_upper,
                      color = contains_zero)) + 
    ylim(c(-1,1)) +
    facet_grid(facet_row ~ facet_col) +
    coord_fixed() +
    scale_color_manual(values = c("black", "grey50", "red4")) + 
    labs(x = "", y = "ACE",
         title = "Bounds on ACE",
         caption = "Assuming P(X = 1 | Z = z, U) is monotonically increasing.",
         color = "",
         linetype = ""
    ) + 
    theme_bw() +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank())

plot_mono

ggsave(plot_mono, filename = here::here("figures/trivariate_bounds_mono_plot.png"),
       width = 10, height = 10, dpi = 300)

