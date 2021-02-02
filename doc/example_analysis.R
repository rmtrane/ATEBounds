## -----------------------------------------------------------------------------
library(TwoSampleMR)
library(tidyverse)
library(magrittr)
library(ACEBounds)


theme_set(theme_bw())

## -----------------------------------------------------------------------------
ao <- available_outcomes()

## -----------------------------------------------------------------------------
smoking_experiment <- ao %>% filter(id == "ukb-d-20116_0")

## -----------------------------------------------------------------------------
smoking_instruments <- extract_instruments(smoking_experiment$id)

## -----------------------------------------------------------------------------
lung_cancer_experiment <- ao %>% filter(id == "ukb-d-40001_C349")
lung_cancer_data <- extract_outcome_data(snps = smoking_instruments$SNP, outcomes = lung_cancer_experiment$id)

## -----------------------------------------------------------------------------
smoking_lung_cancer_harmonized <- harmonise_data(smoking_instruments, lung_cancer_data) %>% as_tibble()

## -----------------------------------------------------------------------------
smoking_lung_cancer_pop_probs <- tibble(regression = c("outcome", "exposure"),
                                       p_outcome = c(lung_cancer_experiment  %$% { ncase / (ncase + ncontrol) },
                                                     smoking_experiment %$% { ncase / (ncase + ncontrol) }))

## -----------------------------------------------------------------------------
find_intercept <- function(beta1, p_outcome, pz, interval = c(-5, 5), zs = 0:2){
  uniroot(f = function(beta0) sum((1+exp(-beta0-beta1*zs))^(-1) * pz) - p_outcome,
          interval = interval)
}

## -----------------------------------------------------------------------------
smoking_lung_cancer_prep_for_bounds <- smoking_lung_cancer_harmonized %>%
  rowwise() %>% 
  mutate(samplesize.exposure = smoking_experiment$ncase + smoking_experiment$ncontrol,
         samplesize.outcome = lung_cancer_experiment$ncase + lung_cancer_experiment$ncontrol,
         ave_eaf = weighted.mean(x = c_across(c(eaf.exposure, eaf.outcome)),
                                 w = c_across(c(samplesize.exposure, samplesize.outcome)))) %>% 
  select(SNP, contains("beta"), ave_eaf) %>% 
  pivot_longer(contains("beta")) %>% 
  separate(name, into = c("variable", "regression"), sep = "\\.") %>% 
  pivot_wider(names_from = variable, values_from = value) %>% 
  mutate(`P(Z = 2)` = (1 - ave_eaf)^2,
         `P(Z = 1)` = 2*ave_eaf*(1 - ave_eaf),
         `P(Z = 0)` = ave_eaf^2)

## -----------------------------------------------------------------------------
smoking_lung_cancer_w_intercept <- smoking_lung_cancer_prep_for_bounds %>% 
  left_join(smoking_lung_cancer_pop_probs) %>% 
  rowwise() %>% 
  mutate(find_intercept = list(find_intercept(beta1 = beta,
                                              p_outcome = p_outcome,
                                              pz = c(`P(Z = 0)`, `P(Z = 1)`, `P(Z = 2)`),
                                              interval = c(-7,7))),
         beta0 = find_intercept$root,
         `thetas/gammas` = list(1/(1+exp(-beta0 - 0:2*beta))))

## ----message = FALSE, warning = FALSE-----------------------------------------
summaries_for_plots <- smoking_lung_cancer_w_intercept %>% 
  ungroup() %>% 
  select(SNP, regression, beta, starts_with("P(Z ="), p_outcome, beta0, `thetas/gammas`) %>% 
  mutate(`thetas/gammas` = map(`thetas/gammas`, ~setNames(.x, c("xy_given_z0", "xy_given_z1", "xy_given_z2")))) %>% 
  unnest_wider(`thetas/gammas`)

write_csv(
  summaries_for_plots,
  here::here("vignettes_data/summary_stats_smoking_on_lung_cancer.csv")
)

## -----------------------------------------------------------------------------
(marginal_Z <- summaries_for_plots %>% 
  filter(regression == "exposure") %>% 
  pivot_longer(cols = starts_with("P(Z ="),
               names_to = "z", values_to = "p") %>% 
  ggplot(aes(x = p)) +
    geom_histogram(boundary = 0, binwidth = 0.05, 
                   color = "black") +
    facet_grid(~ z) +
    xlim(c(0,1)) + 
    theme_bw())

ggsave(
  plot = marginal_Z,
  filename = here::here("figures/example_analyses", 
                        "smoking_lung_cancer_marginal_Z.png"),
  dpi = 300, width = 6.5, height = 2
)

## -----------------------------------------------------------------------------
(lung_cancer_coefficients <- summaries_for_plots %>% 
  pivot_longer(cols = c(beta, beta0),
               names_to = "coef", values_to = "value") %>% 
  mutate(real_coefs = case_when(coef == "beta" & regression == "exposure" ~ "beta[1]", 
                                coef == "beta0" & regression == "exposure" ~ "beta[0]",
                                coef == "beta" & regression == "outcome" ~ "gamma[1]", 
                                coef == "beta0" & regression == "outcome" ~ "gamma[0]")) %>% 
  ggplot(aes(x = value)) +
    geom_histogram(boundary = 0, binwidth = 0.001, 
                   color = "black") +
    facet_wrap(~ real_coefs,
               scales = "free",
               labeller = label_parsed) +
    theme_bw())

ggsave(
  plot = lung_cancer_coefficients,
  filename = here::here("figures/example_analyses", 
                        "smoking_lung_cancer_coefficients.png"),
  dpi = 300, width = 5, height = 3
)

## -----------------------------------------------------------------------------
(marginal_conditionals <- summaries_for_plots %>% 
  pivot_longer(cols = starts_with("xy_given_z"),
               names_to = "vars", values_to = "p") %>% 
  mutate(regression = if_else(regression == "exposure",
                              "P(X = 1 | Z = z)", "P(Y = 1 | Z = z)"),
         vars = paste("z =", str_extract(vars, "[0-2]"))) %>% 
  ggplot(aes(x = p)) + 
    geom_histogram(bins = 20, boundary = 0.2, 
                   color = "black") + 
    facet_grid(vars ~ regression,
               scales = "free") +
    theme_bw())

ggsave(
  plot = marginal_conditionals,
  filename = here::here("figures/example_analyses",
                        "smoking_lung_cancer_marginal_conditionals.png"),
  dpi = 300, width = 5, height = 3
)

## -----------------------------------------------------------------------------
smoking_lung_cancer_bounds <- smoking_lung_cancer_w_intercept %>% 
  select(-find_intercept, -beta0, -beta, -p_outcome, -ave_eaf) %>% 
  mutate(regression = if_else(regression == "outcome", "gammas", "thetas")) %>% 
  pivot_wider(names_from = regression, values_from = `thetas/gammas`) %>% 
  rowwise() %>% 
  mutate(bivariate_bounds = list(get_bounds(thetas = thetas, gammas = gammas, 
                                            stop = FALSE, warning = FALSE)),
         interval = list(bivariate_bounds$interval),
         constraints_violated = bivariate_bounds$constraints_violated, 
         bivariate_bounds_mono = list(get_bounds(thetas = thetas, gammas = gammas, 
                                                 stop = FALSE, warning = FALSE,
                                                 x_mono = TRUE)),
         interval_mono = list(bivariate_bounds_mono$interval),
         constraints_violated_mono = bivariate_bounds_mono$constraints_violated) %>% 
  ungroup() %>% 
  unnest_wider(col = interval_mono) %>% 
  rename(upper_mono = upper, lower_mono = lower) %>% 
  unnest_wider(interval)

## -----------------------------------------------------------------------------
smoking_lung_cancer_bivariate_bounds <- smoking_lung_cancer_bounds %>% 
  arrange(lower) %>% 
  ggplot(aes(y = SNP)) +
    geom_errorbar(aes(xmax = upper, xmin = lower)) +
    theme_bw() +
    facet_wrap(~ lower < median(lower),
               ncol = 2, scales = "free_y") +
    theme(strip.text = element_blank(),
          strip.background = element_blank()) + 
    labs(
      y = "",
      x = "ATE"
    ) +
    lims(x = c(-1, 1))

smoking_lung_cancer_bivariate_bounds +
      labs(title = "Bivariate bounds",
         subtitle = paste0("Exposure: Smoking (id: ", smoking_experiment$id, "). ",
                           "Outcome: Lung Cancer (id: ", lung_cancer_experiment$id, ")."))

ggsave(
  plot = smoking_lung_cancer_bivariate_bounds,
  filename = here::here("figures/example_analyses", 
                        paste0("smoking_lung_cancer_bivariate_bounds_", smoking_experiment$id, "_", lung_cancer_experiment$id,".png")),
  width = 8, height = 10, dpi = 300
)

## -----------------------------------------------------------------------------
smoking_SNP_subset <- unique(smoking_lung_cancer_bounds$SNP)[c(2, 5, 17, 25, 34, 52, 65, 78)]

smoking_lung_cancer_bivariate_bounds_subset <- smoking_lung_cancer_bounds %>% 
  filter(SNP %in% smoking_SNP_subset) %>% 
  arrange(lower) %>% 
  ggplot(aes(y = SNP)) +
    geom_errorbar(aes(xmax = upper, xmin = lower)) +
    theme_bw() +
    theme(strip.text = element_blank(),
          strip.background = element_blank()) + 
    labs(
      y = "",
      x = "ATE"
    ) +
    lims(x = c(-1, 1))

ggsave(
  plot = smoking_lung_cancer_bivariate_bounds_subset,
  filename = here::here("figures/example_analyses", 
                        paste0("smoking_lung_cancer_bivariate_bounds_subset_", smoking_experiment$id, "_", lung_cancer_experiment$id,".png")),
  width = 4, height = 4, dpi = 300
)

## -----------------------------------------------------------------------------
strength_histogram <- smoking_lung_cancer_bounds %>% 
  rowwise() %>% 
  mutate(strength = max(abs(outer(thetas, thetas, `-`)))) %>% 
  ungroup() %>%
  ggplot(aes(x = strength)) + 
    geom_histogram(color = "black", binwidth = 0.001, boundary = 0) + 
    scale_x_continuous(limits = c(0, 0.012), 
                       breaks = seq(0, 0.012, by = 0.002)) +
    labs(x = "Strength", y = "Count") +
    theme_bw()

strength_histogram +
  labs(title = "Strength of IVs on Smoking")

ggsave(
  plot = strength_histogram,
  filename = here::here("figures/example_analyses",
                        "strength_histogram.png"),
  dpi = 300, width = 8, height = 4
)

## -----------------------------------------------------------------------------
smoking_lung_cancer_filename <- paste0("samples_of_joints_", smoking_experiment$id, "_", lung_cancer_experiment$id, ".Rds")

if(file.exists(here::here("vignettes_data/example_analyses", smoking_lung_cancer_filename))){
  smoking_lung_cancer_samples_of_joints <- read_rds(here::here("vignettes_data/example_analyses", smoking_lung_cancer_filename))
} else {
  set.seed(7892312)
  
  smoking_lung_cancer_samples_of_joints <- smoking_lung_cancer_bounds %>% 
    rowwise() %>% 
    mutate(pot_covs = list(potential_covs(thetas = thetas, gammas = gammas)),
           sample_joint = list(sample_joint_probs(pot_covs, n = 500, return_bounds = TRUE)))
  
  write_rds(smoking_lung_cancer_samples_of_joints, here::here("vignettes_data/example_analyses", smoking_lung_cancer_filename))
}

## -----------------------------------------------------------------------------
smoking_lung_cancer_trivariate_bounds <- smoking_lung_cancer_samples_of_joints %>% 
  ungroup() %>% 
  select(SNP, bivariate_lower = lower, bivariate_upper = upper, sample_joint) %>% 
  unnest(sample_joint) %>% 
  unnest(joint) %>% 
  unnest_wider(bounds) %>% 
  mutate(contains_zero = lower < 0 & upper > 0)

## ----fig.width = 8------------------------------------------------------------
smoking_lung_cancer_individual_SNPs_plot <- smoking_lung_cancer_trivariate_bounds %>% 
  group_by(SNP) %>% 
  arrange(lower) %>% 
  mutate(id = row_number()) %>% 
  ungroup() %>% 
  ggplot(aes(y = id, color = contains_zero)) +
    geom_vline(aes(xintercept = bivariate_lower),
               linetype = "dashed") +  
    geom_vline(aes(xintercept = bivariate_upper),
               linetype = "dashed") +
    geom_vline(xintercept = 0) +
    geom_errorbar(aes(xmin = lower,
                      xmax = upper)) +
    geom_point(data = data.frame(x = 0, y = 10,
                                 contains_zero = FALSE),
               aes(x = x, y = y),
               alpha = 0) +
    scale_color_manual("",
                       values = c("grey50", "red4"),
                       labels = c("Overlaps Zero", "Does not Overlap Zero"),
                       breaks = c("TRUE", "FALSE"),
                       drop = FALSE) +
    lims(x = c(-1, 1)) +
    labs(
        x = "ATE", y = ""
    ) +
    coord_fixed(ratio = 1/350) + 
    facet_wrap(~SNP, ncol = 12) + 
    theme(legend.position = "top",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.ticks.x = element_line())

smoking_lung_cancer_individual_SNPs_plot + 
    labs(title = "Possible Trivariate Bounds from Bivariate Data",
         subtitle = paste0("Exposure: Smoking (id: ", smoking_experiment$id, "). ",
                           "Outcome: Lung Cancer (id: ", lung_cancer_experiment$id, ")."))

ggsave(
  plot = smoking_lung_cancer_individual_SNPs_plot,
  filename = here::here("figures/example_analyses", 
                        paste0("smoking_lung_cancer_individual_SNPs_plot_", smoking_experiment$id, "_", lung_cancer_experiment$id,".png")),
  height = 7, width = 10, dpi = 300
)

## -----------------------------------------------------------------------------
smoking_lung_cancer_individual_SNPs_plot_subset <- smoking_lung_cancer_trivariate_bounds %>% 
  filter(SNP %in% smoking_SNP_subset) %>% 
  group_by(SNP) %>% 
  arrange(lower) %>% 
  mutate(id = row_number()) %>% 
  ungroup() %>% 
  ggplot(aes(y = id, color = contains_zero)) +
    geom_errorbar(aes(xmin = lower,
                      xmax = upper)) +
    geom_point(data = data.frame(x = 0, y = 10,
                                 contains_zero = FALSE),
               aes(x = x, y = y),
               alpha = 0) +
    geom_vline(aes(xintercept = bivariate_lower),
               linetype = "dashed") +  
    geom_vline(aes(xintercept = bivariate_upper),
               linetype = "dashed") +
    geom_vline(xintercept = 0) +
    scale_color_manual("",
                       values = c("grey50", "red4"),
                       labels = c("Overlaps Zero", "Does not Overlap Zero"),
                       breaks = c("TRUE", "FALSE"),
                       drop = FALSE) +
    lims(x = c(-1, 1)) +
    labs(
        x = "ATE", y = ""
    ) +
    coord_fixed(ratio = 1/350) + 
    facet_wrap(~SNP, ncol = 2) +
    theme(legend.position = "top",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(),
          axis.ticks.x = element_line())

ggsave(
  plot = smoking_lung_cancer_individual_SNPs_plot_subset,
  filename = here::here("figures/example_analyses", 
                        paste0("smoking_lung_cancer_individual_SNPs_plot_subset_", smoking_experiment$id, "_", lung_cancer_experiment$id,".png")),
  width = 4, height = 6, dpi = 300
)

## ----fig.height = 8, fig.width = 13-------------------------------------------
smoking_lung_cancer_filename_mono <- paste0("samples_of_joints_mono_", smoking_experiment$id, "_", lung_cancer_experiment$id, ".Rds")

if(file.exists(here::here("vignettes_data/example_analyses", smoking_lung_cancer_filename_mono))){
  smoking_lung_cancer_samples_of_joints_mono <- read_rds(here::here("vignettes_data/example_analyses",
                                                                   smoking_lung_cancer_filename_mono))
} else {
  set.seed(6925896)
  smoking_lung_cancer_samples_of_joints_mono <- smoking_lung_cancer_bounds %>% 
    filter(!constraints_violated_mono) %>% 
    mutate(i = row_number()) %>% 
    rowwise() %>% 
    mutate(pot_covs = list(potential_covs(thetas = thetas, gammas = gammas, x_mono = TRUE)),
           sample_joint = list(sample_joint_probs(pot_covs, n = 500, return_bounds = TRUE,
                                                  x_mono = TRUE, print_progress = FALSE,
                                                  print_as_progress = NULL)))
  
  write_rds(smoking_lung_cancer_samples_of_joints_mono,
            here::here("vignettes_data/example_analyses",
                       smoking_lung_cancer_filename_mono))
}

smoking_lung_cancer_trivariate_mono_bounds <- smoking_lung_cancer_samples_of_joints_mono %>% 
  ungroup() %>% 
  select(SNP, bivariate_lower = lower_mono, bivariate_upper = upper_mono, sample_joint) %>% 
  unnest(sample_joint) %>% 
  unnest(joint) %>% 
  unnest_wider(bounds) %>% 
  mutate(contains_zero = lower < 0 & upper > 0)

smoking_lung_cancer_individual_SNPs_mono_plot <- smoking_lung_cancer_trivariate_mono_bounds %>% 
  group_by(SNP) %>% 
  arrange(lower) %>% 
  mutate(id = row_number()) %>% 
  ungroup() %>% 
  ggplot(aes(x = id)) +
    geom_hline(aes(yintercept = bivariate_lower),
               linetype = "dashed") +  
    geom_hline(aes(yintercept = bivariate_upper),
               linetype = "dashed") +
    geom_hline(yintercept = 0) +
    geom_errorbar(aes(ymin = lower,
                      ymax = upper,
                      color = contains_zero)) +
    scale_color_manual(values = c("red4", "grey50"),
                       breaks = c(FALSE, TRUE)) + 
    lims(y = c(-1, 1)) + 
    coord_fixed(ratio = 500/2) +
    facet_wrap(~SNP, nrow = 7) +
    labs(title = "Possible Trivariate Bounds from Bivariate Data assuming Monotonicity",
         subtitle = paste0("Exposure: Smoking (id: ", smoking_experiment$id, "). ",
                           "Outcome: Lung Cancer (id: ", lung_cancer_experiment$id, ")."))

smoking_lung_cancer_individual_SNPs_mono_plot

## -----------------------------------------------------------------------------
smoking_lung_cancer_pairs_of_snps <- expand_grid(SNP1 = unique(smoking_lung_cancer_trivariate_bounds$SNP),
                                                 SNP2 = unique(smoking_lung_cancer_trivariate_bounds$SNP)) %>% 
  filter(SNP1 != SNP2) %>% 
  rowwise() %>% 
  summarize(SNPs = list(sort(c(SNP1, SNP2)))) %>% 
  filter(duplicated(SNPs)) %>% 
  mutate(SNPs = map(SNPs, ~tibble_row(SNP1 = .x[1], SNP2 = .x[2]))) %>% 
  unnest_wider(SNPs)

smoking_lung_cancer_pairs_of_snps

## -----------------------------------------------------------------------------
set.seed(493637)

smoking_lung_cancer_chosen_pairs <- smoking_lung_cancer_pairs_of_snps %>% 
  sample_n(8)

pander::pander(smoking_lung_cancer_chosen_pairs)

## -----------------------------------------------------------------------------
smoking_lung_cancer_both_bounds <- smoking_lung_cancer_chosen_pairs %>% 
  rowwise() %>% 
  mutate(
    SNP1_bounds = list(
      smoking_lung_cancer_trivariate_bounds %>% 
        filter(SNP == SNP1) %>% 
        select(SNP1_bivariate_lower = bivariate_lower, SNP1_bivariate_upper = bivariate_upper, 
               SNP1_lower = lower, SNP1_upper = upper)
    ),
    SNP2_bounds = list(
      smoking_lung_cancer_trivariate_bounds %>% 
        filter(SNP == SNP2) %>% 
        select(SNP2_bivariate_lower = bivariate_lower, SNP2_bivariate_upper = bivariate_upper,
               SNP2_lower = lower, SNP2_upper = upper)
    )
  ) %>% 
  ungroup() %>% 
  unnest(cols = c(SNP1_bounds, SNP2_bounds))

smoking_lung_cancer_both_bounds

## -----------------------------------------------------------------------------
smoking_lung_cancer_intersection_bounds <- smoking_lung_cancer_both_bounds %>% 
  rowwise() %>% 
  mutate(intersection_lower = max(SNP1_lower, SNP2_lower),
         intersection_upper = min(SNP1_upper, SNP2_upper),
         bivariate_intersection_lower = max(SNP1_bivariate_lower, SNP2_bivariate_lower),
         bivariate_intersection_upper = min(SNP1_bivariate_upper, SNP2_bivariate_upper),
  ) %>% 
  ungroup() %>% 
  mutate(contains_zero = intersection_lower < 0 & intersection_upper > 0)

## ----fig.width = 10-----------------------------------------------------------
smoking_lung_cancer_intersection_bounds_plot <- smoking_lung_cancer_intersection_bounds %>% 
  group_by(SNP1, SNP2) %>% 
  arrange(intersection_lower) %>% 
  mutate(i = row_number()) %>% 
  ungroup() %>% 
  ggplot(aes(y = i)) + 
    geom_vline(aes(xintercept = bivariate_intersection_lower),
               linetype = "dashed") +
    geom_vline(aes(xintercept = bivariate_intersection_upper),
               linetype = "dashed") +
    geom_errorbar(aes(xmin = intersection_lower,
                      xmax = intersection_upper,
                      color = contains_zero),
                  alpha = 0.7) +
    geom_point(data = data.frame(x = 0, y = 10,
                                 contains_zero = FALSE),
               aes(x = x, y = y, color = contains_zero),
               alpha = 0) +
    scale_color_manual("",
                       values = c("grey50", "red4"),
                       labels = c("Overlaps Zero", "Does not Overlap Zero"),
                       breaks = c("TRUE", "FALSE"),
                       drop = FALSE) +
    lims(x = c(-1, 1)) + 
    coord_fixed(ratio = 1/350) + 
    geom_vline(xintercept = 0) +
    facet_wrap(~SNP1 + SNP2,
               nrow = 2) +
    labs(
      y = "", x = "ATE"
    ) +
    theme_bw() +
    theme(legend.position = "top",
          axis.text.x = element_text(),
          axis.ticks.x = element_line(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

ggsave(
  plot = smoking_lung_cancer_intersection_bounds_plot,
  filename = here::here("figures/example_analyses", 
                        paste0("smoking_lung_cancer_intersection_bounds_plot_", smoking_experiment$id, "_", lung_cancer_experiment$id,".png")),
  width = 8, height = 6, dpi = 300
)


smoking_lung_cancer_intersection_bounds_plot +
    labs(title = "Intersections of Possible Trivariate Bounds from Bivariate Data",
         subtitle = paste0("Exposure: Smoking (id: ", smoking_experiment$id, "). ",
                           "Outcome: Lung Cancer (id: ", lung_cancer_experiment$id, ")."))

## -----------------------------------------------------------------------------
cholesterol_experiment <- ao %>% filter(id == "ukb-a-108") 
cholesterol_instruments <- extract_instruments(cholesterol_experiment$id)

heart_attack_experiment <- ao %>% filter(id == "ukb-a-434")

heart_attack_data <- extract_outcome_data(snps = cholesterol_instruments$SNP, outcomes = heart_attack_experiment$id)

cholesterol_heart_attack_harmonized <- harmonise_data(cholesterol_instruments, heart_attack_data) %>% as_tibble()

cholesterol_heart_attack_pop_probs <- tibble(regression = c("outcome", "exposure"),
                                             p_outcome = c(cholesterol_experiment  %$% { ncase / (ncase + ncontrol) },
                                                           heart_attack_experiment %$% { ncase / (ncase + ncontrol) }))

cholesterol_heart_attack_prep_for_bounds <- cholesterol_heart_attack_harmonized %>%
  rowwise() %>% 
  mutate(samplesize.exposure = cholesterol_experiment$ncase + cholesterol_experiment$ncontrol,
         samplesize.outcome = heart_attack_experiment$ncase + heart_attack_experiment$ncontrol,
         ave_eaf = weighted.mean(x = c_across(c(eaf.exposure, eaf.outcome)),
                                 w = c_across(c(samplesize.exposure, samplesize.outcome)))) %>% 
  select(SNP, contains("beta"), ave_eaf) %>% 
  pivot_longer(contains("beta")) %>% 
  separate(name, into = c("variable", "regression"), sep = "\\.") %>% 
  pivot_wider(names_from = variable, values_from = value) %>% 
  mutate(`P(Z = 2)` = (1 - ave_eaf)^2,
         `P(Z = 1)` = 2*ave_eaf*(1 - ave_eaf),
         `P(Z = 0)` = ave_eaf^2)

cholesterol_heart_attack_w_intercept <- cholesterol_heart_attack_prep_for_bounds %>% 
  left_join(cholesterol_heart_attack_pop_probs) %>% 
  rowwise() %>% 
  mutate(find_intercept = list(find_intercept(beta1 = beta,
                                              p_outcome = p_outcome,
                                              pz = c(`P(Z = 0)`, `P(Z = 1)`, `P(Z = 2)`),
                                              interval = c(-7,7))),
         beta0 = find_intercept$root,
         `thetas/gammas` = list(1/(1+exp(-beta0 - 0:2*beta))))

## ----message = FALSE, warning = FALSE-----------------------------------------
summaries_for_plots_cholesterol_HA <- cholesterol_heart_attack_w_intercept %>% 
  ungroup() %>% 
  select(SNP, regression, beta, starts_with("P(Z ="), p_outcome, beta0, `thetas/gammas`) %>% 
  mutate(`thetas/gammas` = map(`thetas/gammas`, ~setNames(.x, c("xy_given_z0", "xy_given_z1", "xy_given_z2")))) %>% 
  unnest_wider(`thetas/gammas`)

write_csv(
  summaries_for_plots_cholesterol_HA,
  here::here("vignettes_data/summary_stats_cholesterol_on_heart_attack.csv")
)

(marginal_Z_cholesterol_HA <- summaries_for_plots_cholesterol_HA %>% 
  filter(regression == "exposure") %>% 
  pivot_longer(cols = starts_with("P(Z ="),
               names_to = "z", values_to = "p") %>% 
  ggplot(aes(x = p)) +
    geom_histogram(boundary = 0, binwidth = 0.05, 
                   color = "black") +
    facet_grid(~ z) +
    xlim(c(0,1)) + 
    theme_bw())

ggsave(
  plot = marginal_Z_cholesterol_HA,
  filename = here::here("figures/example_analyses", 
                        "cholesterol_heart_attack_marginal_Z.png"),
  dpi = 300, width = 6.5, height = 2
)

(heart_attack_coefficients <- summaries_for_plots_cholesterol_HA %>% 
  pivot_longer(cols = c(beta, beta0),
               names_to = "coef", values_to = "value") %>% 
  mutate(real_coefs = case_when(coef == "beta" & regression == "exposure" ~ "beta[1]", 
                                coef == "beta0" & regression == "exposure" ~ "beta[0]",
                                coef == "beta" & regression == "outcome" ~ "gamma[1]", 
                                coef == "beta0" & regression == "outcome" ~ "gamma[0]")) %>% 
  ggplot(aes(x = value)) +
    geom_histogram(boundary = 0, binwidth = 0.001, 
                   color = "black") +
    facet_wrap(~ real_coefs,
               scales = "free",
               labeller = label_parsed) +
    theme_bw())

ggsave(
  plot = heart_attack_coefficients,
  filename = here::here("figures/example_analyses", 
                        "cholesterol_heart_attack_coefficients.png"),
  dpi = 300, width = 5, height = 3
)

(cholesterol_heart_attack_marginal_conditionals <- summaries_for_plots_cholesterol_HA %>% 
  pivot_longer(cols = starts_with("xy_given_z"),
               names_to = "vars", values_to = "p") %>% 
  mutate(regression = if_else(regression == "exposure",
                              "P(X = 1 | Z = z)", "P(Y = 1 | Z = z)"),
         vars = paste("z =", str_extract(vars, "[0-2]"))) %>% 
  ggplot(aes(x = p)) + 
    geom_histogram(bins = 20, boundary = 0.2, 
                   color = "black") + 
    facet_grid(vars ~ regression,
               scales = "free") +
    theme_bw())

ggsave(
  plot = cholesterol_heart_attack_marginal_conditionals,
  filename = here::here("figures/example_analyses",
                        "cholesterol_heart_attack_marginal_conditionals.png"),
  dpi = 300, width = 5, height = 3
)

cholesterol_heart_attack_bounds <- cholesterol_heart_attack_w_intercept %>% 
  select(-find_intercept, -beta0, -beta, -p_outcome, -ave_eaf) %>% 
  mutate(regression = if_else(regression == "outcome", "gammas", "thetas")) %>% 
  pivot_wider(names_from = regression, values_from = `thetas/gammas`) %>% 
  rowwise() %>% 
  mutate(bivariate_bounds = list(get_bounds(thetas = thetas, gammas = gammas, 
                                            stop = FALSE, warning = FALSE)),
         interval = list(bivariate_bounds$interval),
         constraints_violated = bivariate_bounds$constraints_violated, 
         bivariate_bounds_mono = list(get_bounds(thetas = thetas, gammas = gammas, 
                                                 stop = FALSE, warning = FALSE,
                                                 x_mono = TRUE)),
         interval_mono = list(bivariate_bounds_mono$interval),
         constraints_violated_mono = bivariate_bounds_mono$constraints_violated) %>% 
  ungroup() %>% 
  unnest_wider(col = interval_mono) %>% 
  rename(upper_mono = upper, lower_mono = lower) %>% 
  unnest_wider(interval)

cholesterol_heart_attack_bivariate_bounds <- cholesterol_heart_attack_bounds %>% 
  arrange(lower) %>% 
  ggplot(aes(y = SNP)) +
    geom_errorbar(aes(xmax = upper, xmin = lower)) +
    theme_bw() +
    facet_wrap(~ lower < median(lower),
               ncol = 2, scales = "free_y") +
    theme(strip.text = element_blank(),
          strip.background = element_blank()) + 
    labs(
      y = "",
      x = "ATE"
    ) +
    lims(x = c(-1, 1))

cholesterol_heart_attack_bivariate_bounds +
      labs(title = "Bivariate bounds",
         subtitle = paste0("Exposure: High Cholesterol (id: ", cholesterol_experiment$id, "). ",
                           "Outcome: Heart Attack (id: ", heart_attack_experiment$id, ")."))

ggsave(
  plot = cholesterol_heart_attack_bivariate_bounds,
  filename = here::here("figures/example_analyses", 
                        paste0("cholesterol_heart_attack_bivariate_bounds_", cholesterol_experiment$id, "_", heart_attack_experiment$id,".png")),
  width = 8, height = 10, dpi = 300
)

cholesterol_SNP_subset <- unique(cholesterol_heart_attack_bounds$SNP)[c(2, 5, 12, 17, 25, 34, 45, 52)]

cholesterol_heart_attack_bivariate_bounds_subset <- cholesterol_heart_attack_bounds %>% 
  filter(SNP %in% cholesterol_SNP_subset) %>% 
  arrange(lower) %>% 
  ggplot(aes(y = SNP)) +
    geom_errorbar(aes(xmax = upper, xmin = lower)) +
    theme_bw() +
    theme(strip.text = element_blank(),
          strip.background = element_blank()) + 
    labs(
      y = "",
      x = "ATE"
    ) +
    lims(x = c(-1, 1))

ggsave(
  plot = cholesterol_heart_attack_bivariate_bounds_subset,
  filename = here::here("figures/example_analyses", 
                        paste0("cholesterol_heart_attack_bivariate_bounds_subset_", cholesterol_experiment$id, "_", heart_attack_experiment$id,".png")),
  width = 4, height = 4, dpi = 300
)

strength_histogram_cholesterol <- cholesterol_heart_attack_bounds %>% 
  rowwise() %>% 
  mutate(strength = max(abs(outer(thetas, thetas, `-`)))) %>% 
  ungroup() %>%
  ggplot(aes(x = strength)) + 
    geom_histogram(color = "black", binwidth = 0.001, boundary = 0) + 
    scale_x_continuous(limits = c(0, 0.012), 
                       breaks = seq(0, 0.012, by = 0.002)) +
    labs(x = "Strength", y = "Count") +
    theme_bw()

strength_histogram_cholesterol +
  labs(title = "Strength of IVs on High Cholesterol")

ggsave(
  plot = strength_histogram_cholesterol,
  filename = here::here("figures/example_analyses",
                        "cholesterol_heart_attack_strength_histogram.png"),
  dpi = 300, width = 8, height = 4
)

## -----------------------------------------------------------------------------
cholesterol_heart_attack_filename <- paste0("samples_of_joints_", cholesterol_experiment$id, "_", heart_attack_experiment$id, ".Rds")

if(file.exists(here::here("vignettes_data/example_analyses", cholesterol_heart_attack_filename))){
  cholesterol_heart_attack_samples_of_joints <- read_rds(here::here("vignettes_data/example_analyses", cholesterol_heart_attack_filename))
} else {
  set.seed(7892312)
  
  cholesterol_heart_attack_samples_of_joints <- cholesterol_heart_attack_bounds %>% 
    rowwise() %>% 
    mutate(pot_covs = list(potential_covs(thetas = thetas, gammas = gammas)),
           sample_joint = list(sample_joint_probs(pot_covs, n = 500, return_bounds = TRUE)))
  
  write_rds(cholesterol_heart_attack_samples_of_joints, here::here("vignettes_data/example_analyses", cholesterol_heart_attack_filename))
}

cholesterol_heart_attack_trivariate_bounds <- cholesterol_heart_attack_samples_of_joints %>% 
  ungroup() %>% 
  select(SNP, bivariate_lower = lower, bivariate_upper = upper, sample_joint) %>% 
  unnest(sample_joint) %>% 
  unnest(joint) %>% 
  unnest_wider(bounds) %>% 
  mutate(contains_zero = lower < 0 & upper > 0)

cholesterol_heart_attack_individual_SNPs_plot <- cholesterol_heart_attack_trivariate_bounds %>% 
  group_by(SNP) %>% 
  arrange(lower) %>% 
  mutate(id = row_number()) %>% 
  ungroup() %>% 
  ggplot(aes(y = id, color = contains_zero)) +
    geom_vline(aes(xintercept = bivariate_lower),
               linetype = "dashed") +  
    geom_vline(aes(xintercept = bivariate_upper),
               linetype = "dashed") +
    geom_vline(xintercept = 0) +
    geom_errorbar(aes(xmin = lower,
                      xmax = upper)) +
    geom_point(data = data.frame(x = 0, y = 10,
                                 contains_zero = FALSE),
               aes(x = x, y = y),
               alpha = 0) +
    scale_color_manual("",
                       values = c("grey50", "red4"),
                       labels = c("Overlaps Zero", "Does not Overlap Zero"),
                       breaks = c("TRUE", "FALSE"),
                       drop = FALSE) +
    lims(x = c(-1, 1)) +
    labs(
        x = "ATE", y = ""
    ) +
    coord_fixed(ratio = 1/350) + 
    facet_wrap(~SNP, ncol = 12) +
    theme(legend.position = "top",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.ticks.x = element_line())

cholesterol_heart_attack_individual_SNPs_plot + 
    labs(title = "Possible Trivariate Bounds from Bivariate Data",
         subtitle = paste0("Exposure: High Cholesterol (id: ", cholesterol_experiment$id, "). ",
                           "Outcome: Heart Attack (id: ", heart_attack_experiment$id, ")."))

ggsave(
  plot = cholesterol_heart_attack_individual_SNPs_plot,
  filename = here::here("figures/example_analyses", 
                        paste0("cholesterol_heart_attack_individual_SNPs_plot_", cholesterol_experiment$id, "_", heart_attack_experiment$id,".png")),
  height = 5.5, width = 10, dpi = 300
)

cholesterol_heart_attack_individual_SNPs_plot_subset <- cholesterol_heart_attack_trivariate_bounds %>% 
  filter(SNP %in% cholesterol_SNP_subset) %>% 
  group_by(SNP) %>% 
  arrange(lower) %>% 
  mutate(id = row_number()) %>% 
  ungroup() %>% 
  ggplot(aes(y = id, color = contains_zero)) +
    geom_errorbar(aes(xmin = lower,
                      xmax = upper)) +
    geom_point(data = data.frame(x = 0, y = 10,
                                 contains_zero = FALSE),
               aes(x = x, y = y),
               alpha = 0) +
    geom_vline(aes(xintercept = bivariate_lower),
               linetype = "dashed") +  
    geom_vline(aes(xintercept = bivariate_upper),
               linetype = "dashed") +
    geom_vline(xintercept = 0) +
    scale_color_manual("",
                       values = c("grey50", "red4"),
                       labels = c("Overlaps Zero", "Does not Overlap Zero"),
                       breaks = c("TRUE", "FALSE"),
                       drop = FALSE) +
    lims(x = c(-1, 1)) +
    labs(
        x = "ATE", y = ""
    ) +
    coord_fixed(ratio = 1/350) + 
    facet_wrap(~SNP, ncol = 2) +
    theme(legend.position = "top",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(),
          axis.ticks.x = element_line())

ggsave(
  plot = cholesterol_heart_attack_individual_SNPs_plot_subset,
  filename = here::here("figures/example_analyses", 
                        paste0("cholesterol_heart_attack_individual_SNPs_plot_subset_", cholesterol_experiment$id, "_", heart_attack_experiment$id,".png")),
  width = 4, height = 6, dpi = 300
)

cholesterol_heart_attack_filename_mono <- paste0("samples_of_joints_mono_", smoking_experiment$id, "_", lung_cancer_experiment$id, ".Rds")

if(file.exists(here::here("vignettes_data/example_analyses", cholesterol_heart_attack_filename_mono))){
  cholesterol_heart_attack_samples_of_joints_mono <- read_rds(here::here("vignettes_data/example_analyses",
                                                                   cholesterol_heart_attack_filename_mono))
} else {
  set.seed(6925896)
  cholesterol_heart_attack_samples_of_joints_mono <- cholesterol_heart_attack_bounds %>% 
    filter(!constraints_violated_mono) %>% 
    mutate(i = row_number()) %>% 
    rowwise() %>% 
    mutate(pot_covs = list(potential_covs(thetas = thetas, gammas = gammas, x_mono = TRUE)),
           sample_joint = list(sample_joint_probs(pot_covs, n = 500, return_bounds = TRUE,
                                                  x_mono = TRUE, print_progress = FALSE,
                                                  print_as_progress = NULL)))
  
  write_rds(cholesterol_heart_attack_samples_of_joints_mono,
            here::here("vignettes_data/example_analyses",
                       cholesterol_heart_attack_filename_mono))
}

cholesterol_heart_attack_trivariate_mono_bounds <- cholesterol_heart_attack_samples_of_joints_mono %>% 
  ungroup() %>% 
  select(SNP, bivariate_lower = lower_mono, bivariate_upper = upper_mono, sample_joint) %>% 
  unnest(sample_joint) %>% 
  unnest(joint) %>% 
  unnest_wider(bounds) %>% 
  mutate(contains_zero = lower < 0 & upper > 0)

cholesterol_heart_attack_individual_SNPs_mono_plot <- cholesterol_heart_attack_trivariate_mono_bounds %>% 
  group_by(SNP) %>% 
  arrange(lower) %>% 
  mutate(id = row_number()) %>% 
  ungroup() %>% 
  ggplot(aes(x = id)) +
    geom_hline(aes(yintercept = bivariate_lower),
               linetype = "dashed") +  
    geom_hline(aes(yintercept = bivariate_upper),
               linetype = "dashed") +
    geom_hline(yintercept = 0) +
    geom_errorbar(aes(ymin = lower,
                      ymax = upper,
                      color = contains_zero)) +
    scale_color_manual(values = c("red4", "grey50"),
                       breaks = c(FALSE, TRUE)) + 
    lims(y = c(-1, 1)) + 
    coord_fixed(ratio = 500/2) +
    facet_wrap(~SNP, nrow = 7) +
    labs(title = "Possible Trivariate Bounds from Bivariate Data assuming Monotonicity",
         subtitle = paste0("Exposure: High Cholesterol (id: ", cholesterol_experiment$id, "). ",
                           "Outcome: Heart Attack (id: ", heart_attack_experiment$id, ")."))

cholesterol_heart_attack_individual_SNPs_mono_plot

## -----------------------------------------------------------------------------
cholesterol_heart_attack_pairs_of_snps <- expand_grid(SNP1 = unique(cholesterol_heart_attack_trivariate_bounds$SNP),
                                                      SNP2 = unique(cholesterol_heart_attack_trivariate_bounds$SNP)) %>% 
  filter(SNP1 != SNP2) %>% 
  rowwise() %>% 
  summarize(SNPs = list(sort(c(SNP1, SNP2)))) %>% 
  filter(duplicated(SNPs)) %>% 
  mutate(SNPs = map(SNPs, ~tibble_row(SNP1 = .x[1], SNP2 = .x[2]))) %>% 
  unnest_wider(SNPs)

cholesterol_heart_attack_pairs_of_snps

set.seed(493637)

cholesterol_heart_attack_chosen_pairs <- cholesterol_heart_attack_pairs_of_snps %>% 
  sample_n(8)

pander::pander(cholesterol_heart_attack_chosen_pairs)

cholesterol_heart_attack_both_bounds <- cholesterol_heart_attack_chosen_pairs %>% 
  rowwise() %>% 
  mutate(
    SNP1_bounds = list(
      cholesterol_heart_attack_trivariate_bounds %>% 
        filter(SNP == SNP1) %>% 
        select(SNP1_bivariate_lower = bivariate_lower, SNP1_bivariate_upper = bivariate_upper, 
               SNP1_lower = lower, SNP1_upper = upper)
    ),
    SNP2_bounds = list(
      cholesterol_heart_attack_trivariate_bounds %>% 
        filter(SNP == SNP2) %>% 
        select(SNP2_bivariate_lower = bivariate_lower, SNP2_bivariate_upper = bivariate_upper,
               SNP2_lower = lower, SNP2_upper = upper)
    )
  ) %>% 
  ungroup() %>% 
  unnest(cols = c(SNP1_bounds, SNP2_bounds))

cholesterol_heart_attack_both_bounds

cholesterol_heart_attack_intersection_bounds <- cholesterol_heart_attack_both_bounds %>% 
  rowwise() %>% 
  mutate(intersection_lower = max(SNP1_lower, SNP2_lower),
         intersection_upper = min(SNP1_upper, SNP2_upper),
         bivariate_intersection_lower = max(SNP1_bivariate_lower, SNP2_bivariate_lower),
         bivariate_intersection_upper = min(SNP1_bivariate_upper, SNP2_bivariate_upper),
  ) %>% 
  ungroup() %>% 
  mutate(contains_zero = intersection_lower < 0 & intersection_upper > 0)

cholesterol_heart_attack_intersection_bounds_plot <- cholesterol_heart_attack_intersection_bounds %>% 
  group_by(SNP1, SNP2) %>% 
  arrange(intersection_lower) %>% 
  mutate(i = row_number()) %>% 
  ungroup() %>% 
  ggplot(aes(y = i)) + 
    geom_vline(aes(xintercept = bivariate_intersection_lower),
               linetype = "dashed") +
    geom_vline(aes(xintercept = bivariate_intersection_upper),
               linetype = "dashed") +
    geom_errorbar(aes(xmin = intersection_lower,
                      xmax = intersection_upper,
                      color = contains_zero),
                  alpha = 0.7) +
    geom_point(data = data.frame(x = 0, y = 10,
                                 contains_zero = FALSE),
               aes(x = x, y = y, color = contains_zero),
               alpha = 0) +
    scale_color_manual("",
                       values = c("grey50", "red4"),
                       labels = c("Overlaps Zero", "Does not Overlap Zero"),
                       breaks = c("TRUE", "FALSE"),
                       drop = FALSE) +
    lims(x = c(-1, 1)) + 
    coord_fixed(ratio = 1/350) + 
    geom_vline(xintercept = 0) +
    facet_wrap(~SNP1 + SNP2,
               nrow = 2) +
    labs(
      y = "", x = "ATE"
    ) +
    theme_bw() +
    theme(legend.position = "top",
          axis.text.x = element_text(),
          axis.ticks.x = element_line(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

ggsave(
  plot = cholesterol_heart_attack_intersection_bounds_plot,
  filename = here::here("figures/example_analyses", 
                        paste0("cholesterol_heart_attack_intersection_bounds_plot_", cholesterol_experiment$id, "_", heart_attack_experiment$id,".png")),
  width = 8, height = 6, dpi = 300
)


cholesterol_heart_attack_intersection_bounds_plot +
    labs(title = "Intersections of Possible Trivariate Bounds from Bivariate Data",
         subtitle = paste0("Exposure: High Cholesterol (id: ", cholesterol_experiment$id, "). ",
                           "Outcome: Heart Attack (id: ", heart_attack_experiment$id, ")."))

## -----------------------------------------------------------------------------
write_csv(x = 
            bind_rows(
              smoking_experiment,
              lung_cancer_experiment,
              cholesterol_experiment,
              heart_attack_experiment
            ),
          file = here::here("vignettes_data/example_analyses/experiment_info.csv")
)

