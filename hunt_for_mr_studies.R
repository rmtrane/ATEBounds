library(tidyverse)
library(TwoSampleMR)
library(magrittr)

find_intercept <- function(beta1, p_outcome, pz, interval = c(-5, 5), zs = 0:2){
  uniroot(f = function(beta0) sum((1+exp(-beta0-beta1*zs))^(-1) * pz) - p_outcome,
          interval = interval)
}

ao <- available_outcomes()

## Tried and failed:
## * obesity (ieu-a-90) and stroke (ieu-a-1109)
## * obesity (ieu-a-90) and heart disease (ukb-d-I9_IHD)
## * cholesterol (ukb-a-108) and heart attack (ukb-b-11590)
## * cholesterol (ukb-a-108) and heart attack (ukb-a-434)
## * cholesterol (ukb-a-108) and obesity (ieu-a-90)
## * obesity (ieu-a-90) and diabetes (ukb-a-306)
## * lipo (ukb-b-17462) and heart attack (ukb-a-434)
## * hypercholesterolaemia (ukb-b-12651) and heart attack (ukb-a-434)
## * cholesterol (ukb-a-108) and stroke (ieu-a-1109)
## * obesity (ieu-a-90) and diabetes (ieu-a-26)
## * cholesterol (ukb-a-108) and heart disease (ukb-d-I9_IHD)

ao %>% filter(str_detect(tolower(trait), "heart attack")) %>%
  select(id, trait, ncase, ncontrol, sample_size) %>% filter(complete.cases(.)) %>%
  print(n = 60)

lipo_id <- "ukb-b-17462" # max(beta) = 0.0015
stroke_id <- "ieu-a-1109" # max(beta) = 0.2198
cholesterol <- "ukb-a-108" # max(beta) = 0.0312
cholesterol_2 <- "ukb-a-488"
cholesterol_3 <- "ukb-b-10912" # max(beta) = 0.0311
hypercholesterolaemia <- "ukb-b-12651" # max(beta) = 0.0185

obesity <- "ieu-a-90" # max(beta) = 0.21
heart_disease <- "ukb-d-I9_IHD" # max(beta) = 0.0245
heart_attack <- "ukb-a-434" # max(beta) = 0.01169
heart_attack_2 <- "ukb-b-15829" # 0.0076
heart_attack_3 <- "ukb-b-11590" # 0.0076

diabetes <- "ukb-a-306" # max(beta) = 0.015
diabetes_2 <- "ieu-a-26" # max(beta) = 0.122
diabetes_3 <- "ukb-b-10753" # max(beta) 0.005
diabetes_4 <- "bbj-a-77" # max(beta) = 0.117

high_blood_pressure <- "ukb-a-437" # max(beta) = 0.0292887
overweight <- "ieu-a-93" # 0.14
vitamin_D <- "ukb-b-4991" #



########
exposure_id <- obesity
outcome_id <- stroke_id

{
  exposure_experiment <- ao %>% filter(id == exposure_id)
  exposure_instruments <- extract_instruments(outcomes = exposure_id)
  exposure_instruments %>% pull(beta.exposure) %>% max() %>% print()

  outcome_experiment <- ao %>% filter(id == outcome_id)
  outcome_data <- extract_outcome_data(snps = exposure_instruments$SNP, outcomes = outcome_experiment$id)
  outcome_data %>% pull(beta.outcome) %>% max() %>% print()
}

{
  harmonized_data <- harmonise_data(exposure_instruments, outcome_data) # %>% as_tibble()

  pop_probs <- tibble(regression = c("outcome", "exposure"),
                      p_outcome = c(outcome_experiment  %$% { ncase / (ncase + ncontrol) },
                                    exposure_experiment %$% { ncase / (ncase + ncontrol) }))

  prep_for_bounds <- harmonized_data %>%
    rowwise() %>%
    mutate(samplesize.exposure = exposure_experiment$ncase + exposure_experiment$ncontrol,
           samplesize.outcome = outcome_experiment$ncase + outcome_experiment$ncontrol,
           ave_eaf = weighted.mean(x = c_across(c(eaf.exposure, eaf.outcome)),
                                   w = c_across(c(samplesize.exposure, samplesize.outcome)),
                                   na.rm = TRUE)) %>%
    select(SNP, contains("beta"), ave_eaf) %>%
    pivot_longer(contains("beta")) %>%
    separate(name, into = c("variable", "regression"), sep = "\\.") %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    mutate(`P(Z = 2)` = (1 - ave_eaf)^2,
           `P(Z = 1)` = 2*ave_eaf*(1 - ave_eaf),
           `P(Z = 0)` = ave_eaf^2)

  include_intercept <- prep_for_bounds %>%
    left_join(pop_probs) %>%
    rowwise() %>%
    mutate(find_intercept = list(find_intercept(beta1 = beta,
                                                p_outcome = p_outcome,
                                                pz = c(`P(Z = 0)`, `P(Z = 1)`, `P(Z = 2)`),
                                                interval = c(-10, 10))),
           beta0 = find_intercept$root,
           `thetas/gammas` = list(1/(1+exp(-beta0 - 0:2*beta))))

  bivariate_bounds <- include_intercept %>%
    select(-find_intercept, -beta0, -beta, -p_outcome, -ave_eaf) %>%
    mutate(regression = if_else(regression == "outcome", "gammas", "thetas")) %>%
    pivot_wider(names_from = regression, values_from = `thetas/gammas`) %>%
    rowwise() %>%
    mutate(bivariate_bounds = list(get_bounds(thetas = thetas, gammas = gammas,
                                              stop = FALSE, warning = FALSE)),
           interval = list(bivariate_bounds$interval),
           constraints_violated = bivariate_bounds$constraints_violated) %>%
    ungroup() %>%
    unnest_wider(interval)

  bivariate_bounds %>%
    arrange(lower) %>%
    ggplot(aes(x = SNP)) +
      geom_errorbar(aes(ymax = upper, ymin = lower)) +
      theme_bw() +
      facet_wrap(~ lower < median(lower),
                 ncol = 1, scales = "free_x") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text = element_blank(),
            strip.background = element_blank()) +
      labs(
        x = "SNP id",
        y = "ATE"
      ) +
      lims(y = c(-1, 1))

# bivariate_bounds_plot +
#   labs(title = "Bivariate bounds",
#        subtitle = paste0("Exposure: ", exposure_experiment, "). ",
#                          "Outcome: Depression (id: ", depression_experiment$id, ")."))
#
# ggsave(
#   plot = smoking_depression_bivariate_bounds,
#   filename = here::here("figures/example_analyses",
#                         paste0("smoking_depression_bivaraite_bounds_", smoking_experiment$id, "_", depression_experiment$id,".png")),
#   width = 8, height = 6, dpi = 300
# )

bivariate_bounds %>%
  rowwise() %>%
  mutate(strength = max(abs(outer(thetas, thetas, `-`)))) %>%
  ungroup() %>%
    ggplot(aes(x = strength, y = upper - lower)) +
    geom_point() +
    theme_bw() +
    geom_abline(slope = -2, intercept = 2) +
    lims(y = c(0, 2), x= c(0, 1))

# obesity_strength_histogram +
#   labs(title = "Strength of IVs on Obesity")
#
# ggsave(
#   plot = strength_histogram,
#   filename = here::here("figures/example_analyses",
#                         "strength_histogram.png"),
#   dpi = 300, width = 8, height = 4
# )

## Random Sampling of Potential Joint Distributions

filename <- paste0("samples_of_joints_", exposure_experiment$id, "_", outcome_experiment$id, ".Rds")

set.seed(3843173)

samples_of_joints <- bivariate_bounds %>%
  ungroup() %>%
  # filter(abs(upper - lower) %in% range(abs(upper - lower))) %>%
  mutate(j = paste(row_number(), "of", n())) %>%
  rowwise() %>%
  mutate(pot_covs = list(potential_covs(thetas = thetas, gammas = gammas)),
         sample_joint = list(sample_joint_probs(pot_covs, n = 200, return_bounds = TRUE, print_progress = FALSE, print_as_progress = j))) %>%
  ungroup()

# write_rds(samples_of_joints, here::here("vignettes_data/example_analyses", filename))


trivariate_bounds <- samples_of_joints %>%
  select(SNP, bivariate_lower = lower, bivariate_upper = upper, sample_joint) %>%
  unnest(sample_joint) %>%
  unnest(joint) %>%
  unnest_wider(bounds) %>%
  mutate(contains_zero = lower < 0 & upper > 0)

trivariate_bounds_plot <- trivariate_bounds %>%
  group_by(SNP) %>%
  arrange(lower) %>%
  mutate(id = row_number()) %>%
  ungroup() %>%
  ggplot(aes(x = id, color = contains_zero)) +
  geom_hline(aes(yintercept = bivariate_lower),
             linetype = "dashed") +
  geom_hline(aes(yintercept = bivariate_upper),
             linetype = "dashed") +
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymin = lower,
                    ymax = upper)) +
  geom_point(data = data.frame(x = 10, y = 0,
                               contains_zero = FALSE),
             aes(x = x, y = y),
             alpha = 0) +
  scale_color_manual("",
                     values = c("grey50", "red4"),
                     labels = c("Overlaps Zero", "Does not Overlap Zero"),
                     breaks = c("TRUE", "FALSE"),
                     drop = FALSE) +
  lims(y = c(-1, 1)) +
  labs(
    x = "", y = "ATE"
  ) +
  #coord_fixed(ratio = 350/2) +
  facet_wrap(~SNP, ncol = 12) +
  theme(legend.position = "top")

trivariate_bounds_plot
}
