n_nodes <- as.numeric(commandArgs(trailingOnly = TRUE)[[1]])
slurm_array_id <- as.numeric(commandArgs(trailingOnly = TRUE)[[2]])

if(!require(V8)){
  Sys.setenv(DOWNLOAD_STATIC_LIBV8 = 1)
  install.packages("V8", repos = "https://cloud.r-project.org", lib = paste("V8", slurm_array_id, sep = "_"))
}

library(tidyverse)
library(ACEBounds)
library(distributions3)
library(furrr)
library(progressr)

ATE_from_simulated_data <- function(from_simulate_data){
  intercept <- filter(from_simulate_data$coefficients, effect == "Yintercept")$coef
  x_beta <- filter(from_simulate_data$coefficients, effect == "X_on_Y")$coef
  u_beta <- filter(from_simulate_data$coefficients, effect == "U_on_Y")$coef

  pY1X0 <- 1 / (1 + exp(-intercept - u_beta*from_simulate_data$simulated_data$U))
  pY1X1 <- 1 / (1 + exp(-intercept - x_beta - u_beta*from_simulate_data$simulated_data$U))

  return(mean(pY1X1 - pY1X0))
}

all_power_sims <- list.files(here::here("data/power_sims"), full.names = TRUE)

plan(multisession,
     workers = n_nodes-1)

to_do <- tibble(data_file = all_power_sims) %>%
  mutate(j = ceiling(row_number()/10)) %>%
  filter(j == slurm_array_id)

with_progress({
  p <- progressor(steps = nrow(to_do))
  bounds_and_ATE <- to_do %>%
    mutate(
      subset = future_map(
        data_file,
        function(x){
          tmp <- read_rds(x)

          out <- tmp %>%
            mutate(
              sum_stats = map(sim_data, ~c(ACEBounds::probs_from_data(.x$simulated_data, X, Y, Z1, data_format = "bivariate"),
                                           ATE = ATE_from_simulated_data(.x)))
            ) %>%
            select(-sim_data) %>%
            mutate(
              get_bounds_res = map(sum_stats, ~ACEBounds::get_bounds(gammas = .x$gammas, thetas = .x$thetas, stop = FALSE, warning = FALSE)),
              bounds = map(get_bounds_res, "interval")
            )

          p()

          return(out)
        }
      )
    )
})

write_rds(bounds_and_ATE, here::here(paste0("data/power_results/power_bounds_and_ATE_", slurm_array_id, ".Rds")))
