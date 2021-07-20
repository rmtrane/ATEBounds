## Specify number of cores for parallelization and run_n. Must be run for run_n = 1,2,3,...,84.
n_cores <- as.numeric(commandArgs(trailingOnly = TRUE)[[1]])
run_n <- as.numeric(commandArgs(trailingOnly = TRUE)[[2]])


## Make sure V8 is available. If not, install in special folder.
## This is mainly a hack to make this work on our department's server, and can
## probably be removed.
if(!require(V8)){
  Sys.setenv(DOWNLOAD_STATIC_LIBV8 = 1)
  install.packages("V8", repos = "https://cloud.r-project.org", lib = paste("V8", run_n, sep = "_"))
}

library(tidyverse)
library(ATEBounds)
library(distributions3)
library(furrr)
library(progressr)


## Helper function to find ATE from simulated data.
ATE_from_simulated_data <- function(from_simulate_data){
  intercept <- filter(from_simulate_data$coefficients, effect == "Yintercept")$coef
  x_beta <- filter(from_simulate_data$coefficients, effect == "X_on_Y")$coef
  u_beta <- filter(from_simulate_data$coefficients, effect == "U_on_Y")$coef

  pY1X0 <- 1 / (1 + exp(-intercept - u_beta*from_simulate_data$simulated_data$U))
  pY1X1 <- 1 / (1 + exp(-intercept - x_beta - u_beta*from_simulate_data$simulated_data$U))

  return(mean(pY1X1 - pY1X0))
}

## Get files with results.
all_power_sims <- list.files(here::here("data/power_sims"), full.names = TRUE)

plan(multisession,
     workers = n_cores-1)

## Which chunk of 10 are we working on?
to_do <- tibble(data_file = all_power_sims) %>%
  mutate(j = ceiling(row_number()/10)) %>%
  filter(j == run_n)

with_progress({
  p <- progressor(steps = nrow(to_do))

  ## For each file...
  bounds_and_ATE <- to_do %>%
    mutate(
      subset = future_map(
        data_file,
        function(x){
          ## ... read in the file
          tmp <- read_rds(x)

          ## ... find summary statistics and ATE
          out <- tmp %>%
            mutate(
              sum_stats = map(sim_data, ~c(ATEBounds::probs_from_data(.x$simulated_data, X, Y, Z1, data_format = "bivariate"),
                                           ATE = ATE_from_simulated_data(.x)))
            ) %>%
            select(-sim_data) %>%
            ## ... find bounds
            mutate(
              get_bounds_res = map(sum_stats, ~ATEBounds::get_bounds(gammas = .x$gammas, thetas = .x$thetas, stop = FALSE, warning = FALSE)),
              bounds = map(get_bounds_res, "interval")
            )

          p()

          return(out)
        }
      )
    )
})

## Save to file
write_rds(bounds_and_ATE, here::here(paste0("data/power_results/power_bounds_and_ATE_", run_n, ".Rds")))
