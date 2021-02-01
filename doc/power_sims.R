## ---- code = readLines(here::here("scripts/power/power_sims_only.R")), eval = FALSE----
#  ## Run for sim_n = 1,2,3,...,840.
#  sim_n <- as.numeric(commandArgs(trailingOnly = TRUE)[[1]])
#  
#  ## Only run if results don't already exist.
#  if(!file.exists(here::here("data/power_sims",
#                             paste0("power_sims_", stringr::str_pad(sim_n, width = 3, side = "left", pad = 0), ".rds")))){
#  
#    ## Make sure V8 is available. If not, install in special folder.
#    ## This is mainly a hack to make this work on our department's server, and can
#    ## probably be removed.
#    if(!require(V8)){
#      Sys.setenv(DOWNLOAD_STATIC_LIBV8 = 1)
#      install.packages("V8", repos = "https://cloud.r-project.org", lib = paste("V8", sim_n, sep = "_"))
#    }
#  
#    library(ACEBounds)
#    library(distributions3)
#    library(tidyverse)
#  
#    set.seed(9866311)
#  
#    ## Create tibble with all combinations of coefficients, ids, and
#    ## seed to use for each combination.
#    all_combos <- expand_grid(indIVs_on_X = seq(0.2, 6, by = 0.2),
#                              X_on_Y = c(0.25, 0.5, 1, 1.5, 2, 4, 6),
#                              U_on_XY = c(0.1, 0.5, 1, 2)) %>%
#      mutate(i = row_number(),
#             seed = round(round(runif(n = n()), digits = 7)*1e7, digits = 0))
#  
#    ## Get the combination we are working on.
#    this_combo <- filter(all_combos, i == sim_n)
#  
#    set.seed(this_combo$seed)
#  
#    ## Simulate data
#    sim_data <- simulate_data(
#      sample_size = 1e7,
#      IVs_ps = list(c(0.25, 0.5, 0.25)),
#      X_intercept = -this_combo$indIVs_on_X,
#      Y_intercept = -this_combo$X_on_Y/2,
#      indIVs_on_X = this_combo$indIVs_on_X,
#      indIVs_on_Y = 0,
#      U = distributions3::Normal(),
#      U_on_X = this_combo$U_on_XY,
#      U_on_Y = this_combo$U_on_XY,
#      X_on_Y = this_combo$X_on_Y
#    )
#  
#    ## Add simulated data to tibble and save to .Rds file.
#    write_rds(this_combo %>% mutate(sim_data = list(sim_data)),
#              file = here::here("data/power_sims",
#                                paste0("power_sims_", str_pad(sim_n, width = 3, side = "left", pad = 0), ".rds")),
#              compress = "gz",
#              compression = 9L)
#  }

## ---- code = readLines(here::here("scripts/power/power_prep_for_plots.R")), eval = FALSE----
#  ## Specify number of cores for parallelization and run_n. Must be run for run_n = 1,2,3,...,84.
#  n_cores <- as.numeric(commandArgs(trailingOnly = TRUE)[[1]])
#  run_n <- as.numeric(commandArgs(trailingOnly = TRUE)[[2]])
#  
#  
#  ## Make sure V8 is available. If not, install in special folder.
#  ## This is mainly a hack to make this work on our department's server, and can
#  ## probably be removed.
#  if(!require(V8)){
#    Sys.setenv(DOWNLOAD_STATIC_LIBV8 = 1)
#    install.packages("V8", repos = "https://cloud.r-project.org", lib = paste("V8", run_n, sep = "_"))
#  }
#  
#  library(tidyverse)
#  library(ACEBounds)
#  library(distributions3)
#  library(furrr)
#  library(progressr)
#  
#  
#  ## Helper function to find ATE from simulated data.
#  ATE_from_simulated_data <- function(from_simulate_data){
#    intercept <- filter(from_simulate_data$coefficients, effect == "Yintercept")$coef
#    x_beta <- filter(from_simulate_data$coefficients, effect == "X_on_Y")$coef
#    u_beta <- filter(from_simulate_data$coefficients, effect == "U_on_Y")$coef
#  
#    pY1X0 <- 1 / (1 + exp(-intercept - u_beta*from_simulate_data$simulated_data$U))
#    pY1X1 <- 1 / (1 + exp(-intercept - x_beta - u_beta*from_simulate_data$simulated_data$U))
#  
#    return(mean(pY1X1 - pY1X0))
#  }
#  
#  ## Get files with results.
#  all_power_sims <- list.files(here::here("data/power_sims"), full.names = TRUE)
#  
#  plan(multisession,
#       workers = n_cores-1)
#  
#  ## Which chunk of 10 are we working on?
#  to_do <- tibble(data_file = all_power_sims) %>%
#    mutate(j = ceiling(row_number()/10)) %>%
#    filter(j == run_n)
#  
#  with_progress({
#    p <- progressor(steps = nrow(to_do))
#  
#    ## For each file...
#    bounds_and_ATE <- to_do %>%
#      mutate(
#        subset = future_map(
#          data_file,
#          function(x){
#            ## ... read in the file
#            tmp <- read_rds(x)
#  
#            ## ... find summary statistics and ATE
#            out <- tmp %>%
#              mutate(
#                sum_stats = map(sim_data, ~c(ACEBounds::probs_from_data(.x$simulated_data, X, Y, Z1, data_format = "bivariate"),
#                                             ATE = ATE_from_simulated_data(.x)))
#              ) %>%
#              select(-sim_data) %>%
#              ## ... find bounds
#              mutate(
#                get_bounds_res = map(sum_stats, ~ACEBounds::get_bounds(gammas = .x$gammas, thetas = .x$thetas, stop = FALSE, warning = FALSE)),
#                bounds = map(get_bounds_res, "interval")
#              )
#  
#            p()
#  
#            return(out)
#          }
#        )
#      )
#  })
#  
#  ## Save to file
#  write_rds(bounds_and_ATE, here::here(paste0("data/power_results/power_bounds_and_ATE_", run_n, ".Rds")))

## ---- code = readLines(here::here("scripts/power/power_combine_results.R")), eval = FALSE----
#  library(tidyverse)
#  
#  map_dfr(list.files(here::here("data/power_results"),
#                     pattern = "power_bounds_and_ATE_[0-9]+.Rds",
#                     full.names = TRUE),
#          read_rds) %>%
#      unnest(subset) %>%
#      write_rds(file = here::here("data/power_combined_results.Rds"))
#  
#  
#  
#  
#  

