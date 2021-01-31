## Run for sim_n = 1,2,3,...,840.
sim_n <- as.numeric(commandArgs(trailingOnly = TRUE)[[1]])

## Only run if results don't already exist.
if(!file.exists(here::here("data/power_sims",
                           paste0("power_sims_", stringr::str_pad(sim_n, width = 3, side = "left", pad = 0), ".rds")))){

  ## Make sure V8 is available. If not, install in special folder.
  ## This is mainly a hack to make this work on our department's server, and can
  ## probably be removed.
  if(!require(V8)){
    Sys.setenv(DOWNLOAD_STATIC_LIBV8 = 1)
    install.packages("V8", repos = "https://cloud.r-project.org", lib = paste("V8", sim_n, sep = "_"))
  }

  library(ACEBounds)
  library(distributions3)
  library(tidyverse)

  set.seed(9866311)

  ## Create tibble with all combinations of coefficients, ids, and
  ## seed to use for each combination.
  all_combos <- expand_grid(indIVs_on_X = seq(0.2, 6, by = 0.2),
                            X_on_Y = c(0.25, 0.5, 1, 1.5, 2, 4, 6),
                            U_on_XY = c(0.1, 0.5, 1, 2)) %>%
    mutate(i = row_number(),
           seed = round(round(runif(n = n()), digits = 7)*1e7, digits = 0))

  ## Get the combination we are working on.
  this_combo <- filter(all_combos, i == sim_n)

  set.seed(this_combo$seed)

  ## Simulate data
  sim_data <- simulate_data(
    sample_size = 1e7,
    IVs_ps = list(c(0.25, 0.5, 0.25)),
    X_intercept = -this_combo$indIVs_on_X,
    Y_intercept = -this_combo$X_on_Y/2,
    indIVs_on_X = this_combo$indIVs_on_X,
    indIVs_on_Y = 0,
    U = distributions3::Normal(),
    U_on_X = this_combo$U_on_XY,
    U_on_Y = this_combo$U_on_XY,
    X_on_Y = this_combo$X_on_Y
  )

  ## Add simulated data to tibble and save to .Rds file.
  write_rds(this_combo %>% mutate(sim_data = list(sim_data)),
            file = here::here("data/power_sims",
                              paste0("power_sims_", str_pad(sim_n, width = 3, side = "left", pad = 0), ".rds")),
            compress = "gz",
            compression = 9L)
}
