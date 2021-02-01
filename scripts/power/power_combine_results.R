library(tidyverse)

map_dfr(list.files(here::here("data/power_results"),
                   pattern = "power_bounds_and_ATE_[0-9]+.Rds",
                   full.names = TRUE),
        read_rds) %>%
    unnest(subset) %>%
    write_rds(file = here::here("data/power_results/power_combined_results.Rds"))





