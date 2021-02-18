library(tidyverse)
library(patchwork)

pip_fig <- read_rds(here::here("vignettes_data/pip_figure.Rds"))
coefs_vs_strength <- read_rds(here::here("vignettes_data/coefs_vs_strength.Rds"))

comb_figure <- pip_fig + coefs_vs_strength +
  plot_annotation(tag_levels = "A")

ggsave(comb_figure,
       filename = here::here("figures/png/pip_and_coefs_vs_strength.png"),
       width = 8, height = 4.7, dpi = 300)

ggsave(comb_figure,
       filename = here::here("figures/tiff/pip_and_coefs_vs_strength.tiff"),
       width = 8, height = 4.7, dpi = 300)
