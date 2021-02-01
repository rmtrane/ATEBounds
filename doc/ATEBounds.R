## ----message = FALSE, warning = FALSE-----------------------------------------
set.seed(184506)

library(ATEBounds)
library(dplyr)

sim_dat <- simulate_data(sample_size = 500,
                         IVs_ps = list(c(1, 1, 1)/3),
                         indIVs_on_X = 0.7)

library(ggdag)
plot_DAG(sim_dat) + scale_x_reverse()

## -----------------------------------------------------------------------------
knitr::kable(head(sim_dat$simulated_data, n = 10))

## -----------------------------------------------------------------------------
one_sample_data <- sim_dat$simulated_data

two_sample_data_X <- sim_dat$simulated_data %>% select(X,Z1)
two_sample_data_Y <- sim_dat$simulated_data %>% select(Y,Z1)

## -----------------------------------------------------------------------------
probs <- probs_from_data(one_sample_data, X = X, Y = Y, Z = Z1)

get_bounds(zetas = probs$zetas)

## -----------------------------------------------------------------------------
two_sample_probs <- probs_from_data(one_sample_data, X = X, Y = Y, Z = Z1, data_format = "bivariate")

get_bounds(thetas = two_sample_probs$thetas,
           gammas = two_sample_probs$gammas)

## -----------------------------------------------------------------------------
X_probs <- probs_from_data(two_sample_data_X, X = X, Z = Z1)
Y_probs <- probs_from_data(two_sample_data_Y, Y = Y, Z = Z1)

get_bounds(thetas = X_probs$thetas,
           gammas = Y_probs$gammas)

## -----------------------------------------------------------------------------
get_bounds(zetas = probs$zetas,
           x_mono = TRUE)

get_bounds(gammas = Y_probs$gammas,
           thetas = X_probs$thetas,
           x_mono = TRUE)

## -----------------------------------------------------------------------------
get_bounds(zetas = probs$zetas,
           x_mono = TRUE,
           y_mono = TRUE,
           stop = FALSE)

get_bounds(gammas = Y_probs$gammas,
           thetas = X_probs$thetas,
           x_mono = TRUE,
           y_mono = TRUE,
           stop = FALSE)

## ----results="asis"-----------------------------------------------------------
counts <- array(data = c(74, 0, 11514, 0,
                         34, 12, 2385, 9663),
                dim = c(2,2,2),
                dimnames = list(x = 0:1,
                                y = 0:1, 
                                z = 0:1))

pander::pandoc.table(counts, style = "rmarkdown")

## ----results="asis"-----------------------------------------------------------
zetas <- counts / rep(apply(counts, 3, sum), each = 4)

pander::pandoc.table(zetas, style = "rmarkdown")

## -----------------------------------------------------------------------------
get_bounds(zetas = zetas)

## -----------------------------------------------------------------------------
thetas <- c(0, 0.8)
gammas <- c(0.994, 0.996)


get_bounds(gammas = gammas,
           thetas = thetas)

