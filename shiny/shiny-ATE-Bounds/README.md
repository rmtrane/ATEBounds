# Welcome!

This shiny application can be used to explore the behavior of two-sample bounds. It consists of two parts:

1. Under "Two-Sample Bounds from Simulated Data" you can specify the parameters in a logistic model for binary exposure and binary outcome. Based on these, data are simulated, and nonparametric bounds on the ATE are found as if two-sample data were observed. I.e. bounds are found using $P(X = 1 | Z = z)$ and $P(Y = 1 | Z = z)$. The model is flexible enough that you can include up to 10 independent (valid or invalid) instruments and a pair of dependent instruments. 

2. Under "One-Sample Bounds from Two-Sample Summary Statistics" you can explore the loss of information associated with using a two-sample instead of a one-sample design. Given a set of two-sample probabilities, i.e. values of $P(X = 1 | Z = z)$ and $P(Y = 1 | Z = z)$, potential one-sample bounds are sampled and displayed within the two-sample bounds. 

The heavy lifting is done using the R package [`ATEBounds`](https://rmtrane.github.io/ATEBounds). If you wish to play around with your own data, you can install it using `devtools::install_github("rmtrane/ATEBounds")`. If you want to run this application locally, you can do so using `shiny::runGitHub("rmtrane/ATEBounds", ref = "ATEBounds", subdir = "shiny/shiny-ATE-bounds")`. (Make sure you have installed the `ATEBounds` package in addition to the following before you do so: `DT`, `tidyverse`, `shiny`, `dagitty`, `ggdag`, `ggrepel`, `shinyjs`, `shinysky`.) If you encounter any bugs or have any questions/comments, please open a GitHub issue [here](https://github.com/rmtrane/ATEBounds/issues). 
