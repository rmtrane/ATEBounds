# Tiny Bit of Background

Bounds constructed using two-sample data generally provides much less information than bounds constructed using one-sample data. This small Shiny application lets you input values for $P(X = 1 | Z = z)$ and $P(Y = 1 | Z = z)$ for which potential one-sample bounds are sampled. This allows us to get a sense of whether an experiment where two-sample bounds are non-informative could have resulted in informative one-sample bound had we used a one-sample study design. 

For more details, see [this vignette](https://rmtrane.github.io/ATEBounds/articles/bounds_from_trivariate.html).

# How to use this application

This application is very simple and hopefully easy to use. On the panel to the left, you have the option to assume monotonicity of $P(X = 1 | Z = z)$, and specify the number of one-sample bounds you want to sample. You can then either put in your own values for $P(X = 1 | Z = z)$ and $P(Y = 1 | Z = z)$, or you can ask for randomly generated values. (The latter is meant for illustrative purposes.) If you choose to generate random values and have selected to assume monotonicity, the randomly generated values will satisfy the monotonicity assumption. Once you hit the "GO!" button, one-sample bounds are randomly drawn and a figure displaying these bounds and the two-sample bounds is generated. The results can be found under the "Potential One-Sample Bounds" tab. 

If you see a figure that looks promising and you'd like to save it for later, you can click "Save plot". This defaults to download a 300 dpi 8in x 4in .png version of the plot. The data used to create that plot is also available through the "Save one-sample data". This saves a `tibble` to a .Rds file. The `tibble` contains the following columns:

* `id`: simple numeric column to keep track of the samples. This is used on the y-axis of the figure.
* `covs`: a list column with the results of using the [`ATEBounds::potential_covs`](https://rmtrane.github.io/ATEBounds/reference/potential_covs.html) function. This function returns the constraints we have on the values of $\text{Cov}(X,Y | Z = z)$. See the vignette linked above for a few more details.
* `joint`: a list column with 2x2x3 matrices giving the sampled probabilities of $P(X = x, Y = y | Z = z)$
* `lower` and `upper`: numeric columns giving the one-sample bounds
* `n_rejected`: number of times a set of covariances were proposed and rejected. 
* `contains_zero`: indicates if the one-sample bounds overlap zero or not. Used to color bounds in figure.
* `two_sample_lower` and `two_sample_upper`: numeric columns giving the two-sample bounds based on the specified values of $P(X = 1 | Z = z)$ and $P(Y = 1 | Z = 1)$. 

# Advanced Options

At the bottom of the "Potential One-Sample Bounds", you find more advanced settings. Here, you can see (or specify) the seed used for generating the one-sample bounds (for reproducibility). You can also modify the height of the plot you see (especially useful if you want to use a larger number of iterations), and change the settings used to save the plot.
