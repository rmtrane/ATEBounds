# Non-parametric Bounds in Two-Sample Summary-Data Mendelian Randomization: Some Cautionary Tales for Practice

This website provides full documentation of the code used for our paper,
and step-by-step guides on how to reproduce every figure and result. For
the curious, it was generated using
[`pkgdown`](https://pkgdown.r-lib.org).

The R-package `ATEBounds` includes every bit of code needed to reproduce
all figures and results included in our paper. Since the main purpose is
to ensure the reproducibility of figures and results from this specific
paper, it has not been thoroughly tested in other settings. The code has
not been optimized for performance, and was mainly written and
documented for internal use only. That being said, many parts of this
can be directly applied to other data sets, or reused in different
settings.

You can install this package using
`remotes::install_github("rmtrane/ATEBounds")`. This site includes an
introduction to finding nonparametric bounds on the ATE using this
package ([Get started](articles/ATEBounds.html)) and documentation of
the functions included. On top of that, there are a few pages with
information on how to reproduce our results; see below for more. One
that might be of particular interest is `vignette("example_analysis")`,
which shows how to use this package to obtain two-sample bounds based on
studies from the MR-Base database.

## Introduction to our work

Nonparametric bounds of the Average Treatment Effect (ATE) are
intriguing for several reasons. Getting bounds on the ATE with no
further assumptions other than what is embedded in the IV model is very
appealing in particular as such assumptions often limit the structure of
unmeasured confounders. By definition, assumptions about unmeasured
confounders cannot be checked, and it is therefore desirable to avoid
such.

Many have considered nonparametric bounds in the past. Swanson et al. (2018) provide a nice overview of different bounds derived from
different sets of assumptions. Arguably the weakest set of assumptions
result in what is known as the Balke-Pearl bounds, first presented in
Balke and Pearl (1993).

The same set of bounds can be obtained using a geometric approach as
presented by Ramsahai (2012). We will rely on this work as it provides a
very general approach to finding bounds. This allows us to obtain bounds
in a setting where we do not have observations of (*X*, *Y*, *Z*), but
rather observations of (*X*, *Z*) and (*Y*, *Z*) from two different data
sources. (We will refer to the former as the one-sample (or trivariate)
setting, and the latter as the two-sample (or bivariate) setting.) This
is relevant today in many mendelian randomization (MR) analyses that
have a natural desire to take advantage of the large databases readily
available. These databases contain results from GWAS in the form of
summary statistics. Each GWAS is essentially exploring associations
between a feature of interest (smoking, obesity, incidence of heart
attack, etc.) and many thousands of SNPs. When the feature of interest
is binary, the summary statistics are often obtained using logistic
regression.

The behavior of nonparametric bounds based on one-sample data is pretty
well-studied, but little is known about the behavior of nonparametric
bounds based on two-sample data. The goal of this work is to shine a
light on when nonparametric bounds can be useful particularly in
two-sample MR studies based on GWAS summary statistics. We do this based
on two aspects of the bounds: (1) the length of the bounds and (2)
whether the bounds cover the null effect of zero.

Our paper contains the following sections:

-   **Section 2** reviews relevant notation, definitions and
    assumptions, and introduces the two-sample IV bounds derived by
    Ramsahai (2012).
-   **Section 3** shows that even when adding two strong monotonicity
    assumptions, the width of two-sample bounds are bounded above by
    2 − 2 ⋅ ST (where
    ST = max<sub>*z*<sub>1</sub> ≠ *z*<sub>2</sub></sub> = |*P*(*X* = 1|*Z* = *z*<sub>1</sub>) − *P*(*X* = 2|*Z* = *z*<sub>2</sub>)|
    is the measure of strength of the instrument we consider). This is
    important as one-sample bounds are guaranteed to have width at
    most 1. We show through simulations that this upper bounds is sharp
    also when the two monotonicity assumptions are not imposed. Figure
    1a of our paper shows this, and is created in
    `vignette("bounds_from_bivariate")`. Since our ST metric isn’t
    common in MR studies, we illustrate the relationship between the ST
    measure and coefficients from a logistic model, and show what kind
    of coefficient is needed for two-sample bounds to exclude the null
    effect 0 for various effect sizes. The latter can to some extend be
    thought of as a “power analysis,” except we work under the
    assumption that we have population based probabilities. Figures 1b
    and 2 are created in `vignette("power_sims_figures")`. Finally, we
    mention the use of multiple instruments. The main part of this
    discussion is presented in eAppendix A.5. Simulations for this
    eAppendix are in `vignette("multiple_IVs_sims")` and figures are
    created in `vignette("multiple_IVs_figures")`.
-   **Section 4** illustrates the information loss we see when using a
    two-sample rather than a one-sample design.
    `vignette("bounds_from_trivariate")` shows how Figure 3 was created.
-   **Section 5** demonstrates our findings on two real MR studies.
    `vignette("example_analysis")` shows how we obtained the data and
    preprocessed it using the `TwoSampleMR` R-package (Hemani et
    al. 2018), and how Figures 4 and 5 were created. It also includes
    the code needed to reproduce eFigures and eTables in eAppendix A.7.

## References

Balke, Alexander, and Judea Pearl. 1993. “Nonparametric Bounds on Causal
Effects from Partial Compliance Data.” JOURNAL OF THE AMERICAN
STATISTICAL ASSOCIATION.

Hemani, Gibran, Jie Zheng, Benjamin Elsworth, Kaitlin H Wade, Valeriia
Haberland, Denis Baird, Charles Laurin, et al. 2018. “The MR-Base
Platform Supports Systematic Causal Inference Across the Human Phenome.”
Edited by Ruth Loos. *eLife* 7 (May): e34408.
<https://doi.org/10.7554/eLife.34408>.

Ramsahai, Roland R. 2012. “Causal Bounds and Observable Constraints for
Non-Deterministic Models.” *J. Mach. Learn. Res.* 13 (March): 829–48.

Swanson, Sonja A., Miguel A. Hern’an, Matthew Miller, James M. Robins,
and Thomas S. Richardson. 2018. “Partial Identification of the Average
Treatment Effect Using Instrumental Variables: Review of Methods for
Binary Instruments, Treatments, and Outcomes.” *Journal of the American
Statistical Association* 113 (522): 933–47.
<https://doi.org/10.1080/01621459.2018.1434530>.
