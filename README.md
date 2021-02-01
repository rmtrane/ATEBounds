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
`remotes::install_github("rmtrane/ATEBounds")`.

## Introduction

Nonparametric bounds of the Average Treatment Effect (ATE) are
intriguing for several reasons. Getting bounds on the ATE with no
further assumptions other than what is embedded in the IV model is very
appealing in particular as such assumptions often limit the structure of
unmeasured confounders. By definition, assumptions about unmeasured
confounders cannot be checked, and it is therefore desirable to avoid
such.

Many have considered nonparametric bounds in the past. Swanson et al.
(2018) provide a nice overview of different bounds derived from
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
    important as one-sample bounds are guaranteed to have width at most
    1. We show through simulations that this upper bounds is sharp also
    when the two monotonicity assumptions are not imposed. Figure 1a of
    our paper shows this, and is created in
    `vignette("bounds_from_bivariate")`. We then illustrate the
    relationship between the ST measure and coefficients from a logistic
    model, and show what kind of coefficient is needed for two-sample
    bounds to exclude the null effect 0.

## Contents

-   [Get started](articles/ATEBounds.html) is a brief introduction to
    this package, and how to find bounds using it.
-   [Reference](reference/index.html) gives a list with links to the
    documentation of all functions included.
-   Under [Articles](articles/index.html) you will find a set of
    documents that show how results and figures were created. Every
    document provides the seed used in order to make the work 100%
    reproducible. In particular, referencing our paper, Figure 1a is
    created [here](articles/bounds_from_bivaraite.html), Figures 1b, 2
    and eFigure 1 are [here](articles/bounds_from_trivariate.html), all
    figures and results relating to real-world data examples are created
    [here](articles/example_analysis.html), and eFigure 2, 3, 4, 5, and
    6 are created [here](articles/multiple_IVs.html).

<!-- ## Width of Two-Sample Bounds -->
<!-- One nice result for bounds obtained from one-sample data tells us that the width will never exceed $1$. This is important since bounds with width greater than $1$ will always include $0$, the null effect, which means direction of the ATE cannot be determined. Unfortunately, two-sample bounds can be as large as $2$. Particularly, for two-sample MR studies, most bounds are in fact close to or exceeding $1$.  -->
<!-- ## Improving Bounds from Two-Sample Data -->
<!-- Since relatively large databases with two-sample data are already available, it seems worthwhile to explore if we can utilize these data in any way to improve the bounds we obtain. We considered two ways of doing so. -->
<!-- ### Multiple IVs -->
<!-- The simplest way of utilizing multiple IVs is by simply taking intersections of the individual intervals. We will see that in most scenarios, the width of an intersection interval is mainly driven by the strongest IV. In fact, if monotonicity is assumed, it can be shown that the lower and upper bounds are monotonically increasing and decreasing, respectively, as a function of the strength of the IV, all else being equal. This means that the intersection will be exactly the bounds of the strongest IV. We show through simulations that even when all else is not equal, the gain in the width of the intersection is miniscule compared to the best individual bounds. -->
<!-- ### Possible One-Sample Distributions -->
<!-- Since we know that the bounds obtained from one-sample data will always have length less than $1$, and generally be narrower than the two-sample bounds, one could ask what distributions of $(X,Y|Z)$ would be possible given the known distributions of $(X|Z)$, $(Y|Z)$, and assuming the IV model holds. Doing so, one can arrive at what can be thought of as a posterior distribution of the one-sample bounds. This enables us to draw a few different conclusions depending on the specific scenario: -->
<!-- * a one-sample study provides no further insights over the two-sample bounds -->
<!-- * a one-sample study might be very likely to tell us the direction of the effect -->
<!-- * a one-sample study might be very unlikely to tell us the direction of the effect -->
<!-- ## Results/Conclusions -->
<!-- * Bounds based on two-sample data are not very useful -- we need very, very strong IVs to be guaranteed any insights. -->
<!-- * Even when utilizing multiple IVs, two-sample data will in most cases (in our experience) not provide enough information to be useful. -->
<!-- * It seems the best use of two-sample data is to get a sense of what one might be able to get from a one-sample study design. This could be used to decide whether or not pursuing such data for a bound based analysis is worth the time. -->
<!-- ## References -->

Balke, Alexander, and Judea Pearl. 1993. “Nonparametric Bounds on Causal
Effects from Partial Compliance Data.” JOURNAL OF THE AMERICAN
STATISTICAL ASSOCIATION.

Ramsahai, Roland R. 2012. “Causal Bounds and Observable Constraints for
Non-Deterministic Models.” *J. Mach. Learn. Res.* 13 (March): 829–48.

Swanson, Sonja A., Miguel A. Hern’an, Matthew Miller, James M. Robins,
and Thomas S. Richardson. 2018. “Partial Identification of the Average
Treatment Effect Using Instrumental Variables: Review of Methods for
Binary Instruments, Treatments, and Outcomes.” *Journal of the American
Statistical Association* 113 (522): 933–47.
<https://doi.org/10.1080/01621459.2018.1434530>.
