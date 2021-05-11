# regmedint <img src="man/figures/hex.png" align="right" height="140"/>

[![Travis-CI Build
Status](https://travis-ci.org/kaz-yos/regmedint.svg?branch=master)](https://travis-ci.org/kaz-yos/regmedint)
[![](http://www.r-pkg.org/badges/version/regmedint)](http://www.r-pkg.org/pkg/regmedint)
[![CRAN RStudio mirror
downloads](http://cranlogs.r-pkg.org/badges/regmedint)](http://www.r-pkg.org/pkg/regmedint)

This is an R re-implementation of the regression-based causal mediation
analysis method, supporting a treatment-mediator interaction term, as
implemented in the SAS macro by Valeri and VanderWeele (2013 and 2015).
The original is found at Dr. VanderWeele’s [Tools and
Tutorials](https://www.hsph.harvard.edu/tyler-vanderweele/tools-and-tutorials/).

This package is meant to be an educational tool. Thanks to R’s
expressibility, the code is likely easier to read than SAS IML. Thus,
the correspondence between the mathematical formulas presented in the
Appendix of Explanation in Causal Inference (V2015) and the code should
be easier to grasp. See the vignette (Article on the package website)
for the code that implements formulas.

To cite this software, please refer to <https://osf.io/6c79f/>.

# Implemented models

This package is currently under development and is not on CRAN, yet.
Following VV2015, the following grid of models are implemented. yreg
refers to the outcome model and mreg refers to the mediator model. See
the [Github
repo](https://github.com/kaz-yos/regmedint/tree/master/tests/reference_results)
for cross-check results against the SAS macro.

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">

<table>
<thead>
<tr class="header">
<th style="text-align: left;">yreg \\ mreg</th>
<th style="text-align: center;">linear</th>
<th style="text-align: center;">logistic</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">linear</td>
<td style="text-align: center;">implemented</td>
<td style="text-align: center;">implemented</td>
</tr>
<tr class="even">
<td style="text-align: left;">logistic (1)</td>
<td style="text-align: center;">implemented</td>
<td style="text-align: center;">implemented</td>
</tr>
<tr class="odd">
<td style="text-align: left;">loglinear</td>
<td style="text-align: center;">implemented (2)</td>
<td style="text-align: center;">implemented (2)</td>
</tr>
<tr class="even">
<td style="text-align: left;">poisson</td>
<td style="text-align: center;">implemented</td>
<td style="text-align: center;">implemented</td>
</tr>
<tr class="odd">
<td style="text-align: left;">negbin</td>
<td style="text-align: center;">implemented</td>
<td style="text-align: center;">implemented</td>
</tr>
<tr class="even">
<td style="text-align: left;">survCox (1)</td>
<td style="text-align: center;">implemented</td>
<td style="text-align: center;">implemented</td>
</tr>
<tr class="odd">
<td style="text-align: left;">survAFT exp</td>
<td style="text-align: center;">implemented</td>
<td style="text-align: center;">implemented</td>
</tr>
<tr class="even">
<td style="text-align: left;">survAFT weibull</td>
<td style="text-align: center;">implemented</td>
<td style="text-align: center;">implemented</td>
</tr>
</tbody>
</table>

1.  Approximation depends on the rare event assumptions.

2.  As oppose to the SAS macro, the outcome model is implemented as a
    modified Poisson model (log link with robust variance) as in Z2004.

See the corresponding vignettes (Articles on the package website) for
how to perform bootstrapping and multiple imputation.

# Installation

For the developmental version on Github, use the following commands to
install the package.

    install.packages("devtools") # If you do not have devtools already.
    devtools::install_github("kaz-yos/regmedint")

Once released on CRAN, the following should install the stable version.

    install.packages("regmedint")

# Example

Here we will analyze the simulated example that is included with the SAS
macro. See the introduction vignette (Article on the package website)
for details of the user interface functions.

## `regmedint()` to fit models

The `regmedint` function is the user interface for constructing a result
object of class `regmedint`. The interface is made similar to the
original SAS macro. For survival outcomes, the indicator variable is an
event indicator (1 for event, 0 for censoring). The `c_cond` vector is
required. This vector is the vector of covariate values at which the
conditional effects are evaluated at.

    library(regmedint)

    ## Error in library(regmedint): there is no package called 'regmedint'

    library(tidyverse)

    ## Example data
    data(vv2015)

    regmedint_obj <- regmedint(data = vv2015,
                               ## Variables
                               yvar = "y",
                               avar = "x",
                               mvar = "m",
                               cvar = c("c"),
                               eventvar = "event",
                               ## Values at which effects are evaluated
                               a0 = 0,
                               a1 = 1,
                               m_cde = 1,
                               c_cond = 0.5,
                               ## Model types
                               mreg = "logistic",
                               yreg = "survAFT_weibull",
                               ## Additional specification
                               interaction = TRUE,
                               casecontrol = FALSE)

    ## Error in regmedint(data = vv2015, yvar = "y", avar = "x", mvar = "m", : could not find function "regmedint"

## `print()` to examine simplified results

Implicit printing prints mreg, yreg, and mediation analysis point
estimates. All effect estimates are on the scale of the link function.

    regmedint_obj

    ## Error in eval(expr, envir, enclos): object 'regmedint_obj' not found

## `summary()` to examine extended results

The `summary` method gives the summary for `mreg`, `yreg`, and mediation
analysis results. The `exponentiate` option will add the exponentiated
estimate and confidence limits if the outcome model is not a linear
model. The pure natural direct effect (`pnde`) is what is typically
called the natural direct effect (NDE). The total natural indirect
effect (`tnie`) is the corresponding natural indirect effect (NIE).

    summary(regmedint_obj, exponentiate = TRUE)

    ## Error in summary(regmedint_obj, exponentiate = TRUE): object 'regmedint_obj' not found

Use `coef` to extract the mediation analysis results only.

    coef(summary(regmedint_obj, exponentiate = TRUE))

    ## Error in summary(regmedint_obj, exponentiate = TRUE): object 'regmedint_obj' not found

Note that the estimates can be re-evaluated at different `m_cde` and
`c_cond` without re-fitting the underlyng models.

    coef(summary(regmedint_obj, exponentiate = TRUE, m_cde = 0, c_cond = 5))

    ## Error in summary(regmedint_obj, exponentiate = TRUE, m_cde = 0, c_cond = 5): object 'regmedint_obj' not found

# Formulas

The point estimate and standard error formulas (multivariate delta
method) were taken from the following references. See
<https://github.com/kaz-yos/regmedint-supplement/blob/master/supplement.pdf>
for type-set formulas.

-   V2015: VanderWeele (2015) Explanation in Causal Inference.
-   VV2013A: Valeri & VanderWeele (2013) Appendix
-   VV2015A: Valeri & VanderWeele (2015) Appendix

As seen below, there are only four patterns.

## Effect formulas

<table>
<colgroup>
<col style="width: 18%" />
<col style="width: 40%" />
<col style="width: 40%" />
</colgroup>
<thead>
<tr class="header">
<th>yreg \\ mreg</th>
<th>linear</th>
<th>logistic</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>linear</td>
<td>V2015 p466 Proposition 2.3</td>
<td>V2015 p471 Proposition 2.5</td>
</tr>
<tr class="even">
<td></td>
<td></td>
<td></td>
</tr>
<tr class="odd">
<td>logistic</td>
<td>V2015 p468 Proposition 2.4</td>
<td>V2015 p473 Proposition 2.6</td>
</tr>
<tr class="even">
<td>loglinear</td>
<td>VV2013A p8 Use Proposition 2.4</td>
<td>VV2013A p8 Use Proposition 2.6</td>
</tr>
<tr class="odd">
<td>poisson</td>
<td>VV2013A p8 Use Proposition 2.4</td>
<td>VV2013A p8 Use Proposition 2.6</td>
</tr>
<tr class="even">
<td>negbin</td>
<td>VV2013A p8 Use Proposition 2.4</td>
<td>VV2013A p8 Use Proposition 2.6</td>
</tr>
<tr class="odd">
<td></td>
<td></td>
<td></td>
</tr>
<tr class="even">
<td>survCox</td>
<td>V2015 p496 Proposition 4.4 (Use 2.4)</td>
<td>V2015 p499 Proposition 4.6 (Use 2.6)</td>
</tr>
<tr class="odd">
<td>survAFT exp</td>
<td>V2015 p494 Proposition 4.1 (Use 2.4)</td>
<td>V2015 p495 Proposition 4.3 (Use 2.6)</td>
</tr>
<tr class="even">
<td>survAFT weibull</td>
<td>V2015 p494 Proposition 4.1 (Use 2.4)</td>
<td>V2015 p495 Proposition 4.3 (Use 2.6)</td>
</tr>
</tbody>
</table>

Note the loglinear outcome model means log link with binomial error.

## Standard error formulas

<table>
<colgroup>
<col style="width: 20%" />
<col style="width: 39%" />
<col style="width: 39%" />
</colgroup>
<thead>
<tr class="header">
<th>yreg \\ mreg</th>
<th>linear</th>
<th>logistic</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>linear</td>
<td>V2015 p466 Proposition 2.3</td>
<td>V2015 p471 Proposition 2.5</td>
</tr>
<tr class="even">
<td></td>
<td></td>
<td></td>
</tr>
<tr class="odd">
<td>logistic</td>
<td>V2015 p468 Proposition 2.4</td>
<td>V2015 p473 Proposition 2.6</td>
</tr>
<tr class="even">
<td>loglinear</td>
<td>VV2013A p8 Use Proposition 2.4</td>
<td>VV2013A p8 Use Proposition 2.6</td>
</tr>
<tr class="odd">
<td>poisson</td>
<td>VV2013A p8 Use Proposition 2.4</td>
<td>VV2013A p8 Use Proposition 2.6</td>
</tr>
<tr class="even">
<td>negbin</td>
<td>VV2013A p8 Use Proposition 2.4</td>
<td>VV2013A p8 Use Proposition 2.6</td>
</tr>
<tr class="odd">
<td></td>
<td></td>
<td></td>
</tr>
<tr class="even">
<td>survCox</td>
<td>V2015 p496 Use Proposition 2.4</td>
<td>V2015 p499 Use Proposition 2.6</td>
</tr>
<tr class="odd">
<td>survAFT exp</td>
<td>V2015 p494 Use Proposition 2.4</td>
<td>V2015 p495 Use Proposition 2.6</td>
</tr>
<tr class="even">
<td>survAFT weibull</td>
<td>V2015 p494 Use Proposition 2.4</td>
<td>V2015 p495 Use Proposition 2.6</td>
</tr>
</tbody>
</table>

# Design

The software design is outlined here for those who may be interested.

-   Call structure
    -   regmedint UI function
        -   new\_regmedint internal constructor
            -   fit\_mreg
            -   fit\_yreg
            -   calc\_myreg calls a specialized worker function, which
                return two functions, one for point estimates and the
                other for standard error estimate.
                -   calc\_myreg\_mreg\_linear\_yreg\_linear
                -   calc\_myreg\_mreg\_linear\_yreg\_logistic
                -   calc\_myreg\_mreg\_logistic\_yreg\_linear
                -   calc\_myreg\_mreg\_logistic\_yreg\_logistic
-   regmedint object structure
    -   mreg\_fit mediator regression model object as is
    -   yreg\_fit outcome regression model object as is
    -   myreg\_funs list
        -   est\_fun: (a0,a1,m\_cde,c\_cond) →
            (cde,pnde,tnie,tnde,pnie,te,pm)
        -   se\_fun: (a0,a1,m\_cde,c\_cond) → se for
            (cde,pnde,tnie,tnde,pnie,te,pm)
        -   args preserves arguments given to the UI
-   User methods for the regmedint object
    -   print.regmedint: prints coefficients for mreg, yreg, and
        mediation analysis
    -   summary.regmedint: regmedint → summary\_regmedint
        -   print.summary\_regmedint: prints summary objects for mreg,
            yreg, and mediation analysis
        -   coef.summary\_regmedint:
    -   coef.regmedint: regmedint → vector
        (cde,pnde,tnie,tnde,pnie,te,pm)
    -   vcov.regmedint: regmedint → matrix
        (cde,pnde,tnie,tnde,pnie,te,pm). Off-diagonals are NA.
    -   confint.regmedint: regmedint → matrix of (lower,upper)

# Similar or related projects for counterfactual causal mediation analysis

## R

-   mediation (simulation-based):
    <https://CRAN.R-project.org/package=mediation>
-   medflex (natural effect model):
    <https://CRAN.R-project.org/package=medflex>
-   intmed (interventional analogue):
    <https://CRAN.R-project.org/package=intmed>
-   mediator (regression-based): <https:/github.com/GerkeLab/mediator>
-   causalMediation (regression-based):
    <https:/github.com/harvard-P01/causalMediation>

## Other statistical environment

-   SAS macro (original regression-based)
    <https://www.hsph.harvard.edu/tyler-vanderweele/tools-and-tutorials/>
-   SAS PROC CAUSALMED (regression-based)
    <https://support.sas.com/rnd/app/stat/procedures/causalmed.html>

# References

-   V2015: VanderWeele (2015) Explanation in Causal Inference.
-   VV2013: Valeri & VanderWeele (2013) Psych Method. 18:137.
-   VV2015: Valeri & VanderWeele (2015) Epidemiology. 26:e23.
-   Z2004: Zou (2004) Am J Epidemiol 159:702.
