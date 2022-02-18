regmedint <img src="man/figures/hex.png" align="right" height="140"/>
=====================================================================

[![R-CMD-check](https://github.com/kaz-yos/regmedint/workflows/R-CMD-check/badge.svg)](https://github.com/kaz-yos/regmedint/actions)
[![](http://www.r-pkg.org/badges/version/regmedint)](http://www.r-pkg.org/pkg/regmedint)
[![CRAN RStudio mirror
downloads](http://cranlogs.r-pkg.org/badges/regmedint)](http://www.r-pkg.org/pkg/regmedint)

This is an extension of the regression-based causal mediation analysis
first proposed by Valeri and VanderWeele (2013) and Valeri and
VanderWeele (2015). The current version supports including effect
measure modification by covariates (treatment-covariate and
mediator-covariate product terms in mediator and outcome regression
models). It also accommodates the original ‘SAS’ macro (can be found at
Dr. VanderWeele’s [Tools and
Tutorials](https://www.hsph.harvard.edu/tyler-vanderweele/tools-and-tutorials/))
and PROC CAUSALMED procedure in ‘SAS’ when there is no effect measure
modification. Linear and logistic models are supported for the mediator
model. Linear, logistic, loglinear, Poisson, negative binomial, Cox, and
accelerated failure time (exponential and Weibull) models are supported
for the outcome model.

To cite this software, please use: *regmedint* (v1.0.0; Yoshida, Li, &
Mathur, 2021)
<!-- cite the r package:  https://easystats.github.io/report/articles/cite_packages.html -->

Implemented models
==================

The following grid of models are implemented. `yreg` refers to the
outcome model and `mreg` refers to the mediator model.

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
<td style="text-align: center;">:heavy_check_mark:</td>
<td style="text-align: center;">:heavy_check_mark:</td>
</tr>
<tr class="even">
<td style="text-align: left;">logistic<span class="math inline"><sup>1</sup></span></td>
<td style="text-align: center;">:heavy_check_mark:</td>
<td style="text-align: center;">:heavy_check_mark:</td>
</tr>
<tr class="odd">
<td style="text-align: left;">loglinear</td>
<td style="text-align: center;">:heavy_check_mark:<span class="math inline"><sup>2</sup></span></td>
<td style="text-align: center;">:heavy_check_mark:<span class="math inline"><sup>2</sup></span></td>
</tr>
<tr class="even">
<td style="text-align: left;">poisson</td>
<td style="text-align: center;">:heavy_check_mark:</td>
<td style="text-align: center;">:heavy_check_mark:</td>
</tr>
<tr class="odd">
<td style="text-align: left;">negbin</td>
<td style="text-align: center;">:heavy_check_mark:</td>
<td style="text-align: center;">:heavy_check_mark:</td>
</tr>
<tr class="even">
<td style="text-align: left;">survCox<span class="math inline"><sup>1</sup></span></td>
<td style="text-align: center;">:heavy_check_mark:</td>
<td style="text-align: center;">:heavy_check_mark:</td>
</tr>
<tr class="odd">
<td style="text-align: left;">survAFT exp</td>
<td style="text-align: center;">:heavy_check_mark:</td>
<td style="text-align: center;">:heavy_check_mark:</td>
</tr>
<tr class="even">
<td style="text-align: left;">survAFT weibull</td>
<td style="text-align: center;">:heavy_check_mark:</td>
<td style="text-align: center;">:heavy_check_mark:</td>
</tr>
</tbody>
</table>

<sup>1</sup> Approximation depends on the rare event assumptions.

<sup>2</sup> Implemented as a modified Poisson model (log link with
robust variance) as in Z2004.

See the corresponding vignettes (Articles on the package website) for
how to perform bootstrapping and multiple imputation.

Installation
============

For the developmental version on Github, use the following commands to
install the package.

    # install.packages("devtools") # If you do not have devtools already.
    devtools::install_github("kaz-yos/regmedint")

    install.packages("regmedint")

Data Example
============

We use `VV2015` dataset for demonstration.

    library(regmedint)
    data(vv2015)

`regmedint()` to fit models
---------------------------

The `regmedint` function is the user interface for constructing a result
object of class `regmedint`. The interface is similar to the original
SAS macro. For survival outcomes, the indicator variable is an event
indicator (1 for event, 0 for censoring). `c_cond` vector is required be
specified. This vector is the vector of covariate values at which the
conditional effects are evaluated at.

1.  When there is no effect measure modification by covariates,
    `emm_ac_mreg = NULL`, `emm_ac_yreg = NULL`, `emm_mc_yreg = NULL`.

<!-- -->

    regmedint_obj1 <- regmedint(data = vv2015,
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
                                c_cond = 3,
                                ## Model types
                                mreg = "logistic",
                                yreg = "survAFT_weibull",
                                ## Additional specification
                                interaction = TRUE,
                                casecontrol = FALSE)
     summary(regmedint_obj1)

    ## ### Mediator model
    ## 
    ## Call:
    ## glm(formula = m ~ x + c, family = binomial(link = "logit"), data = data)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.5143  -1.1765   0.9177   1.1133   1.4602  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)
    ## (Intercept)  -0.3545     0.3252  -1.090    0.276
    ## x             0.3842     0.4165   0.922    0.356
    ## c             0.2694     0.2058   1.309    0.191
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 138.59  on 99  degrees of freedom
    ## Residual deviance: 136.08  on 97  degrees of freedom
    ## AIC: 142.08
    ## 
    ## Number of Fisher Scoring iterations: 4
    ## 
    ## ### Outcome model
    ## 
    ## Call:
    ## survival::survreg(formula = Surv(y, event) ~ x + m + x:m + c, 
    ##     data = data, dist = "weibull")
    ##               Value Std. Error     z           p
    ## (Intercept) -1.0424     0.1903 -5.48 0.000000043
    ## x            0.4408     0.3008  1.47        0.14
    ## m            0.0905     0.2683  0.34        0.74
    ## c           -0.0669     0.0915 -0.73        0.46
    ## x:m          0.1003     0.4207  0.24        0.81
    ## Log(scale)  -0.0347     0.0810 -0.43        0.67
    ## 
    ## Scale= 0.966 
    ## 
    ## Weibull distribution
    ## Loglik(model)= -11.4   Loglik(intercept only)= -14.5
    ##  Chisq= 6.31 on 4 degrees of freedom, p= 0.18 
    ## Number of Newton-Raphson Iterations: 5 
    ## n= 100 
    ## 
    ## ### Mediation analysis 
    ##              est         se         Z          p       lower      upper
    ## cde  0.541070807 0.29422958 1.8389409 0.06592388 -0.03560858 1.11775019
    ## pnde 0.505391952 0.21797147 2.3186151 0.02041591  0.07817572 0.93260819
    ## tnie 0.015988820 0.03171597 0.5041252 0.61417338 -0.04617334 0.07815098
    ## tnde 0.513662425 0.22946248 2.2385465 0.02518544  0.06392423 0.96340062
    ## pnie 0.007718348 0.02398457 0.3218047 0.74760066 -0.03929055 0.05472725
    ## te   0.521380773 0.22427066 2.3247837 0.02008353  0.08181835 0.96094319
    ## pm   0.039039346 0.07444080 0.5244348 0.59997616 -0.10686194 0.18494063
    ## 
    ## Evaluated at:
    ## avar: x
    ##  a1 (intervened value of avar) = 1
    ##  a0 (reference value of avar)  = 0
    ## mvar: m
    ##  m_cde (intervend value of mvar for cde) = 1
    ## cvar: c
    ##  c_cond (covariate vector value) = 3
    ## 
    ## Note that effect estimates can vary over m_cde and c_cond values when interaction = TRUE.

1.  When there is effect measure modification by covariates,
    `emm_ac_mreg`, `emm_ac_yreg` and `emm_mc_yreg` can take a sub-vector
    of covariates in `cvar`.

<!-- -->

    regmedint_obj2 <- regmedint(data = vv2015,
                                ## Variables
                                yvar = "y",
                                avar = "x",
                                mvar = "m",
                                cvar = c("c"),
                                emm_ac_mreg = c("c"),
                                emm_ac_yreg = c("c"),
                                emm_mc_yreg = c("c"),
                                eventvar = "event",
                                ## Values at which effects are evaluated
                                a0 = 0,
                                a1 = 1,
                                m_cde = 1,
                                c_cond = 3,
                                ## Model types
                                mreg = "logistic",
                                yreg = "survAFT_weibull",
                                ## Additional specification
                                interaction = TRUE,
                                casecontrol = FALSE)

    ## Error in regmedint(data = vv2015, yvar = "y", avar = "x", mvar = "m", : unused arguments (emm_ac_mreg = c("c"), emm_ac_yreg = c("c"), emm_mc_yreg = c("c"))

     summary(regmedint_obj2)

    ## Error in summary(regmedint_obj2): object 'regmedint_obj2' not found

<!-- ## `print()` to examine simplified results -->
<!-- Implicit printing prints `mreg`, `yreg`, and mediation analysis point estimates. All effect estimates are on the scale of the link function. -->

`summary()` to examine extended results
---------------------------------------

The `summary` method gives the summary for `mreg`, `yreg`, and mediation
analysis results. The `exponentiate` option will add the exponentiated
estimate and confidence limits if the outcome model is not a linear
model. The pure natural direct effect (`pnde`) is what is typically
called the natural direct effect (NDE). The total natural indirect
effect (`tnie`) is the corresponding natural indirect effect (NIE).

    summary(regmedint_obj2, exponentiate = TRUE)

    ## Error in summary(regmedint_obj2, exponentiate = TRUE): object 'regmedint_obj2' not found

Use `coef` to extract the mediation analysis results only.

    coef(summary(regmedint_obj2, exponentiate = TRUE))

    ## Error in summary(regmedint_obj2, exponentiate = TRUE): object 'regmedint_obj2' not found

Note that the estimates can be re-evaluated at different `m_cde` and
`c_cond` without re-fitting the underlyng models.

    coef(summary(regmedint_obj2, exponentiate = TRUE,
                 m_cde = min(heart.disease$bp),
                 c_cond = c(min(heart.disease$age), 0, 1,
                                           max(heart.disease$bloodsugar_T),
                                           min(heart.disease$maximumHR))))

    ## Error in summary(regmedint_obj2, exponentiate = TRUE, m_cde = min(heart.disease$bp), : object 'regmedint_obj2' not found

Formulas
========

See [here](https://osf.io/d4brv/) for the following formulas.

Effect formulas in the supplementary document
---------------------------------------------

<table>
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
<td>Formulas (1) - (5)</td>
<td>Formulas (11) - (15)</td>
</tr>
<tr class="even">
<td></td>
<td></td>
<td></td>
</tr>
<tr class="odd">
<td>logistic</td>
<td>Formulas (21) - (25)</td>
<td>Formulas (31) - (35)</td>
</tr>
<tr class="even">
<td>loglinear</td>
<td>Formulas (21) - (25)</td>
<td>Formulas (31) - (35)</td>
</tr>
<tr class="odd">
<td>poisson</td>
<td>Formulas (21) - (25)</td>
<td>Formulas (31) - (35)</td>
</tr>
<tr class="even">
<td>negbin</td>
<td>Formulas (21) - (25)</td>
<td>Formulas (31) - (35)</td>
</tr>
<tr class="odd">
<td></td>
<td></td>
<td></td>
</tr>
<tr class="even">
<td>survCox</td>
<td>Formulas (21) - (25)</td>
<td>Formulas (31) - (35)</td>
</tr>
<tr class="odd">
<td>survAFT exp</td>
<td>Formulas (21) - (25)</td>
<td>Formulas (31) - (35)</td>
</tr>
<tr class="even">
<td>survAFT weibull</td>
<td>Formulas (21) - (25)</td>
<td>Formulas (31) - (35)</td>
</tr>
</tbody>
</table>

Standard error formulas in the supplementary document
-----------------------------------------------------

<table>
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
<td>Formulas (6) - (10)</td>
<td>Formulas (16) - (20)</td>
</tr>
<tr class="even">
<td></td>
<td></td>
<td></td>
</tr>
<tr class="odd">
<td>logistic</td>
<td>Formulas (26) - (30)</td>
<td>Formulas (36) - (40)</td>
</tr>
<tr class="even">
<td>loglinear</td>
<td>Formulas (26) - (30)</td>
<td>Formulas (36) - (40)</td>
</tr>
<tr class="odd">
<td>poisson</td>
<td>Formulas (26) - (30)</td>
<td>Formulas (36) - (40)</td>
</tr>
<tr class="even">
<td>negbin</td>
<td>Formulas (26) - (30)</td>
<td>Formulas (36) - (40)</td>
</tr>
<tr class="odd">
<td></td>
<td></td>
<td></td>
</tr>
<tr class="even">
<td>survCox</td>
<td>Formulas (26) - (30)</td>
<td>Formulas (36) - (40)</td>
</tr>
<tr class="odd">
<td>survAFT exp</td>
<td>Formulas (26) - (30)</td>
<td>Formulas (36) - (40)</td>
</tr>
<tr class="even">
<td>survAFT weibull</td>
<td>Formulas (26) - (30)</td>
<td>Formulas (36) - (40)</td>
</tr>
</tbody>
</table>

Note: The point estimate and standard error formulas (multivariate delta
method) were derived based on the following references.

-   V2015: VanderWeele (2015) Explanation in Causal Inference.
-   VV2013A: Valeri & VanderWeele (2013) Appendix
-   VV2015A: Valeri & VanderWeele (2015) Appendix

Effect formulas are based on the following propositions
-------------------------------------------------------

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

Standard error formulas are based on the following propositions
---------------------------------------------------------------

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

<!-- # Design -->
<!-- The software design is outlined here for those who may be interested. -->
<!-- - Call structure -->
<!--   - regmedint UI function -->
<!--     - new_regmedint internal constructor -->
<!--       - fit_mreg -->
<!--       - fit_yreg -->
<!--       - calc_myreg calls a specialized worker function, which return two functions, one for point estimates and the other for standard error estimate. -->
<!--         - calc_myreg_mreg_linear_yreg_linear -->
<!--         - calc_myreg_mreg_linear_yreg_logistic -->
<!--         - calc_myreg_mreg_logistic_yreg_linear -->
<!--         - calc_myreg_mreg_logistic_yreg_logistic -->
<!-- - regmedint object structure -->
<!--   - mreg_fit mediator regression model object as is -->
<!--   - yreg_fit outcome regression model object as is -->
<!--   - myreg_funs list -->
<!--     - est_fun: (a0,a1,m_cde,c_cond) → (cde,pnde,tnie,tnde,pnie,te,pm) -->
<!--     - se_fun: (a0,a1,m_cde,c_cond) → se for (cde,pnde,tnie,tnde,pnie,te,pm) -->
<!--     - args preserves arguments given to the UI -->
<!-- - User methods for the regmedint object -->
<!--   - print.regmedint: prints coefficients for mreg, yreg, and mediation analysis -->
<!--   - summary.regmedint: regmedint → summary_regmedint -->
<!--     - print.summary_regmedint: prints summary objects for mreg, yreg, and mediation analysis -->
<!--     - coef.summary_regmedint: -->
<!--   - coef.regmedint: regmedint → vector (cde,pnde,tnie,tnde,pnie,te,pm) -->
<!--   - vcov.regmedint: regmedint → matrix (cde,pnde,tnie,tnde,pnie,te,pm). Off-diagonals are NA. -->
<!--   - confint.regmedint: regmedint → matrix of (lower,upper) -->

Similar or related projects for counterfactual-based causal mediation analysis
==============================================================================

R
-

-   mediation (simulation-based):
    <https://CRAN.R-project.org/package=mediation>
-   medflex (natural effect model):
    <https://CRAN.R-project.org/package=medflex>
-   intmed (interventional analogue):
    <https://CRAN.R-project.org/package=intmed>
-   CMAverse (regression-based approach, weighting-based approach,
    inverse odd-ratio weighting, natural effect model, marginal
    structural model, g-formula approach):
    <https://bs1125.github.io/CMAverse/>
-   mediator (regression-based): <https://github.com/GerkeLab/mediator>
-   causalMediation (regression-based):
    <https://github.com/harvard-P01/causalMediation>

Other statistical environment
-----------------------------

-   SAS macro (original regression-based)
    <https://www.hsph.harvard.edu/tyler-vanderweele/tools-and-tutorials/>
-   SAS PROC CAUSALMED (regression-based)
    <https://support.sas.com/rnd/app/stat/procedures/causalmed.html>

References
==========

-   V2015: VanderWeele (2015) Explanation in Causal Inference.
-   VV2013: Valeri & VanderWeele (2013) Psych Method. 18:137.
-   VV2015: Valeri & VanderWeele (2015) Epidemiology. 26:e23.
-   Z2004: Zou (2004) Am J Epidemiol 159:702.
