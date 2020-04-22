regmedint <img src="man/figures/hex.png" align="right" height="140"/>
=====================================================================

[![Travis-CI Build
Status](https://travis-ci.org/kaz-yos/regmedint.svg?branch=master)](https://travis-ci.org/kaz-yos/regmedint)
[![](http://www.r-pkg.org/badges/version/regmedint)](http://www.r-pkg.org/pkg/regmedint)
[![CRAN RStudio mirror
downloads](http://cranlogs.r-pkg.org/badges/regmedint)](http://www.r-pkg.org/pkg/regmedint)

This is an R reimplementation of the regression-based causal mediation
analysis methods as implemented in the SAS macro by Valeri and
VanderWeele (2013 and 2015). The original is found at Dr. VanderWeele’s
[Tools and
Tutorials](https://www.hsph.harvard.edu/tyler-vanderweele/tools-and-tutorials/).

This package is meant to be an educational tool. Thanks to R’s
expressibility, the code is likely easier to read than SAS IML. Thus,
the correspondence between the formulas presented in the Appendix of
Explanation in Causal Inference (VanderWeele 2015) and the code should
be easier to grasp.

Example
=======

Here we will analyze the simulated example that is included with the SAS
macro.

`regmedint()` to fit models
---------------------------

The `regmedint` function is the user interface for constructing a result
object of class `regmedint`. The interface is made similar to the
original SAS macro. For survival outcomes, the indicator variable is an
event indicator (1 for event, 0 for censoring). The `c_cond` vector is
required. This vector is the vector of covariate values at which the
conditional effects are evaluated at.

    ## install.packages("devtools") # If you do not have it.
    ## devtools::install_github("kaz-yos/regmedint") # If you have not installed it, yet.
    library(regmedint)
    library(tidyverse)

    vv2015 <- read_delim(file = "./tests/reference_results/data-valeri-vanderweele-2015.txt",
                         delim = " ") %>%
        ## Following R convention, an event indicator is used.
        mutate(event = if_else(cens == 0, 1L, 0L))

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

`print()` to examine simplified results
---------------------------------------

Implicit printing prints mreg, yreg, and mediation analysis point
estimates. All effect estimates are on the scale of the link function.

    regmedint_obj

    ## ### Mediator model
    ## 
    ## Call:  glm(formula = m ~ x + c, family = binomial(link = "logit"), data = data)
    ## 
    ## Coefficients:
    ## (Intercept)            x            c  
    ##     -0.3545       0.3842       0.2694  
    ## 
    ## Degrees of Freedom: 99 Total (i.e. Null);  97 Residual
    ## Null Deviance:       138.6 
    ## Residual Deviance: 136.1     AIC: 142.1
    ## ### Outcome model
    ## Call:
    ## survival::survreg(formula = Surv(y, event) ~ x * m + c, data = data, 
    ##     dist = "weibull")
    ## 
    ## Coefficients:
    ## (Intercept)           x           m           c         x:m 
    ## -1.04244118  0.44075656  0.09053705 -0.06689165  0.10031424 
    ## 
    ## Scale= 0.9658808 
    ## 
    ## Loglik(model)= -11.4   Loglik(intercept only)= -14.5
    ##  Chisq= 6.31 on 4 degrees of freedom, p= 0.177 
    ## n= 100 
    ## ### Mediation analysis 
    ##         cde        pnde        tnie        tnde        pnie          te          pm 
    ## 0.541070807 0.488930417 0.018240025 0.498503455 0.008666987 0.507170442 0.045436278

`summary()` to examine extended results
---------------------------------------

The `summary` method gives the summary for `mreg`, `yreg`, and mediation
analysis results. The `exponentiate` option will add the exponentiated
estimate and confidence limits if the outcome model is not a linear
model. The pure natural direct effect (`pnde`) is what is typically
called the natural direct effect (NDE). The total natural indirect
effect (`tnie`) is the corresponding natural indirect effect (NIE).

    summary(regmedint_obj, exponentiate = TRUE)

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
    ## survival::survreg(formula = Surv(y, event) ~ x * m + c, data = data, 
    ##     dist = "weibull")
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
    ##              est         se         Z           p       lower      upper exp(est) exp(lower) exp(upper)
    ## cde  0.541070807 0.29422958 1.8389409 0.065923882 -0.03560858 1.11775019 1.717845  0.9650179   3.057967
    ## pnde 0.488930417 0.21049248 2.3227928 0.020190284  0.07637274 0.90148809 1.630571  1.0793648   2.463266
    ## tnie 0.018240025 0.03706111 0.4921608 0.622605663 -0.05439841 0.09087846 1.018407  0.9470547   1.095136
    ## tnde 0.498503455 0.21209540 2.3503737 0.018754573  0.08280410 0.91420281 1.646256  1.0863290   2.494786
    ## pnie 0.008666987 0.02730994 0.3173565 0.750973092 -0.04485951 0.06219348 1.008705  0.9561318   1.064168
    ## te   0.507170442 0.21090051 2.4047853 0.016181972  0.09381303 0.92052785 1.660586  1.0983544   2.510615
    ## pm   0.045436278 0.01725694 2.6329276 0.008465238  0.01161329 0.07925926       NA         NA         NA
    ## 
    ## Evaluated at:
    ## avar: x
    ##  a1 (intervened value of avar) = 1
    ##  a0 (reference value of avar)  = 0
    ## mvar: m
    ##  m_cde (intervend value of mvar for cde) = 1
    ## cvar: c
    ##  c_cond (covariate vector value) = 0.5
    ## 
    ## Note that effect estimates can vary over m_cde and c_cond values when interaction = TRUE.

Use `coef` to extract the mediation analysis results only.

    coef(summary(regmedint_obj, exponentiate = TRUE))

    ##              est         se         Z           p       lower      upper exp(est) exp(lower) exp(upper)
    ## cde  0.541070807 0.29422958 1.8389409 0.065923882 -0.03560858 1.11775019 1.717845  0.9650179   3.057967
    ## pnde 0.488930417 0.21049248 2.3227928 0.020190284  0.07637274 0.90148809 1.630571  1.0793648   2.463266
    ## tnie 0.018240025 0.03706111 0.4921608 0.622605663 -0.05439841 0.09087846 1.018407  0.9470547   1.095136
    ## tnde 0.498503455 0.21209540 2.3503737 0.018754573  0.08280410 0.91420281 1.646256  1.0863290   2.494786
    ## pnie 0.008666987 0.02730994 0.3173565 0.750973092 -0.04485951 0.06219348 1.008705  0.9561318   1.064168
    ## te   0.507170442 0.21090051 2.4047853 0.016181972  0.09381303 0.92052785 1.660586  1.0983544   2.510615
    ## pm   0.045436278 0.01725694 2.6329276 0.008465238  0.01161329 0.07925926       NA         NA         NA

Note that the estimates can be re-evaluated at different `m_cde` and
`c_cond` without re-fitting the underlyng models.

    coef(summary(regmedint_obj, exponentiate = TRUE, m_cde = 0, c_cond = 5))

    ##              est         se         Z          p        lower      upper exp(est) exp(lower) exp(upper)
    ## cde  0.440756562 0.30083077 1.4651313 0.14288511 -0.148860903 1.03037403 1.553882  0.8616890   2.802114
    ## pnde 0.516631282 0.23499226 2.1985034 0.02791325  0.056054918 0.97720765 1.676371  1.0576558   2.657027
    ## tnie 0.012479419 0.02489241 0.5013344 0.61613580 -0.036308799 0.06126764 1.012558  0.9643425   1.063183
    ## tnde 0.523024107 0.24782597 2.1104492 0.03481969  0.037294136 1.00875408 1.687122  1.0379983   2.742182
    ## pnie 0.006086594 0.01889371 0.3221493 0.74733960 -0.030944391 0.04311758 1.006105  0.9695295   1.044061
    ## te   0.529110701 0.24194191 2.1869328 0.02874743  0.054913270 1.00330813 1.697422  1.0564490   2.727289
    ## pm   0.030184325 0.01276593 2.3644440 0.01805716  0.005163563 0.05520509       NA         NA         NA

Implementation status
=====================

This package is currently under development. The eAppendix for SAS macro
for causal mediation analysis with survival data (Valeri & VanderWeele
2015) describes the following grid of models.

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">

<table>
<thead>
<tr class="header">
<th style="text-align: left;">yreg  mreg</th>
<th style="text-align: center;">linear</th>
<th style="text-align: center;">logistic</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">linear</td>
<td style="text-align: center;">3</td>
<td style="text-align: center;">3</td>
</tr>
<tr class="even">
<td style="text-align: left;">logistic</td>
<td style="text-align: center;">3</td>
<td style="text-align: center;">3</td>
</tr>
<tr class="odd">
<td style="text-align: left;">loglinear</td>
<td style="text-align: center;">1</td>
<td style="text-align: center;">1</td>
</tr>
<tr class="even">
<td style="text-align: left;">poisson</td>
<td style="text-align: center;">3</td>
<td style="text-align: center;">3</td>
</tr>
<tr class="odd">
<td style="text-align: left;">negbin</td>
<td style="text-align: center;">3</td>
<td style="text-align: center;">3</td>
</tr>
<tr class="even">
<td style="text-align: left;">survCox</td>
<td style="text-align: center;">3</td>
<td style="text-align: center;">3</td>
</tr>
<tr class="odd">
<td style="text-align: left;">survAFT exp</td>
<td style="text-align: center;">3</td>
<td style="text-align: center;">3</td>
</tr>
<tr class="even">
<td style="text-align: left;">survAFT weibull</td>
<td style="text-align: center;">3</td>
<td style="text-align: center;">3</td>
</tr>
</tbody>
</table>

The current implementation status: 0. None (left blank) 1. SAS reference
results created 2. R implementation ongoing 3. R implementation complete

Formulas
========

The point estimate and standard error formulas (multivariate delta
method) were taken from the following references.

-   V2015: VanderWeele (2015) Explanation in Causal Inference.
-   VV2013A: Valeri & VanderWeele (2013) Appendix
-   VV2015A: Valeri & VanderWeele (2015) Appendix

As seen below, there are only four patterns.

Effect formulas
---------------

<table>
<colgroup>
<col style="width: 18%" />
<col style="width: 40%" />
<col style="width: 40%" />
</colgroup>
<thead>
<tr class="header">
<th>yreg  mreg</th>
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

Standard error formulas
-----------------------

<table>
<colgroup>
<col style="width: 20%" />
<col style="width: 39%" />
<col style="width: 39%" />
</colgroup>
<thead>
<tr class="header">
<th>yreg  mreg</th>
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

Design
======

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

Similar or related R projects
=============================

-   mediation:
    <a href="https:/cran.r-project.org/web/packages/mediation/index.html/" class="uri">https:/cran.r-project.org/web/packages/mediation/index.html/</a>
-   medflex:
    <a href="https:/cran.r-project.org/web/packages/medflex/index.html/" class="uri">https:/cran.r-project.org/web/packages/medflex/index.html/</a>
-   intmed:
    <a href="https:/cran.r-project.org/web/packages/intmed/index.html/" class="uri">https:/cran.r-project.org/web/packages/intmed/index.html/</a>
-   mediator:
    <a href="https:/github.com/GerkeLab/mediator/" class="uri">https:/github.com/GerkeLab/mediator/</a>
-   causalMediation:
    <a href="https:/github.com/harvard-P01/causalMediation" class="uri">https:/github.com/harvard-P01/causalMediation</a>
