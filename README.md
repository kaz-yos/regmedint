regmedint (developmental repo) <img src="hex.png" align="right" height="139"/>
==============================================================================

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

    ## install.packages("devtools") # If you do not have it.
    ## devtools::install_github("kaz-yos/regmedint") # If you have not installed it, yet.
    library(regmedint)
    library(tidyverse)

    ## FIXME: Find a meaningful data example within R.

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
