regmedint:
================
Kazuki Yoshida
2020-03-29

# regmedint (developmental repo) <img src="man/figures/hex.png" align="right" height="140"/>

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

# Example

``` r
## install.packages("devtools") # If you do not have it.
## devtools::install_github("kaz-yos/regmedint") # If you have not installed it, yet.
library(regmedint)
library(tidyverse)

## FIXME: Find a meaningful data example within R.
regmedint_obj <- regmedint(data,
                           ##
                           yvar,
                           avar,
                           mvar,
                           cvar,
                           ##
                           a0,
                           a1,
                           m_cde,
                           c_cond,
                           ##
                           mreg,
                           yreg,
                           ##
                           interaction = TRUE,
                           casecontrol = FALSE,
                           eventvar = NULL)

## Implicit printing prints mreg, yreg, and mediation analysis point estimates.
regmedint_obj
## A fuller result is obtained.
summary(regmedint_obj)
## Add exponentiated results.
summary(regmedint_obj, exponentiate = TRUE)
## The estimates can be re-evaluated without model refitting at a different mvar value.
summary(regmedint_obj, m_cde = )
## The estimates can be re-evaluated without model refitting at a different c_cond value.
summary(regmedint_obj, c_cond = )
```

# Implementation status

This package is currently under development. The eAppendix for SAS macro
for causal mediation analysis with survival data (Valeri & VanderWeele
2015) describes the following grid of
models.

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">

| yreg  mreg      | linear | logistic |
| :-------------- | :----: | :------: |
| linear          |   3    |    3     |
| logistic        |   3    |    3     |
| loglinear       |   1    |    1     |
| poisson         |   3    |    3     |
| negbin          |   3    |    3     |
| survCox         |   3    |    3     |
| survAFT exp     |   3    |    3     |
| survAFT weibull |   3    |    3     |

The current implementation status: 0. None (left blank) 1. SAS reference
results created 2. R implementation ongoing 3. R implementation complete

# Formulas

The point estimate and standard error formulas (multivariate delta
method) were taken from the following references.

  - V2015: VanderWeele (2015) Explanation in Causal Inference.
  - VV2013A: Valeri & VanderWeele (2013) Appendix
  - VV2015A: Valeri & VanderWeele (2015) Appendix

As seen below, there are only four
patterns.

## Effect formulas

| yreg  mreg      | linear                               | logistic                             |
| --------------- | ------------------------------------ | ------------------------------------ |
| linear          | V2015 p466 Proposition 2.3           | V2015 p471 Proposition 2.5           |
|                 |                                      |                                      |
| logistic        | V2015 p468 Proposition 2.4           | V2015 p473 Proposition 2.6           |
| loglinear       | VV2013A p8 Use Proposition 2.4       | VV2013A p8 Use Proposition 2.6       |
| poisson         | VV2013A p8 Use Proposition 2.4       | VV2013A p8 Use Proposition 2.6       |
| negbin          | VV2013A p8 Use Proposition 2.4       | VV2013A p8 Use Proposition 2.6       |
|                 |                                      |                                      |
| survCox         | V2015 p496 Proposition 4.4 (Use 2.4) | V2015 p499 Proposition 4.6 (Use 2.6) |
| survAFT exp     | V2015 p494 Proposition 4.1 (Use 2.4) | V2015 p495 Proposition 4.3 (Use 2.6) |
| survAFT weibull | V2015 p494 Proposition 4.1 (Use 2.4) | V2015 p495 Proposition 4.3 (Use 2.6) |

Note the loglinear outcome model means log link with binomial
error.

## Standard error formulas

| yreg  mreg      | linear                         | logistic                       |
| --------------- | ------------------------------ | ------------------------------ |
| linear          | V2015 p466 Proposition 2.3     | V2015 p471 Proposition 2.5     |
|                 |                                |                                |
| logistic        | V2015 p468 Proposition 2.4     | V2015 p473 Proposition 2.6     |
| loglinear       | VV2013A p8 Use Proposition 2.4 | VV2013A p8 Use Proposition 2.6 |
| poisson         | VV2013A p8 Use Proposition 2.4 | VV2013A p8 Use Proposition 2.6 |
| negbin          | VV2013A p8 Use Proposition 2.4 | VV2013A p8 Use Proposition 2.6 |
|                 |                                |                                |
| survCox         | V2015 p496 Use Proposition 2.4 | V2015 p499 Use Proposition 2.6 |
| survAFT exp     | V2015 p494 Use Proposition 2.4 | V2015 p495 Use Proposition 2.6 |
| survAFT weibull | V2015 p494 Use Proposition 2.4 | V2015 p495 Use Proposition 2.6 |

# Design

The software design is outlined here for those who may be interested.

  - Call structure
      - regmedint UI function
          - new\_regmedint internal constructor
              - fit\_mreg
              - fit\_yreg
              - calc\_myreg calls a specialized worker function, which
                return two functions, one for point estimates and the
                other for standard error estimate.
                  - calc\_myreg\_mreg\_linear\_yreg\_linear
                  - calc\_myreg\_mreg\_linear\_yreg\_logistic
                  - calc\_myreg\_mreg\_logistic\_yreg\_linear
                  - calc\_myreg\_mreg\_logistic\_yreg\_logistic
  - regmedint object structure
      - mreg\_fit mediator regression model object as is
      - yreg\_fit outcome regression model object as is
      - myreg\_funs list
          - est\_fun: (a0,a1,m\_cde,c\_cond) →
            (cde,pnde,tnie,tnde,pnie,te,pm)
          - se\_fun: (a0,a1,m\_cde,c\_cond) → se for
            (cde,pnde,tnie,tnde,pnie,te,pm)
          - args preserves arguments given to the UI
  - User methods for the regmedint object
      - print.regmedint: prints coefficients for mreg, yreg, and
        mediation analysis
      - summary.regmedint: regmedint → summary\_regmedint
          - print.summary\_regmedint: prints summary objects for mreg,
            yreg, and mediation analysis
          - coef.summary\_regmedint:
      - coef.regmedint: regmedint → vector
        (cde,pnde,tnie,tnde,pnie,te,pm)
      - vcov.regmedint: regmedint → matrix
        (cde,pnde,tnie,tnde,pnie,te,pm). Off-diagonals are NA.
      - confint.regmedint: regmedint → matrix of (lower,upper)

# Similar or related R projects

  - mediation:
    <https:/cran.r-project.org/web/packages/mediation/index.html/>
  - medflex:
    <https:/cran.r-project.org/web/packages/medflex/index.html/>
  - intmed: <https:/cran.r-project.org/web/packages/intmed/index.html/>
  - mediator: <https:/github.com/GerkeLab/mediator/>
  - causalMediation: <https:/github.com/harvard-P01/causalMediation>
