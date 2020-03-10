################################################################################
### Tests for internal functions
##
## Created on: 2020-03-09
## Author: Kazuki Yoshida
################################################################################

## Load testthat in case this is run in isolation.
library(testthat)
library(survival)
library(tidyverse)


###
### Internal function for myreg model fitting
################################################################################

test_that("fit_myreg fit linear / Weibull AFT models correctly", {

    data(pbc)

    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L),
               ## Combine transplant and death for testing purpose
               status = if_else(status == 0, 0L, 1L))

    ## No covariates
    mreg_fit0 <- fit_mreg(mreg = "linear",
                          data = pbc_cc,
                          avar = "trt",
                          mvar = "bili",
                          cvar = NULL)
    yreg_fit0 <- fit_yreg(yreg = "survAFT_weibull",
                          data = pbc_cc,
                          yvar = "time",
                          avar = "trt",
                          mvar = "bili",
                          cvar = NULL,
                          interaction = FALSE,
                          eventvar = "status")
    myreg_fit0 <- fit_myreg(mreg_fit0, yreg_fit0)
    ref_fit0 <- tibble(beta_0 = coef(mreg_fit0)["(Intercept)"],
                       beta_A = coef(mreg_fit0)["trt"],
                       beta_C = list(0),
                       sigma = sigma(mreg_fit0),
                       theta_A = coef(yreg_fit0)["trt"],
                       theta_M = coef(yreg_fit0)["bili"],
                       theta_AM = 0,
                       a1 = 1,
                       a0 = 0,
                       m_cde = 0.6,
                       c_cond = list(1.1),
                       ## c part of linear predictor by inner product of beta_C and c_cond
                       clp = map2_dbl(beta_C, c_cond, function(a,b) {sum(a*b)})) %>%
        ## VanderWeele 2015 p494
        ## natural effects on log scale
        mutate(pnde = (theta_A +
                       theta_AM * ((beta_0 + beta_A * a0 + clp) +
                                   sigma^2 * (theta_M + 1/2 * theta_AM * (a1 + a0)))) * (a1 - a0),
               tnie = beta_A * (theta_M + theta_AM * a1) * (a1 - a0),
               cde_m = (theta_A + theta_AM * m_cde) * (a1 - a0),
               ## FIXME: Made upg
               tnde = 999,
               pnie = 888)
    expect_equal(myreg_fit0$myreg_fit["pnde"],
                 ref_fit0$pnde)
    expect_equal(myreg_fit0$myreg_fit["tnie"],
                 ref_fit0$tnie)
    expect_equal(myreg_fit0$myreg_fit["cde_m"],
                 ref_fit0$cde_m)
    expect_equal(myreg_fit0$myreg_fit["tnde"],
                 ref_fit0$tnde)
    expect_equal(myreg_fit0$myreg_fit["pnie"],
                 ref_fit0$pnie)

    ## One covariates
    yreg_fit1 <- fit_yreg(yreg = "survAFT_weibull",
                          data = pbc_cc,
                          yvar = "time",
                          avar = "trt",
                          mvar = "bili",
                          cvar = c("age"),
                          interaction = FALSE,
                          eventvar = "status")
    mreg_fit1 <- fit_mreg(mreg = "linear",
                          data = pbc_cc,
                          avar = "trt",
                          mvar = "bili",
                          cvar = c("age"))
    expect_equal(TRUE,
                 FALSE)

    ## Three covariates
    yreg_fit3 <- fit_yreg(yreg = "survAFT_weibull",
                          data = pbc_cc,
                          yvar = "time",
                          avar = "trt",
                          mvar = "bili",
                          cvar = c("age","male","stage"),
                          interaction = FALSE,
                          eventvar = "status")
    mreg_fit3 <- fit_mreg(mreg = "linear",
                          data = pbc_cc,
                          avar = "trt",
                          mvar = "bili",
                          cvar = c("age","male","stage"))
    expect_equal(TRUE,
                 FALSE)

})


test_that("fit_myreg fit linear / Weibull AFT models correctly with interaction", {

    data(pbc)

    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L),
               ## Combine transplant and death for testing purpose
               status = if_else(status == 0, 0L, 1L))

    ## No covariates
    yreg_fit0 <- fit_yreg(yreg = "survAFT_weibull",
                          data = pbc_cc,
                          yvar = "time",
                          avar = "trt",
                          mvar = "bili",
                          cvar = NULL,
                          interaction = TRUE,
                          eventvar = "status")
    expect_equal(TRUE,
                 FALSE)

    ## One covariates
    yreg_fit1 <- fit_yreg(yreg = "survAFT_weibull",
                          data = pbc_cc,
                          yvar = "time",
                          avar = "trt",
                          mvar = "bili",
                          cvar = c("age"),
                          interaction = TRUE,
                          eventvar = "status")
    expect_equal(TRUE,
                 FALSE)

    ## Three covariates
    yreg_fit3 <- fit_yreg(yreg = "survAFT_weibull",
                          data = pbc_cc,
                          yvar = "time",
                          avar = "trt",
                          mvar = "bili",
                          cvar = c("age","male","stage"),
                          interaction = TRUE,
                          eventvar = "status")
    expect_equal(TRUE,
                 FALSE)


})
