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

test_that("calc_myreg fit linear / Weibull AFT models correctly", {

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
    myreg_fit0 <- calc_myreg(mreg = "linear",
                             mreg_fit = mreg_fit0,
                             yreg = "survAFT_weibull",
                             yreg_fit = yreg_fit0,
                             interaction = FALSE)
    ## Point estimates
    ref_est0 <- tibble(beta_0 = coef(mreg_fit0)["(Intercept)"],
                       beta1 = coef(mreg_fit0)["trt"],
                       beta2 = list(0),
                       sigma = sigma(mreg_fit0),
                       theta1 = coef(yreg_fit0)["trt"],
                       theta2 = coef(yreg_fit0)["bili"],
                       theta3 = 0,
                       a1 = 1,
                       a0 = 0,
                       m_cde = 0.6,
                       c_cond = list(1.1),
                       ## c part of linear predictor by inner product of beta2 and c_cond
                       clp = map2_dbl(beta_C, c_cond, function(a,b) {sum(a*b)})) %>%
        ## VanderWeele 2015 p494
        ## natural effects on log scale
        mutate(pnde = (theta_A +
                       theta3 * ((beta_0 + beta1 * a0 + clp) +
                                   sigma^2 * (theta_M + 1/2 * theta3 * (a1 + a0)))) * (a1 - a0),
               tnie = beta1 * (theta_M + theta3 * a1) * (a1 - a0),
               cde_m = (theta_A + theta3 * m_cde) * (a1 - a0),
               ## FIXME: Made upg
               tnde = 999,
               pnie = 888)
    expect_equal(myreg_fit0$myreg_fit["pnde"],
                 ref_est0$pnde)
    expect_equal(myreg_fit0$myreg_fit["tnie"],
                 ref_est0$tnie)
    expect_equal(myreg_fit0$myreg_fit["cde_m"],
                 ref_est0$cde_m)
    expect_equal(myreg_fit0$myreg_fit["tnde"],
                 ref_est0$tnde)
    expect_equal(myreg_fit0$myreg_fit["pnie"],
                 ref_est0$pnie)
    ## Standard error estimates
    ref_se0 <- tibble(Sigma_beta = vcov(mreg_fit0),
                      ## FIXME: This contains dimension for the scale parameter?
                      Sigma_theta = vcov(yreg_fit0),
                      ## FIXME: Variance of sigma^2 estimate
                      Sigma_sigma = 0,
                      a1 = 1,
                      a0 = 0,
                      m_cde = 0.6,
                      c_cond = list(1.1),
                      ## c part of linear predictor by inner product of beta2 and c_cond
                      clp = map2_dbl(beta_C, c_cond, function(a,b) {sum(a*b)}),
                      ## VanderWeele 2015 p468
                      Gamma_cde_m =
                          list(c(0,0,
                                 0,1,0,m_cde,
                                 0)),
                      Gamma_pnde =
                          list(c(theta_AM,theta_AM*a0,
                                 0,1,theta_AM * sigma^2, beta_0 + beta1 * a0 + clp + theta2 * sigma^2 + theta3 * sigma^2 * (a1 + a0),
                                 theta3 * theta2 + 1/2 * theta3^2 * (a1 + a0))),
                      Gamma_tnie =
                          list(c(0,theta_M + theta3 * a,
                                 0,0,beta1,beta1 * a,
                                 0)))

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


test_that("calc_myreg fit linear / Weibull AFT models correctly with interaction", {

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
