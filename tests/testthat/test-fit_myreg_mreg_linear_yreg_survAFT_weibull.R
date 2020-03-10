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
### Internal function for yreg model fitting (logistic)
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
    myreg_fit0 <- fit_myreg()


    ## One covariates
    yreg_fit1 <- fit_yreg(yreg = "survAFT_weibull",
                          data = pbc_cc,
                          yvar = "time",
                          avar = "trt",
                          mvar = "bili",
                          cvar = c("age"),
                          interaction = FALSE,
                          eventvar = "status")
    ref_fit1 <- survreg(formula = Surv(time,status) ~ trt + bili + age,
                        dist = "weibull",
                        data = pbc_cc)
    ## Same classes
    expect_equal(class(yreg_fit1),
                 class(ref_fit1))
    ## Same formula
    expect_equal(as.character(yreg_fit1$call$formula),
                 as.character(ref_fit1$call$formula))
    ## Same coef
    expect_equal(coef(yreg_fit1),
                 coef(ref_fit1))
    ## Same vcov
    expect_equal(vcov(yreg_fit1),
                 vcov(ref_fit1))

    ## Three covariates
    yreg_fit3 <- fit_yreg(yreg = "survAFT_weibull",
                          data = pbc_cc,
                          yvar = "time",
                          avar = "trt",
                          mvar = "bili",
                          cvar = c("age","male","stage"),
                          interaction = FALSE,
                          eventvar = "status")
    ref_fit3 <- survreg(formula = Surv(time,status) ~ trt + bili + age + male + stage,
                        dist = "weibull",
                        data = pbc_cc)
    ## Same classes
    expect_equal(class(yreg_fit3),
                 class(ref_fit3))
    ## Same formula
    expect_equal(as.character(yreg_fit3$call$formula),
                 as.character(ref_fit3$call$formula))
    ## Same coef
    expect_equal(coef(yreg_fit3),
                 coef(ref_fit3))
    ## Same vcov
    expect_equal(vcov(yreg_fit3),
                 vcov(ref_fit3))

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
    ref_fit0 <- survreg(formula = Surv(time,status) ~ trt*bili,
                        dist = "weibull",
                        data = pbc_cc)
    ## Same classes
    expect_equal(class(yreg_fit0),
                 class(ref_fit0))
    ## Same formula
    expect_equal(as.character(yreg_fit0$call$formula),
                 as.character(ref_fit0$call$formula))
    ## Same coef
    expect_equal(coef(yreg_fit0),
                 coef(ref_fit0))
    ## Same vcov
    expect_equal(vcov(yreg_fit0),
                 vcov(ref_fit0))

    ## One covariates
    yreg_fit1 <- fit_yreg(yreg = "survAFT_weibull",
                          data = pbc_cc,
                          yvar = "time",
                          avar = "trt",
                          mvar = "bili",
                          cvar = c("age"),
                          interaction = TRUE,
                          eventvar = "status")
    ref_fit1 <- survreg(formula = Surv(time,status) ~ trt*bili + age,
                        dist = "weibull",
                        data = pbc_cc)
    ## Same classes
    expect_equal(class(yreg_fit1),
                 class(ref_fit1))
    ## Same formula
    expect_equal(as.character(yreg_fit1$call$formula),
                 as.character(ref_fit1$call$formula))
    ## Same coef
    expect_equal(coef(yreg_fit1),
                 coef(ref_fit1))
    ## Same vcov
    expect_equal(vcov(yreg_fit1),
                 vcov(ref_fit1))

    ## Three covariates
    yreg_fit3 <- fit_yreg(yreg = "survAFT_weibull",
                          data = pbc_cc,
                          yvar = "time",
                          avar = "trt",
                          mvar = "bili",
                          cvar = c("age","male","stage"),
                          interaction = TRUE,
                          eventvar = "status")
    ref_fit3 <- survreg(formula = Surv(time,status) ~ trt*bili + age + male + stage,
                        dist = "weibull",
                        data = pbc_cc)
    ## Same classes
    expect_equal(class(yreg_fit3),
                 class(ref_fit3))
    ## Same formula
    expect_equal(as.character(yreg_fit3$call$formula),
                 as.character(ref_fit3$call$formula))
    ## Same coef
    expect_equal(coef(yreg_fit3),
                 coef(ref_fit3))
    ## Same vcov
    expect_equal(vcov(yreg_fit3),
                 vcov(ref_fit3))

})
