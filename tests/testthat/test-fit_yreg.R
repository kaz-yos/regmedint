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
### Internal function for yreg model fitting (linear)
################################################################################

test_that("fit_yreg fit linear models with lm correctly", {

    data(pbc)

    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L))

    ## No covariates
    yreg_linear_fit0 <- fit_yreg(yreg = "linear",
                                 data = pbc_cc,
                                 yvar = "alk.phos",
                                 avar = "trt",
                                 mvar = "bili",
                                 cvar = NULL,
                                 interaction = FALSE,
                                 eventvar = NULL)
    lm0 <- lm(formula = alk.phos ~ trt + bili,
              data = pbc_cc)
    ## Same classes
    expect_equal(class(yreg_linear_fit0),
                 class(lm0))
    ## Same formula
    expect_equal(as.character(yreg_linear_fit0$call$formula),
                 as.character(lm0$call$formula))
    ## Same coef
    expect_equal(coef(yreg_linear_fit0),
                 coef(lm0))
    ## Same vcov
    expect_equal(vcov(yreg_linear_fit0),
                 vcov(lm0))

    ## One covariates
    yreg_linear_fit1 <- fit_yreg(yreg = "linear",
                                 data = pbc_cc,
                                 yvar = "alk.phos",
                                 avar = "trt",
                                 mvar = "bili",
                                 cvar = c("age"),
                                 interaction = FALSE,
                                 eventvar = NULL)
    lm1 <- lm(formula = alk.phos ~ trt + bili + age,
              data = pbc_cc)
    ## Same classes
    expect_equal(class(yreg_linear_fit1),
                 class(lm1))
    ## Same formula
    expect_equal(as.character(yreg_linear_fit1$call$formula),
                 as.character(lm1$call$formula))
    ## Same coef
    expect_equal(coef(yreg_linear_fit1),
                 coef(lm1))
    ## Same vcov
    expect_equal(vcov(yreg_linear_fit1),
                 vcov(lm1))

    ## Three covariates
    yreg_linear_fit3 <- fit_yreg(yreg = "linear",
                                 data = pbc_cc,
                                 yvar = "alk.phos",
                                 avar = "trt",
                                 mvar = "bili",
                                 cvar = c("age","male","stage"),
                                 interaction = FALSE,
                                 eventvar = NULL)
    lm3 <- lm(formula = alk.phos ~ trt + bili + age + male + stage,
              data = pbc_cc)
    ## Same classes
    expect_equal(class(yreg_linear_fit3),
                 class(lm3))
    ## Same formula
    expect_equal(as.character(yreg_linear_fit3$call$formula),
                 as.character(lm3$call$formula))
    ## Same coef
    expect_equal(coef(yreg_linear_fit3),
                 coef(lm3))
    ## Same vcov
    expect_equal(vcov(yreg_linear_fit3),
                 vcov(lm3))

})


test_that("fit_yreg fit linear models with lm correctly with interaction", {

    data(pbc)

    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L))

    ## No covariates
    yreg_linear_fit0 <- fit_yreg(yreg = "linear",
                                 data = pbc_cc,
                                 yvar = "alk.phos",
                                 avar = "trt",
                                 mvar = "bili",
                                 cvar = NULL,
                                 interaction = TRUE,
                                 eventvar = NULL)
    lm0 <- lm(formula = alk.phos ~ trt*bili,
              data = pbc_cc)
    ## Same classes
    expect_equal(class(yreg_linear_fit0),
                 class(lm0))
    ## Same formula
    expect_equal(as.character(yreg_linear_fit0$call$formula),
                 as.character(lm0$call$formula))
    ## Same coef
    expect_equal(coef(yreg_linear_fit0),
                 coef(lm0))
    ## Same vcov
    expect_equal(vcov(yreg_linear_fit0),
                 vcov(lm0))

    ## One covariates
    yreg_linear_fit1 <- fit_yreg(yreg = "linear",
                                 data = pbc_cc,
                                 yvar = "alk.phos",
                                 avar = "trt",
                                 mvar = "bili",
                                 cvar = c("age"),
                                 interaction = TRUE,
                                 eventvar = NULL)
    lm1 <- lm(formula = alk.phos ~ trt*bili + age,
              data = pbc_cc)
    ## Same classes
    expect_equal(class(yreg_linear_fit1),
                 class(lm1))
    ## Same formula
    expect_equal(as.character(yreg_linear_fit1$call$formula),
                 as.character(lm1$call$formula))
    ## Same coef
    expect_equal(coef(yreg_linear_fit1),
                 coef(lm1))
    ## Same vcov
    expect_equal(vcov(yreg_linear_fit1),
                 vcov(lm1))

    ## Three covariates
    yreg_linear_fit3 <- fit_yreg(yreg = "linear",
                                 data = pbc_cc,
                                 yvar = "alk.phos",
                                 avar = "trt",
                                 mvar = "bili",
                                 cvar = c("age","male","stage"),
                                 interaction = TRUE,
                                 eventvar = NULL)
    lm3 <- lm(formula = alk.phos ~ trt*bili + age + male + stage,
              data = pbc_cc)
    ## Same classes
    expect_equal(class(yreg_linear_fit3),
                 class(lm3))
    ## Same formula
    expect_equal(as.character(yreg_linear_fit3$call$formula),
                 as.character(lm3$call$formula))
    ## Same coef
    expect_equal(coef(yreg_linear_fit3),
                 coef(lm3))
    ## Same vcov
    expect_equal(vcov(yreg_linear_fit3),
                 vcov(lm3))

})


###
### Internal function for yreg model fitting (logistic)
################################################################################

test_that("fit_yreg fit logistic models with glm correctly", {

    data(pbc)

    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L))

    ## No covariates
    yreg_linear_fit0 <- fit_yreg(yreg = "logistic",
                                 data = pbc_cc,
                                 yvar = "spiders",
                                 avar = "trt",
                                 mvar = "bili",
                                 cvar = NULL,
                                 interaction = FALSE,
                                 eventvar = NULL)
    lm0 <- lm(formula = spiders ~ trt + bili,
              data = pbc_cc)
    ## Same classes
    expect_equal(class(yreg_linear_fit0),
                 class(lm0))
    ## Same formula
    expect_equal(as.character(yreg_linear_fit0$call$formula),
                 as.character(lm0$call$formula))
    ## Same coef
    expect_equal(coef(yreg_linear_fit0),
                 coef(lm0))
    ## Same vcov
    expect_equal(vcov(yreg_linear_fit0),
                 vcov(lm0))

    ## One covariates
    yreg_linear_fit1 <- fit_yreg(yreg = "logistic",
                                 data = pbc_cc,
                                 yvar = "spiders",
                                 avar = "trt",
                                 mvar = "bili",
                                 cvar = c("age"),
                                 interaction = FALSE,
                                 eventvar = NULL)
    lm1 <- lm(formula = spiders ~ trt + bili + age,
              data = pbc_cc)
    ## Same classes
    expect_equal(class(yreg_linear_fit1),
                 class(lm1))
    ## Same formula
    expect_equal(as.character(yreg_linear_fit1$call$formula),
                 as.character(lm1$call$formula))
    ## Same coef
    expect_equal(coef(yreg_linear_fit1),
                 coef(lm1))
    ## Same vcov
    expect_equal(vcov(yreg_linear_fit1),
                 vcov(lm1))

    ## Three covariates
    yreg_linear_fit3 <- fit_yreg(yreg = "logistic",
                                 data = pbc_cc,
                                 yvar = "spiders",
                                 avar = "trt",
                                 mvar = "bili",
                                 cvar = c("age","male","stage"),
                                 interaction = FALSE,
                                 eventvar = NULL)
    lm3 <- lm(formula = spiders ~ trt + bili + age + male + stage,
              data = pbc_cc)
    ## Same classes
    expect_equal(class(yreg_linear_fit3),
                 class(lm3))
    ## Same formula
    expect_equal(as.character(yreg_linear_fit3$call$formula),
                 as.character(lm3$call$formula))
    ## Same coef
    expect_equal(coef(yreg_linear_fit3),
                 coef(lm3))
    ## Same vcov
    expect_equal(vcov(yreg_linear_fit3),
                 vcov(lm3))

})


test_that("fit_yreg fit logistic models with glm correctly with interaction", {

    data(pbc)

    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L))

    ## No covariates
    yreg_linear_fit0 <- fit_yreg(yreg = "logistic",
                                 data = pbc_cc,
                                 yvar = "spiders",
                                 avar = "trt",
                                 mvar = "bili",
                                 cvar = NULL,
                                 interaction = TRUE,
                                 eventvar = NULL)
    lm0 <- lm(formula = spiders ~ trt*bili,
              data = pbc_cc)
    ## Same classes
    expect_equal(class(yreg_linear_fit0),
                 class(lm0))
    ## Same formula
    expect_equal(as.character(yreg_linear_fit0$call$formula),
                 as.character(lm0$call$formula))
    ## Same coef
    expect_equal(coef(yreg_linear_fit0),
                 coef(lm0))
    ## Same vcov
    expect_equal(vcov(yreg_linear_fit0),
                 vcov(lm0))

    ## One covariates
    yreg_linear_fit1 <- fit_yreg(yreg = "logistic",
                                 data = pbc_cc,
                                 yvar = "spiders",
                                 avar = "trt",
                                 mvar = "bili",
                                 cvar = c("age"),
                                 interaction = TRUE,
                                 eventvar = NULL)
    lm1 <- lm(formula = spiders ~ trt*bili + age,
              data = pbc_cc)
    ## Same classes
    expect_equal(class(yreg_linear_fit1),
                 class(lm1))
    ## Same formula
    expect_equal(as.character(yreg_linear_fit1$call$formula),
                 as.character(lm1$call$formula))
    ## Same coef
    expect_equal(coef(yreg_linear_fit1),
                 coef(lm1))
    ## Same vcov
    expect_equal(vcov(yreg_linear_fit1),
                 vcov(lm1))

    ## Three covariates
    yreg_linear_fit3 <- fit_yreg(yreg = "logistic",
                                 data = pbc_cc,
                                 yvar = "spiders",
                                 avar = "trt",
                                 mvar = "bili",
                                 cvar = c("age","male","stage"),
                                 interaction = TRUE,
                                 eventvar = NULL)
    lm3 <- lm(formula = spiders ~ trt*bili + age + male + stage,
              data = pbc_cc)
    ## Same classes
    expect_equal(class(yreg_linear_fit3),
                 class(lm3))
    ## Same formula
    expect_equal(as.character(yreg_linear_fit3$call$formula),
                 as.character(lm3$call$formula))
    ## Same coef
    expect_equal(coef(yreg_linear_fit3),
                 coef(lm3))
    ## Same vcov
    expect_equal(vcov(yreg_linear_fit3),
                 vcov(lm3))

})
