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

test_that("fit_yreg fit logistic models with glm correctly", {

    data(pbc)

    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L))

    ## No covariates
    yreg_fit0 <- fit_yreg(yreg = "logistic",
                          data = pbc_cc,
                          yvar = "spiders",
                          avar = "trt",
                          mvar = "bili",
                          cvar = NULL,
                          interaction = FALSE,
                          eventvar = NULL)
    glm0 <- glm(formula = spiders ~ trt + bili,
                family = binomial(link = "logit"),
                data = pbc_cc)
    ## Same classes
    expect_equal(class(yreg_fit0),
                 class(glm0))
    ## Same formula
    expect_equal(as.character(yreg_fit0$call$formula),
                 as.character(glm0$call$formula))
    ## Same coef
    expect_equal(coef(yreg_fit0),
                 coef(glm0))
    ## Same vcov
    expect_equal(vcov(yreg_fit0),
                 vcov(glm0))

    ## One covariates
    yreg_fit1 <- fit_yreg(yreg = "logistic",
                          data = pbc_cc,
                          yvar = "spiders",
                          avar = "trt",
                          mvar = "bili",
                          cvar = c("age"),
                          interaction = FALSE,
                          eventvar = NULL)
    glm1 <- glm(formula = spiders ~ trt + bili + age,
                family = binomial(link = "logit"),
                data = pbc_cc)
    ## Same classes
    expect_equal(class(yreg_fit1),
                 class(glm1))
    ## Same formula
    expect_equal(as.character(yreg_fit1$call$formula),
                 as.character(glm1$call$formula))
    ## Same coef
    expect_equal(coef(yreg_fit1),
                 coef(glm1))
    ## Same vcov
    expect_equal(vcov(yreg_fit1),
                 vcov(glm1))

    ## Three covariates
    yreg_fit3 <- fit_yreg(yreg = "logistic",
                          data = pbc_cc,
                          yvar = "spiders",
                          avar = "trt",
                          mvar = "bili",
                          cvar = c("age","male","stage"),
                          interaction = FALSE,
                          eventvar = NULL)
    glm3 <- glm(formula = spiders ~ trt + bili + age + male + stage,
                family = binomial(link = "logit"),
                data = pbc_cc)
    ## Same classes
    expect_equal(class(yreg_fit3),
                 class(glm3))
    ## Same formula
    expect_equal(as.character(yreg_fit3$call$formula),
                 as.character(glm3$call$formula))
    ## Same coef
    expect_equal(coef(yreg_fit3),
                 coef(glm3))
    ## Same vcov
    expect_equal(vcov(yreg_fit3),
                 vcov(glm3))

})


test_that("fit_yreg fit logistic models with glm correctly with interaction", {

    data(pbc)

    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L))

    ## No covariates
    yreg_fit0 <- fit_yreg(yreg = "logistic",
                          data = pbc_cc,
                          yvar = "spiders",
                          avar = "trt",
                          mvar = "bili",
                          cvar = NULL,
                          interaction = TRUE,
                          eventvar = NULL)
    glm0 <- glm(formula = spiders ~ trt*bili,
                family = binomial(link = "logit"),
                data = pbc_cc)
    ## Same classes
    expect_equal(class(yreg_fit0),
                 class(glm0))
    ## Same formula
    expect_equal(as.character(yreg_fit0$call$formula),
                 as.character(glm0$call$formula))
    ## Same coef
    expect_equal(coef(yreg_fit0),
                 coef(glm0))
    ## Same vcov
    expect_equal(vcov(yreg_fit0),
                 vcov(glm0))

    ## One covariates
    yreg_fit1 <- fit_yreg(yreg = "logistic",
                          data = pbc_cc,
                          yvar = "spiders",
                          avar = "trt",
                          mvar = "bili",
                          cvar = c("age"),
                          interaction = TRUE,
                          eventvar = NULL)
    glm1 <- glm(formula = spiders ~ trt*bili + age,
                family = binomial(link = "logit"),
                data = pbc_cc)
    ## Same classes
    expect_equal(class(yreg_fit1),
                 class(glm1))
    ## Same formula
    expect_equal(as.character(yreg_fit1$call$formula),
                 as.character(glm1$call$formula))
    ## Same coef
    expect_equal(coef(yreg_fit1),
                 coef(glm1))
    ## Same vcov
    expect_equal(vcov(yreg_fit1),
                 vcov(glm1))

    ## Three covariates
    yreg_fit3 <- fit_yreg(yreg = "logistic",
                          data = pbc_cc,
                          yvar = "spiders",
                          avar = "trt",
                          mvar = "bili",
                          cvar = c("age","male","stage"),
                          interaction = TRUE,
                          eventvar = NULL)
    glm3 <- glm(formula = spiders ~ trt*bili + age + male + stage,
                family = binomial(link = "logit"),
                data = pbc_cc)
    ## Same classes
    expect_equal(class(yreg_fit3),
                 class(glm3))
    ## Same formula
    expect_equal(as.character(yreg_fit3$call$formula),
                 as.character(glm3$call$formula))
    ## Same coef
    expect_equal(coef(yreg_fit3),
                 coef(glm3))
    ## Same vcov
    expect_equal(vcov(yreg_fit3),
                 vcov(glm3))

})
