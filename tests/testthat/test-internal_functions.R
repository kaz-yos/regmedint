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
### Internal function for mreg string formula creation
################################################################################

test_that("string_mreg_formula create sound formula strings", {

    expect_equal(string_mreg_formula("M","A",NULL),
                 "M ~ A")

    expect_equal(string_mreg_formula("M","A",c("C")),
                 "M ~ A + C")

    expect_equal(string_mreg_formula("M","A",c("C1","C2","C3")),
                 "M ~ A + C1 + C2 + C3")

})


###
### Internal function for yreg string formula creation
################################################################################

test_that("string_yreg_formula create sound formula strings for non-survival outcomes", {

    ## Zero covariates
    expect_equal(string_yreg_formula("Y","A","M",NULL, interaction = FALSE,
                                     eventvar = NULL),
                 "Y ~ A + M")
    expect_equal(string_yreg_formula("Y","A","M",NULL, interaction = TRUE,
                                     eventvar = NULL),
                 "Y ~ A*M")

    ## One covariate
    expect_equal(string_yreg_formula("Y","A","M",c("C"), interaction = FALSE,
                                     eventvar = NULL),
                 "Y ~ A + M + C")
    expect_equal(string_yreg_formula("Y","A","M",c("C"), interaction = TRUE,
                                     eventvar = NULL),
                 "Y ~ A*M + C")

    ## Three covariates
    expect_equal(string_yreg_formula("Y","A","M",c("C1","C2","C3"), interaction = FALSE,
                                     eventvar = NULL),
                 "Y ~ A + M + C1 + C2 + C3")
    expect_equal(string_yreg_formula("Y","A","M",c("C1","C2","C3"), interaction = TRUE,
                                     eventvar = NULL),
                 "Y ~ A*M + C1 + C2 + C3")
})


test_that("string_yreg_formula create sound formula strings for survival outcomes", {

    ## Zero covariates
    expect_equal(string_yreg_formula("time","A","M",NULL,"event", interaction = FALSE,
                                     eventvar = NULL),
                 "Surv(time, event) ~ A + M")
    expect_equal(string_yreg_formula("time","A","M",NULL,"event", interaction = TRUE,
                                     eventvar = NULL),
                 "Surv(time, event) ~ A*M")

    ## One covariate
    expect_equal(string_yreg_formula("time","A","M",c("C"),"event", interaction = FALSE,
                                     eventvar = NULL),
                 "Surv(time, event) ~ A + M + C")
    expect_equal(string_yreg_formula("time","A","M",c("C"),"event", interaction = TRUE,
                                     eventvar = NULL),
                 "Surv(time, event) ~ A*M + C")

    ## Three covariates
    expect_equal(string_yreg_formula("time","A","M",c("C1","C2","C3"), interaction = FALSE,
                                     eventvar = "event"),
                 "Surv(time, event) ~ A + M + C1 + C2 + C3")
    expect_equal(string_yreg_formula("time","A","M",c("C1","C2","C3"), interaction = TRUE,
                                     eventvar = "event"),
                 "Surv(time, event) ~ A*M + C1 + C2 + C3")
})


###
### Internal function for mreg model fitting
################################################################################

test_that("fit_mreg fit linear models with lm correctly", {

    data(pbc)

    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L))

    ## No covariates
    mreg_linear_fit0 <- fit_mreg(mreg = "linear",
                                 data = pbc_cc,
                                 avar = "trt",
                                 mvar = "bili",
                                 cvar = NULL)
    lm0 <- lm(formula = bili ~ trt,
              data = pbc_cc)
    ## Same classes
    expect_equal(class(mreg_linear_fit0),
                 class(lm0))
    ## Same formula
    expect_equal(as.character(mreg_linear_fit0$call$formula),
                 as.character(lm0$call$formula))
    ## Same coef
    expect_equal(coef(mreg_linear_fit0),
                 coef(lm0))
    ## Same vcov
    expect_equal(vcov(mreg_linear_fit0),
                 vcov(lm0))

    ## One covariates
    mreg_linear_fit1 <- fit_mreg(mreg = "linear",
                                 data = pbc_cc,
                                 avar = "trt",
                                 mvar = "bili",
                                 cvar = c("age"))
    lm1 <- lm(formula = bili ~ trt + age,
              data = pbc_cc)
    ## Same classes
    expect_equal(class(mreg_linear_fit1),
                 class(lm1))
    ## Same formula
    expect_equal(as.character(mreg_linear_fit1$call$formula),
                 as.character(lm1$call$formula))
    ## Same coef
    expect_equal(coef(mreg_linear_fit1),
                 coef(lm1))
    ## Same vcov
    expect_equal(vcov(mreg_linear_fit1),
                 vcov(lm1))

    ## Three covariates
    mreg_linear_fit3 <- fit_mreg(mreg = "linear",
                                 data = pbc_cc,
                                 avar = "trt",
                                 mvar = "bili",
                                 cvar = c("age","male","stage"))
    lm3 <- lm(formula = bili ~ trt + age + male + stage,
              data = pbc_cc)
    ## Same classes
    expect_equal(class(mreg_linear_fit3),
                 class(lm3))
    ## Same formula
    expect_equal(as.character(mreg_linear_fit3$call$formula),
                 as.character(lm3$call$formula))
    ## Same coef
    expect_equal(coef(mreg_linear_fit3),
                 coef(lm3))
    ## Same vcov
    expect_equal(vcov(mreg_linear_fit3),
                 vcov(lm3))

})

test_that("fit_mreg fit logistic models with glm correctly", {

    data(pbc)

    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L))

    ## No covariates
    mreg_logistic_fit0 <- fit_mreg(mreg = "logistic",
                                   data = pbc_cc,
                                   avar = "trt",
                                   mvar = "hepato",
                                   cvar = NULL)
    glm0 <- glm(formula = hepato ~ trt,
                family = binomial(link = "logit"),
                data = pbc_cc)
    ## Same classes
    expect_equal(class(mreg_logistic_fit0),
                 class(glm0))
    ## Same formula
    expect_equal(as.character(mreg_logistic_fit0$call$formula),
                 as.character(glm0$call$formula))
    ## Same coef
    expect_equal(coef(mreg_logistic_fit0),
                 coef(glm0))
    ## Same vcov
    expect_equal(vcov(mreg_logistic_fit0),
                 vcov(glm0))

    ## One covariates
    mreg_logistic_fit1 <- fit_mreg(mreg = "logistic",
                                   data = pbc_cc,
                                   avar = "trt",
                                   mvar = "hepato",
                                   cvar = c("age"))
    glm1 <- glm(formula = hepato ~ trt + age,
                family = binomial(link = "logit"),
                data = pbc_cc)
    ## Same classes
    expect_equal(class(mreg_logistic_fit1),
                 class(glm1))
    ## Same formula
    expect_equal(as.character(mreg_logistic_fit1$call$formula),
                 as.character(glm1$call$formula))
    ## Same coef
    expect_equal(coef(mreg_logistic_fit1),
                 coef(glm1))
    ## Same vcov
    expect_equal(vcov(mreg_logistic_fit1),
                 vcov(glm1))

    ## Three covariates
    mreg_logistic_fit3 <- fit_mreg(mreg = "logistic",
                                   data = pbc_cc,
                                   avar = "trt",
                                   mvar = "hepato",
                                   cvar = c("age","male","stage"))
    glm3 <- glm(formula = hepato ~ trt + age + male + stage,
                family = binomial(link = "logit"),
                data = pbc_cc)
    ## Same classes
    expect_equal(class(mreg_logistic_fit3),
                 class(glm3))
    ## Same formula
    expect_equal(as.character(mreg_logistic_fit3$call$formula),
                 as.character(glm3$call$formula))
    ## Same coef
    expect_equal(coef(mreg_logistic_fit3),
                 coef(glm3))
    ## Same vcov
    expect_equal(vcov(mreg_logistic_fit3),
                 vcov(glm3))

})


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
