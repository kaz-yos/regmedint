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
### Internal function for mreg model fitting
################################################################################

describe("fit_mreg", {

    data(pbc)
    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L))

    describe("fit_mreg = linear (lm)", {
        it("works ok with cvar = NULL", {
            ## No covariates
            mreg_linear_fit0 <- fit_mreg(mreg = "linear",
                                         data = pbc_cc,
                                         avar = "trt",
                                         mvar = "bili",
                                         cvar = NULL,
                                         emm_ac_mreg = NULL)
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
        })
        it("works ok with one cvar", {
            ## One covariates
            mreg_linear_fit1 <- fit_mreg(mreg = "linear",
                                         data = pbc_cc,
                                         avar = "trt",
                                         mvar = "bili",
                                         cvar = c("age"),
                                         emm_ac_mreg = NULL)
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
        })
        it("works ok with three cvar", {
            ## Three covariates
            mreg_linear_fit3 <- fit_mreg(mreg = "linear",
                                         data = pbc_cc,
                                         avar = "trt",
                                         mvar = "bili",
                                         cvar = c("age","male","stage"),
                                         emm_ac_mreg = NULL)
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
        it("works ok with three cvar, and non-null 'emm_ac_mreg'", {
            ## Three covariates
            mreg_linear_fit4 <- fit_mreg(mreg = "linear",
                                         data = pbc_cc,
                                         avar = "trt",
                                         mvar = "bili",
                                         cvar = c("age","male","stage"),
                                         emm_ac_mreg = c("male", "stage"))
            lm4 <- lm(formula = bili ~ trt + age + male + stage + trt:male + trt:stage,
                      data = pbc_cc)
            ## Same classes
            expect_equal(class(mreg_linear_fit4),
                         class(lm4))
            ## Same formula
            expect_equal(as.character(mreg_linear_fit4$call$formula),
                         as.character(lm4$call$formula))
            ## Same coef
            expect_equal(coef(mreg_linear_fit4),
                         coef(lm4))
            ## Same vcov
            expect_equal(vcov(mreg_linear_fit4),
                         vcov(lm4))
        })
    })
    ##
    describe("fit_mreg = logistic (glm)", {
        it("works ok with cvar = NULL", {
            ## No covariates
            mreg_logistic_fit0 <- fit_mreg(mreg = "logistic",
                                           data = pbc_cc,
                                           avar = "trt",
                                           mvar = "hepato",
                                           cvar = NULL,
                                           emm_ac_mreg = NULL)
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
        })
        it("works ok with one cvar", {
            ## One covariates
            mreg_logistic_fit1 <- fit_mreg(mreg = "logistic",
                                           data = pbc_cc,
                                           avar = "trt",
                                           mvar = "hepato",
                                           cvar = c("age"),
                                           emm_ac_mreg = NULL)
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
        })
        it("works ok with three cvar", {
            ## Three covariates
            mreg_logistic_fit3 <- fit_mreg(mreg = "logistic",
                                           data = pbc_cc,
                                           avar = "trt",
                                           mvar = "hepato",
                                           cvar = c("age","male","stage"),
                                           emm_ac_mreg = NULL)
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
        it("works ok with three cvar, and non-null emm_ac_mreg", {
            ## Three covariates
            mreg_logistic_fit4 <- fit_mreg(mreg = "logistic",
                                          data = pbc_cc,
                                          avar = "trt",
                                          mvar = "hepato",
                                          cvar = c("age","male","stage"),
                                          emm_ac_mreg = c("male", "stage"))
            glm4 <- glm(formula = hepato ~ trt + age + male + stage + trt:male + trt:stage,
                        family = binomial(link = "logit"),
                        data = pbc_cc)
            ## Same classes
            expect_equal(class(mreg_logistic_fit4),
                         class(glm4))
            ## Same formula
            expect_equal(as.character(mreg_logistic_fit4$call$formula),
                         as.character(glm4$call$formula))
            ## Same coef
            expect_equal(coef(mreg_logistic_fit4),
                         coef(glm4))
            ## Same vcov
            expect_equal(vcov(mreg_logistic_fit4),
                         vcov(glm4))
        })
    })
})
