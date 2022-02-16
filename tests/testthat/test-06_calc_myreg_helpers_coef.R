################################################################################
### Tests for regression-based causal mediation analysis
##
## Created on: 2020-03-11
## Author: Kazuki Yoshida
################################################################################

## Load testthat in case this is run in isolation.
library(testthat)
library(survival)
library(tidyverse)


###
### Tests for estimate extractors
################################################################################

## BDD-style
## https://github.com/r-lib/testthat/issues/747
describe("beta_hat", {

    data(pbc)
    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L))

    describe("beta_hat (zero covariates)", {
        ##
        it("extracts coef from linear models correctly", {

            lm0 <- lm(formula = bili ~ trt,
                      data = pbc_cc)
            expect_equal(beta_hat(mreg = "linear",
                                  mreg_fit = lm0,
                                  avar = c("trt"),
                                  cvar = NULL,
                                  emm_ac_mreg = NULL),
                         coef(lm0))
            vars <- c("(Intercept)","trt")
            expect_equal(beta_hat(mreg = "linear",
                                  mreg_fit = lm0,
                                  avar = c("trt"),
                                  cvar = NULL,
                                  emm_ac_mreg = NULL),
                         coef(lm0)[vars])
        })
        ##
        it("extracts coef from logistic models correctly", {

            glm0 <- glm(formula = hepato ~ trt,
                        family = binomial(link = "logit"),
                        data = pbc_cc)
            expect_equal(beta_hat(mreg = "logistic",
                                  mreg_fit = glm0,
                                  avar = c("trt"),
                                  cvar = NULL,
                                  emm_ac_mreg = NULL),
                         coef(glm0))
            vars <- c("(Intercept)","trt")
            expect_equal(beta_hat(mreg = "logistic",
                                  mreg_fit = glm0,
                                  avar = c("trt"),
                                  cvar = NULL,
                                  emm_ac_mreg = NULL),
                         coef(glm0)[vars])
        })
    })

    describe("beta_hat (1 covariate)", {
        ##
        it("extracts coef from linear models correctly", {

            lm1 <- lm(formula = bili ~ trt + age,
                      data = pbc_cc)
            expect_equal(beta_hat(mreg = "linear",
                                  mreg_fit = lm1,
                                  avar = c("trt"),
                                  cvar = c("age")),
                         coef(lm1))
            vars <- c("(Intercept)","trt","age")
            expect_equal(beta_hat(mreg = "linear",
                                  mreg_fit = lm1,
                                  avar = c("trt"),
                                  cvar = c("age")),
                         coef(lm1)[vars])
        })
        ##
        it("extracts coef from logistic models correctly", {

            glm1 <- glm(formula = hepato ~ trt + age,
                        family = binomial(link = "logit"),
                        data = pbc_cc)
            expect_equal(beta_hat(mreg = "logistic",
                                  mreg_fit = glm1,
                                  avar = c("trt"),
                                  cvar = c("age")),
                         coef(glm1))
            vars <- c("(Intercept)","trt","age")
            expect_equal(beta_hat(mreg = "logistic",
                                  mreg_fit = glm1,
                                  avar = c("trt"),
                                  cvar = c("age")),
                         coef(glm1)[vars])
        })
    })

    describe("beta_hat (3 covariates)", {
        ##
        it("extracts coef from linear models correctly", {

            lm3 <- lm(formula = bili ~ trt + age + male + stage,
                      data = pbc_cc)
            expect_equal(beta_hat(mreg = "linear",
                                  mreg_fit = lm3,
                                  avar = c("trt"),
                                  cvar = c("age","male","stage")),
                         coef(lm3))
            vars <- c("(Intercept)","trt","age","stage","male")
            expect_equal(beta_hat(mreg = "linear",
                                  mreg_fit = lm3,
                                  avar = c("trt"),
                                  cvar = c("age","stage","male")),
                         coef(lm3)[vars])
        })
        ##
        it("extracts coef from logistic models correctly", {

            glm3 <- glm(formula = hepato ~ trt + age + male + stage,
                        family = binomial(link = "logit"),
                        data = pbc_cc)
            expect_equal(beta_hat(mreg = "logistic",
                                  mreg_fit = glm3,
                                  avar = c("trt"),
                                  cvar = c("age","male","stage")),
                         coef(glm3))
            vars <- c("(Intercept)","trt","age","stage","male")
            expect_equal(beta_hat(mreg = "logistic",
                                  mreg_fit = glm3,
                                  avar = c("trt"),
                                  cvar = c("age","stage","male")),
                         coef(glm3)[vars])
        })
    })
})


describe("sigma_hat_sq", {

    data(pbc)
    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L))

    it("extracts the estimate for sigma^2", {

        ## lm
        lm_fit <- lm(formula = alk.phos ~ trt + bili,
                     data    = pbc_cc)
        expect_equal(sigma_hat_sq(lm_fit),
                     sigma(lm_fit)^2)
    })
})


describe("theta_hat", {
    data(pbc)
    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L),
               status = if_else(status == 0, 0L, 1L))
    ##
    describe("theta_hat (NULL cvar)", {
        describe("theta_hat (NULL cvar) for yreg linear", {
            it("extracts coef correctly when there is no interaction (add zero)", {
                yreg_fit0 <- fit_yreg(yreg = "linear",
                                      data = pbc_cc,
                                      yvar = "alk.phos",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- NULL
                vars <- c("(Intercept)","trt","bili","trt:bili")
                ref_coef <- c(coef(yreg_fit0)[vars1],
                              "trt:bili" = 0)
                expect_equal(theta_hat(yreg = "linear",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE),
                             ref_coef[vars])
                expect_equal(theta_hat(yreg = "linear",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE) %>% length(),
                             length(coef(yreg_fit0)) + 1)
            })
            it("extracts coef correctly when there is an interaction", {
                yreg_fit0 <- fit_yreg(yreg = "linear",
                                      data = pbc_cc,
                                      yvar = "alk.phos",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili")
                expect_equal(theta_hat(yreg = "linear",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE),
                             coef(yreg_fit0)[vars])
                expect_equal(theta_hat(yreg = "linear",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE) %>% length(),
                             length(coef(yreg_fit0)) + 0)
            })
        })
        describe("theta_hat (NULL cvar) for yreg logistic", {
            it("extracts coef correctly when there is no interaction (add zero)", {
                yreg_fit0 <- fit_yreg(yreg = "logistic",
                                      data = pbc_cc,
                                      yvar = "spiders",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- NULL
                vars <- c("(Intercept)","trt","bili","trt:bili")
                ref_coef <- c(coef(yreg_fit0)[vars1],
                              "trt:bili" = 0)
                expect_equal(theta_hat(yreg = "logistic",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE),
                             ref_coef[vars])
                expect_equal(theta_hat(yreg = "logistic",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE) %>% length(),
                             length(coef(yreg_fit0)) + 1)
            })
            it("extracts coef correctly when there is an interaction", {
                yreg_fit0 <- fit_yreg(yreg = "logistic",
                                      data = pbc_cc,
                                      yvar = "spiders",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili")
                expect_equal(theta_hat(yreg = "logistic",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE),
                             coef(yreg_fit0)[vars])
                expect_equal(theta_hat(yreg = "poisson",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE) %>% length(),
                             length(coef(yreg_fit0)) + 0)
            })
        })
        describe("theta_hat (NULL cvar) for yreg loglinear", {
            it("extracts coef correctly when there is no interaction (add zero)", {
                yreg_fit0 <- fit_yreg(yreg = "loglinear",
                                      data = pbc_cc,
                                      yvar = "spiders",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- NULL
                vars <- c("(Intercept)","trt","bili","trt:bili")
                ref_coef <- c(coef(yreg_fit0)[vars1],
                              "trt:bili" = 0)
                expect_equal(theta_hat(yreg = "loglinear",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE),
                             ref_coef[vars])
                expect_equal(theta_hat(yreg = "loglinear",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE) %>% length(),
                             length(coef(yreg_fit0)) + 1)
            })
            it("extracts coef correctly when there is an interaction", {
                yreg_fit0 <- fit_yreg(yreg = "loglinear",
                                      data = pbc_cc,
                                      yvar = "spiders",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili")
                expect_equal(theta_hat(yreg = "poisson",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE),
                             coef(yreg_fit0)[vars])
                expect_equal(theta_hat(yreg = "loglinear",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE) %>% length(),
                             length(coef(yreg_fit0)) + 0)
            })
        })
        describe("theta_hat (NULL cvar) for yreg poisson", {
            ## Use platelet as a fake count variable
            it("extracts coef correctly when there is no interaction (add zero)", {
                yreg_fit0 <- fit_yreg(yreg = "poisson",
                                      data = pbc_cc,
                                      yvar = "platelet",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- NULL
                vars <- c("(Intercept)","trt","bili","trt:bili")
                ref_coef <- c(coef(yreg_fit0)[vars1],
                              "trt:bili" = 0)
                expect_equal(theta_hat(yreg = "poisson",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE),
                             ref_coef[vars])
                expect_equal(theta_hat(yreg = "poisson",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE) %>% length(),
                             length(coef(yreg_fit0)) + 1)
            })
            it("extracts coef correctly when there is an interaction", {
                yreg_fit0 <- fit_yreg(yreg = "poisson",
                                      data = pbc_cc,
                                      yvar = "platelet",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili")
                expect_equal(theta_hat(yreg = "poisson",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE),
                             coef(yreg_fit0)[vars])
                expect_equal(theta_hat(yreg = "poisson",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE) %>% length(),
                             length(coef(yreg_fit0)) + 0)
            })
        })
        describe("theta_hat (NULL cvar) for yreg negbin", {
            ## Use platelet as a fake count variable
            it("extracts coef correctly when there is no interaction (add zero)", {
                yreg_fit0 <- fit_yreg(yreg = "negbin",
                                      data = pbc_cc,
                                      yvar = "platelet",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- NULL
                vars <- c("(Intercept)","trt","bili","trt:bili")
                ref_coef <- c(coef(yreg_fit0)[vars1],
                              "trt:bili" = 0)
                expect_equal(theta_hat(yreg = "negbin",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE),
                             ref_coef[vars])
                expect_equal(theta_hat(yreg = "negbin",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE) %>% length(),
                             length(coef(yreg_fit0)) + 1)
            })
            it("extracts coef correctly when there is an interaction", {
                yreg_fit0 <- fit_yreg(yreg = "negbin",
                                      data = pbc_cc,
                                      yvar = "platelet",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili")
                expect_equal(theta_hat(yreg = "negbin",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE),
                             coef(yreg_fit0)[vars])
                expect_equal(theta_hat(yreg = "negbin",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE) %>% length(),
                             length(coef(yreg_fit0)) + 0)
            })
        })
        describe("theta_hat (NULL cvar) for yreg survCox", {
            it("extracts coef correctly when there is no interaction (add two zeros)", {
                yreg_fit0 <- fit_yreg(yreg = "survCox",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = FALSE,
                                      eventvar = "status")
                vars1 <- c("trt","bili")
                vars2 <- NULL
                vars <- c("(Intercept)","trt","bili","trt:bili")
                ref_coef <- c("(Intercept)" = 0,
                              coef(yreg_fit0)[vars1],
                              "trt:bili" = 0)
                expect_equal(theta_hat(yreg = "survCox",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE),
                             ref_coef[vars])
                expect_equal(theta_hat(yreg = "survCox",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE) %>% length(),
                             length(coef(yreg_fit0)) + 2)
            })
            it("extracts coef correctly when there is an interaction (add zero for Intercept)", {
                yreg_fit0 <- fit_yreg(yreg = "survCox",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = TRUE,
                                      eventvar = "status")
                vars <- c("trt","bili","trt:bili")
                ref_coef <- c("(Intercept)" = 0,
                              coef(yreg_fit0)[vars])
                expect_equal(theta_hat(yreg = "survCox",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE),
                             ref_coef)
                expect_equal(theta_hat(yreg = "survCox",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE) %>% length(),
                             length(coef(yreg_fit0)) + 1)
            })
        })
        describe("theta_hat (NULL cvar) for yreg survAFT_exp", {
            it("extracts coef correctly when there is no interaction (add zero)", {
                yreg_fit0 <- fit_yreg(yreg = "survAFT_exp",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = FALSE,
                                      eventvar = "status")
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- NULL
                vars <- c("(Intercept)","trt","bili","trt:bili")
                ref_coef <- c(coef(yreg_fit0)[vars1],
                              "trt:bili" = 0)
                expect_equal(theta_hat(yreg = "survAFT_exp",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE),
                             ref_coef[vars])
                expect_equal(theta_hat(yreg = "survAFT_exp",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE) %>% length(),
                             length(coef(yreg_fit0)) + 1)
            })
            it("extracts coef correctly when there is an interaction", {
                yreg_fit0 <- fit_yreg(yreg = "survAFT_exp",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = TRUE,
                                      eventvar = "status")
                vars <- c("(Intercept)","trt","bili","trt:bili")
                expect_equal(theta_hat(yreg = "survAFT_exp",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE),
                             coef(yreg_fit0)[vars])
                expect_equal(theta_hat(yreg = "survAFT_exp",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE) %>% length(),
                             length(coef(yreg_fit0)) + 0)
            })
        })
        describe("theta_hat (NULL cvar) for yreg survAFT_weibull", {
            it("extracts coef correctly when there is no interaction (add zero)", {
                yreg_fit0 <- fit_yreg(yreg = "survAFT_weibull",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = FALSE,
                                      eventvar = "status")
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- NULL
                vars <- c("(Intercept)","trt","bili","trt:bili")
                ref_coef <- c(coef(yreg_fit0)[vars1],
                              "trt:bili" = 0)
                expect_equal(theta_hat(yreg = "survAFT_weibull",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE),
                             ref_coef[vars])
                expect_equal(theta_hat(yreg = "survAFT_weibull",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE) %>% length(),
                             length(coef(yreg_fit0)) + 1)
            })
            it("extracts coef correctly when there is an interaction", {
                yreg_fit0 <- fit_yreg(yreg = "survAFT_weibull",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = TRUE,
                                      eventvar = "status")
                vars <- c("(Intercept)","trt","bili","trt:bili")
                expect_equal(theta_hat(yreg = "survAFT_weibull",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE),
                             coef(yreg_fit0)[vars])
                expect_equal(theta_hat(yreg = "survAFT_weibull",
                                       yreg_fit = yreg_fit0,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE) %>% length(),
                             length(coef(yreg_fit0)) + 0)
            })
        })
    })
    ##
    describe("theta_hat (1 cvar)", {
        describe("theta_hat (1 cvar) for yreg linear", {
            it("extracts coef correctly when there is no interaction (add zero)", {
                yreg_fit1 <- fit_yreg(yreg = "linear",
                                      data = pbc_cc,
                                      yvar = "alk.phos",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age")
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                ref_coef <- c(coef(yreg_fit1)[vars1],
                              "trt:bili" = 0,
                              coef(yreg_fit1)[vars2])
                expect_equal(theta_hat(yreg = "linear",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE),
                             ref_coef[vars])
                expect_equal(theta_hat(yreg = "linear",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE) %>% length(),
                             length(coef(yreg_fit1)) + 1)
            })
            it("extracts coef correctly when there is an interaction", {
                yreg_fit1 <- fit_yreg(yreg = "linear",
                                      data = pbc_cc,
                                      yvar = "alk.phos",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                expect_equal(theta_hat(yreg = "linear",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE),
                             coef(yreg_fit1)[vars])
                expect_equal(theta_hat(yreg = "linear",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE) %>% length(),
                             length(coef(yreg_fit1)) + 0)
            })
        })
        describe("theta_hat (1 cvar) for yreg logistic", {
            it("extracts coef correctly when there is no interaction (add zero)", {
                yreg_fit1 <- fit_yreg(yreg = "logistic",
                                      data = pbc_cc,
                                      yvar = "spiders",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age")
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                ref_coef <- c(coef(yreg_fit1)[vars1],
                              "trt:bili" = 0,
                              coef(yreg_fit1)[vars2])
                expect_equal(theta_hat(yreg = "logistic",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE),
                             ref_coef[vars])
                expect_equal(theta_hat(yreg = "logistic",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE) %>% length(),
                             length(coef(yreg_fit1)) + 1)
            })
            it("extracts coef correctly when there is an interaction", {
                yreg_fit1 <- fit_yreg(yreg = "logistic",
                                      data = pbc_cc,
                                      yvar = "spiders",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                expect_equal(theta_hat(yreg = "logistic",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE),
                             coef(yreg_fit1)[vars])
                expect_equal(theta_hat(yreg = "poisson",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE) %>% length(),
                             length(coef(yreg_fit1)) + 0)
            })
        })
        describe("theta_hat (1 cvar) for yreg loglinear", {
            it("extracts coef correctly when there is no interaction (add zero)", {
                yreg_fit1 <- fit_yreg(yreg = "loglinear",
                                      data = pbc_cc,
                                      yvar = "spiders",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age")
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                ref_coef <- c(coef(yreg_fit1)[vars1],
                              "trt:bili" = 0,
                              coef(yreg_fit1)[vars2])
                expect_equal(theta_hat(yreg = "loglinear",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE),
                             ref_coef[vars])
                expect_equal(theta_hat(yreg = "loglinear",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE) %>% length(),
                             length(coef(yreg_fit1)) + 1)
            })
            it("extracts coef correctly when there is an interaction", {
                yreg_fit1 <- fit_yreg(yreg = "loglinear",
                                      data = pbc_cc,
                                      yvar = "spiders",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                expect_equal(theta_hat(yreg = "poisson",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE),
                             coef(yreg_fit1)[vars])
                expect_equal(theta_hat(yreg = "loglinear",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE) %>% length(),
                             length(coef(yreg_fit1)) + 0)
            })
        })
        describe("theta_hat (1 cvar) for yreg poisson", {
            ## Use platelet as a fake count variable
            it("extracts coef correctly when there is no interaction (add zero)", {
                yreg_fit1 <- fit_yreg(yreg = "poisson",
                                      data = pbc_cc,
                                      yvar = "platelet",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age")
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                ref_coef <- c(coef(yreg_fit1)[vars1],
                              "trt:bili" = 0,
                              coef(yreg_fit1)[vars2])
                expect_equal(theta_hat(yreg = "poisson",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE),
                             ref_coef[vars])
                expect_equal(theta_hat(yreg = "poisson",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE) %>% length(),
                             length(coef(yreg_fit1)) + 1)
            })
            it("extracts coef correctly when there is an interaction", {
                yreg_fit1 <- fit_yreg(yreg = "poisson",
                                      data = pbc_cc,
                                      yvar = "platelet",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                expect_equal(theta_hat(yreg = "poisson",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE),
                             coef(yreg_fit1)[vars])
                expect_equal(theta_hat(yreg = "poisson",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE) %>% length(),
                             length(coef(yreg_fit1)) + 0)
            })
        })
        describe("theta_hat (1 cvar) for yreg negbin", {
            ## Use platelet as a fake count variable
            it("extracts coef correctly when there is no interaction (add zero)", {
                yreg_fit1 <- fit_yreg(yreg = "negbin",
                                      data = pbc_cc,
                                      yvar = "platelet",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age")
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                ref_coef <- c(coef(yreg_fit1)[vars1],
                              "trt:bili" = 0,
                              coef(yreg_fit1)[vars2])
                expect_equal(theta_hat(yreg = "negbin",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE),
                             ref_coef[vars])
                expect_equal(theta_hat(yreg = "negbin",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE) %>% length(),
                             length(coef(yreg_fit1)) + 1)
            })
            it("extracts coef correctly when there is an interaction", {
                yreg_fit1 <- fit_yreg(yreg = "negbin",
                                      data = pbc_cc,
                                      yvar = "platelet",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                expect_equal(theta_hat(yreg = "negbin",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE),
                             coef(yreg_fit1)[vars])
                expect_equal(theta_hat(yreg = "negbin",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE) %>% length(),
                             length(coef(yreg_fit1)) + 0)
            })
        })
        describe("theta_hat (1 cvar) for yreg survCox", {
            it("extracts coef correctly when there is no interaction (add two zeros)", {
                yreg_fit1 <- fit_yreg(yreg = "survCox",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = FALSE,
                                      eventvar = "status")
                vars1 <- c("trt","bili")
                vars2 <- c("age")
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                ref_coef <- c("(Intercept)" = 0,
                              coef(yreg_fit1)[vars1],
                              "trt:bili" = 0,
                              coef(yreg_fit1)[vars2])
                expect_equal(theta_hat(yreg = "survCox",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE),
                             ref_coef[vars])
                expect_equal(theta_hat(yreg = "survCox",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE) %>% length(),
                             length(coef(yreg_fit1)) + 2)
            })
            it("extracts coef correctly when there is an interaction (add zero for Intercept)", {
                yreg_fit1 <- fit_yreg(yreg = "survCox",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = TRUE,
                                      eventvar = "status")
                vars <- c("trt","bili","trt:bili","age")
                ref_coef <- c("(Intercept)" = 0,
                              coef(yreg_fit1)[vars])
                expect_equal(theta_hat(yreg = "survCox",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE),
                             ref_coef)
                expect_equal(theta_hat(yreg = "survCox",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE) %>% length(),
                             length(coef(yreg_fit1)) + 1)
            })
        })
        describe("theta_hat (1 cvar) for yreg survAFT_exp", {
            it("extracts coef correctly when there is no interaction (add zero)", {
                yreg_fit1 <- fit_yreg(yreg = "survAFT_exp",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = FALSE,
                                      eventvar = "status")
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age")
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                ref_coef <- c(coef(yreg_fit1)[vars1],
                              "trt:bili" = 0,
                              coef(yreg_fit1)[vars2])
                expect_equal(theta_hat(yreg = "survAFT_exp",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE),
                             ref_coef[vars])
                expect_equal(theta_hat(yreg = "survAFT_exp",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE) %>% length(),
                             length(coef(yreg_fit1)) + 1)
            })
            it("extracts coef correctly when there is an interaction", {
                yreg_fit1 <- fit_yreg(yreg = "survAFT_exp",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = TRUE,
                                      eventvar = "status")
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                expect_equal(theta_hat(yreg = "survAFT_exp",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE),
                             coef(yreg_fit1)[vars])
                expect_equal(theta_hat(yreg = "survAFT_exp",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE) %>% length(),
                             length(coef(yreg_fit1)) + 0)
            })
        })
        describe("theta_hat (1 cvar) for yreg survAFT_weibull", {
            it("extracts coef correctly when there is no interaction (add zero)", {
                yreg_fit1 <- fit_yreg(yreg = "survAFT_weibull",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = FALSE,
                                      eventvar = "status")
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age")
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                ref_coef <- c(coef(yreg_fit1)[vars1],
                              "trt:bili" = 0,
                              coef(yreg_fit1)[vars2])
                expect_equal(theta_hat(yreg = "survAFT_weibull",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE),
                             ref_coef[vars])
                expect_equal(theta_hat(yreg = "survAFT_weibull",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE) %>% length(),
                             length(coef(yreg_fit1)) + 1)
            })
            it("extracts coef correctly when there is an interaction", {
                yreg_fit1 <- fit_yreg(yreg = "survAFT_weibull",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = TRUE,
                                      eventvar = "status")
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                expect_equal(theta_hat(yreg = "survAFT_weibull",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE),
                             coef(yreg_fit1)[vars])
                expect_equal(theta_hat(yreg = "survAFT_weibull",
                                       yreg_fit = yreg_fit1,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE) %>% length(),
                             length(coef(yreg_fit1)) + 0)
            })
        })
    })
    ##
    describe("theta_hat (3 cvar)", {
        describe("theta_hat (3 cvar) for yreg linear", {
            it("extracts coef correctly when there is no interaction (add zero)", {
                yreg_fit3 <- fit_yreg(yreg = "linear",
                                      data = pbc_cc,
                                      yvar = "alk.phos",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age","male","stage")
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                ref_coef <- c(coef(yreg_fit3)[vars1],
                              "trt:bili" = 0,
                              coef(yreg_fit3)[vars2])
                expect_equal(theta_hat(yreg = "linear",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE),
                             ref_coef[vars])
                expect_equal(theta_hat(yreg = "linear",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE) %>% length(),
                             length(coef(yreg_fit3)) + 1)
            })
            it("extracts coef correctly when there is an interaction", {
                yreg_fit3 <- fit_yreg(yreg = "linear",
                                      data = pbc_cc,
                                      yvar = "alk.phos",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                expect_equal(theta_hat(yreg = "linear",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE),
                             coef(yreg_fit3)[vars])
                expect_equal(theta_hat(yreg = "linear",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE) %>% length(),
                             length(coef(yreg_fit3)) + 0)
            })
        })
        describe("theta_hat (3 cvar) for yreg logistic", {
            it("extracts coef correctly when there is no interaction (add zero)", {
                yreg_fit3 <- fit_yreg(yreg = "logistic",
                                      data = pbc_cc,
                                      yvar = "spiders",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age","male","stage")
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                ref_coef <- c(coef(yreg_fit3)[vars1],
                              "trt:bili" = 0,
                              coef(yreg_fit3)[vars2])
                expect_equal(theta_hat(yreg = "logistic",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE),
                             ref_coef[vars])
                expect_equal(theta_hat(yreg = "logistic",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE) %>% length(),
                             length(coef(yreg_fit3)) + 1)
            })
            it("extracts coef correctly when there is an interaction", {
                yreg_fit3 <- fit_yreg(yreg = "logistic",
                                      data = pbc_cc,
                                      yvar = "spiders",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                expect_equal(theta_hat(yreg = "logistic",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE),
                             coef(yreg_fit3)[vars])
                expect_equal(theta_hat(yreg = "poisson",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE) %>% length(),
                             length(coef(yreg_fit3)) + 0)
            })
        })
        describe("theta_hat (3 cvar) for yreg loglinear", {
            it("extracts coef correctly when there is no interaction (add zero)", {
                yreg_fit3 <- fit_yreg(yreg = "loglinear",
                                      data = pbc_cc,
                                      yvar = "spiders",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age","male","stage")
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                ref_coef <- c(coef(yreg_fit3)[vars1],
                              "trt:bili" = 0,
                              coef(yreg_fit3)[vars2])
                expect_equal(theta_hat(yreg = "loglinear",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE),
                             ref_coef[vars])
                expect_equal(theta_hat(yreg = "loglinear",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE) %>% length(),
                             length(coef(yreg_fit3)) + 1)
            })
            it("extracts coef correctly when there is an interaction", {
                yreg_fit3 <- fit_yreg(yreg = "loglinear",
                                      data = pbc_cc,
                                      yvar = "spiders",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                expect_equal(theta_hat(yreg = "poisson",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE),
                             coef(yreg_fit3)[vars])
                expect_equal(theta_hat(yreg = "loglinear",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE) %>% length(),
                             length(coef(yreg_fit3)) + 0)
            })
        })
        describe("theta_hat (3 cvar) for yreg poisson", {
            ## Use platelet as a fake count variable
            it("extracts coef correctly when there is no interaction (add zero)", {
                yreg_fit3 <- fit_yreg(yreg = "poisson",
                                      data = pbc_cc,
                                      yvar = "platelet",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age","male","stage")
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                ref_coef <- c(coef(yreg_fit3)[vars1],
                              "trt:bili" = 0,
                              coef(yreg_fit3)[vars2])
                expect_equal(theta_hat(yreg = "poisson",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE),
                             ref_coef[vars])
                expect_equal(theta_hat(yreg = "poisson",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE) %>% length(),
                             length(coef(yreg_fit3)) + 1)
            })
            it("extracts coef correctly when there is an interaction", {
                yreg_fit3 <- fit_yreg(yreg = "poisson",
                                      data = pbc_cc,
                                      yvar = "platelet",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                expect_equal(theta_hat(yreg = "poisson",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE),
                             coef(yreg_fit3)[vars])
                expect_equal(theta_hat(yreg = "poisson",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE) %>% length(),
                             length(coef(yreg_fit3)) + 0)
            })
        })
        describe("theta_hat (3 cvar) for yreg negbin", {
            ## Use platelet as a fake count variable
            it("extracts coef correctly when there is no interaction (add zero)", {
                yreg_fit3 <- fit_yreg(yreg = "negbin",
                                      data = pbc_cc,
                                      yvar = "platelet",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age","male","stage")
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                ref_coef <- c(coef(yreg_fit3)[vars1],
                              "trt:bili" = 0,
                              coef(yreg_fit3)[vars2])
                expect_equal(theta_hat(yreg = "negbin",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE),
                             ref_coef[vars])
                expect_equal(theta_hat(yreg = "negbin",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE) %>% length(),
                             length(coef(yreg_fit3)) + 1)
            })
            it("extracts coef correctly when there is an interaction", {
                yreg_fit3 <- fit_yreg(yreg = "negbin",
                                      data = pbc_cc,
                                      yvar = "platelet",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                expect_equal(theta_hat(yreg = "negbin",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE),
                             coef(yreg_fit3)[vars])
                expect_equal(theta_hat(yreg = "negbin",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE) %>% length(),
                             length(coef(yreg_fit3)) + 0)
            })
        })
        describe("theta_hat (3 cvar) for yreg survCox", {
            it("extracts coef correctly when there is no interaction (add two zeros)", {
                yreg_fit3 <- fit_yreg(yreg = "survCox",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = FALSE,
                                      eventvar = "status")
                vars1 <- c("trt","bili")
                vars2 <- c("age","male","stage")
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                ref_coef <- c("(Intercept)" = 0,
                              coef(yreg_fit3)[vars1],
                              "trt:bili" = 0,
                              coef(yreg_fit3)[vars2])
                expect_equal(theta_hat(yreg = "survCox",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE),
                             ref_coef[vars])
                expect_equal(theta_hat(yreg = "survCox",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE) %>% length(),
                             length(coef(yreg_fit3)) + 2)
            })
            it("extracts coef correctly when there is an interaction (Add zero for Intercept)", {
                yreg_fit3 <- fit_yreg(yreg = "survCox",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = TRUE,
                                      eventvar = "status")
                vars <- c("trt","bili","trt:bili","age","male","stage")
                ref_coef <- c("(Intercept)" = 0,
                              coef(yreg_fit3)[vars])
                expect_equal(theta_hat(yreg = "survCox",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE),
                             ref_coef)
                expect_equal(theta_hat(yreg = "survCox",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE) %>% length(),
                             length(coef(yreg_fit3)) + 1)
            })
        })
        describe("theta_hat (3 cvar) for yreg survAFT_exp", {
            it("extracts coef correctly when there is no interaction (add zero)", {
                yreg_fit3 <- fit_yreg(yreg = "survAFT_exp",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = FALSE,
                                      eventvar = "status")
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age","male","stage")
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                ref_coef <- c(coef(yreg_fit3)[vars1],
                              "trt:bili" = 0,
                              coef(yreg_fit3)[vars2])
                expect_equal(theta_hat(yreg = "survAFT_exp",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE),
                             ref_coef[vars])
                expect_equal(theta_hat(yreg = "survAFT_exp",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE) %>% length(),
                             length(coef(yreg_fit3)) + 1)
            })
            it("extracts coef correctly when there is an interaction", {
                yreg_fit3 <- fit_yreg(yreg = "survAFT_exp",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = TRUE,
                                      eventvar = "status")
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                expect_equal(theta_hat(yreg = "survAFT_exp",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE),
                             coef(yreg_fit3)[vars])
                expect_equal(theta_hat(yreg = "survAFT_exp",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE) %>% length(),
                             length(coef(yreg_fit3)) + 0)
            })
        })
        describe("theta_hat (3 cvar) for yreg survAFT_weibull", {
            it("extracts coef correctly when there is no interaction (add zero)", {
                yreg_fit3 <- fit_yreg(yreg = "survAFT_weibull",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = FALSE,
                                      eventvar = "status")
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age","male","stage")
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                ref_coef <- c(coef(yreg_fit3)[vars1],
                              "trt:bili" = 0,
                              coef(yreg_fit3)[vars2])
                expect_equal(theta_hat(yreg = "survAFT_weibull",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE),
                             ref_coef[vars])
                expect_equal(theta_hat(yreg = "survAFT_weibull",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = FALSE) %>% length(),
                             length(coef(yreg_fit3)) + 1)
            })
            it("extracts coef correctly when there is an interaction", {
                yreg_fit3 <- fit_yreg(yreg = "survAFT_weibull",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      emm_ac_yreg = NULL,
                                      emm_mc_yreg = NULL,
                                      interaction = TRUE,
                                      eventvar = "status")
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                expect_equal(theta_hat(yreg = "survAFT_weibull",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE),
                             coef(yreg_fit3)[vars])
                expect_equal(theta_hat(yreg = "survAFT_weibull",
                                       yreg_fit = yreg_fit3,
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = c("age","male","stage"),
                                       emm_ac_yreg = NULL,
                                       emm_mc_yreg = NULL,
                                       interaction = TRUE) %>% length(),
                             length(coef(yreg_fit3)) + 0)
            })
        })
    })
})
