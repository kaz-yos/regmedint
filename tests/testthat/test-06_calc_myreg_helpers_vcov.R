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
### Tests for vcov extractors
################################################################################

## BDD-style
## https://github.com/r-lib/testthat/issues/747
describe("Sigma_beta_hat", {
    ##
    data(pbc)
    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L))
    ##
    describe("Sigma_beta_hat (NULL cvar)", {
        it("extracts vcov from linear models correctly without cvar", {
            lm3 <- lm(formula = bili ~ trt,
                      data = pbc_cc)
            expect_equal(Sigma_beta_hat(mreg = "linear",
                                        mreg_fit = lm3,
                                        avar = c("trt"),
                                        cvar = NULL),
                         vcov(lm3))
            vars <- c("(Intercept)","trt")
            expect_equal(Sigma_beta_hat(mreg = "linear",
                                        mreg_fit = lm3,
                                        avar = c("trt"),
                                        cvar = NULL),
                         vcov(lm3)[vars,vars])
        })
        it("extracts vcov from logistic models correctly without cvar", {
            glm3 <- glm(formula = hepato ~ trt,
                        family = binomial(link = "logit"),
                        data = pbc_cc)
            expect_equal(Sigma_beta_hat(mreg = "logistic",
                                        mreg_fit = glm3,
                                        avar = c("trt"),
                                        cvar = NULL),
                         vcov(glm3))
            vars <- c("(Intercept)","trt")
            expect_equal(Sigma_beta_hat(mreg = "logistic",
                                        mreg_fit = glm3,
                                        avar = c("trt"),
                                        cvar = NULL),
                         vcov(glm3)[vars,vars])
        })
    })
    ##
    describe("Sigma_beta_hat (1 cvar)", {
        it("extracts vcov from linear models correctly with one cvar", {
            lm3 <- lm(formula = bili ~ trt + age,
                      data = pbc_cc)
            expect_equal(Sigma_beta_hat(mreg = "linear",
                                        mreg_fit = lm3,
                                        avar = c("trt"),
                                        cvar = c("age")),
                         vcov(lm3))
            vars <- c("(Intercept)","trt","age")
            expect_equal(Sigma_beta_hat(mreg = "linear",
                                        mreg_fit = lm3,
                                        avar = c("trt"),
                                        cvar = c("age")),
                         vcov(lm3)[vars,vars])
        })
        it("extracts vcov from logistic models correctly", {
            glm3 <- glm(formula = hepato ~ trt + age,
                        family = binomial(link = "logit"),
                        data = pbc_cc)
            expect_equal(Sigma_beta_hat(mreg = "logistic",
                                        mreg_fit = glm3,
                                        avar = c("trt"),
                                        cvar = c("age")),
                         vcov(glm3))
            vars <- c("(Intercept)","trt","age")
            expect_equal(Sigma_beta_hat(mreg = "logistic",
                                        mreg_fit = glm3,
                                        avar = c("trt"),
                                        cvar = c("age")),
                         vcov(glm3)[vars,vars])
        })
    })
    ##
    describe("Sigma_beta_hat (3 cvar)", {
        it("extracts vcov from linear models correctly", {
            lm3 <- lm(formula = bili ~ trt + age + male + stage,
                      data = pbc_cc)
            expect_equal(Sigma_beta_hat(mreg = "linear",
                                        mreg_fit = lm3,
                                        avar = c("trt"),
                                        cvar = c("age","male","stage")),
                         vcov(lm3))
            vars <- c("(Intercept)","trt","age","stage","male")
            expect_equal(Sigma_beta_hat(mreg = "linear",
                                        mreg_fit = lm3,
                                        avar = c("trt"),
                                        cvar = c("age","stage","male")),
                         vcov(lm3)[vars,vars])
        })
        it("extracts vcov from logistic models correctly", {
            glm3 <- glm(formula = hepato ~ trt + age + male + stage,
                        family = binomial(link = "logit"),
                        data = pbc_cc)
            expect_equal(Sigma_beta_hat(mreg = "logistic",
                                        mreg_fit = glm3,
                                        avar = c("trt"),
                                        cvar = c("age","male","stage")),
                         vcov(glm3))
            vars <- c("(Intercept)","trt","age","stage","male")
            expect_equal(Sigma_beta_hat(mreg = "logistic",
                                        mreg_fit = glm3,
                                        avar = c("trt"),
                                        cvar = c("age","stage","male")),
                         vcov(glm3)[vars,vars])
        })
    })
})


describe("Sigma_sigma_sq_hat", {

    data(pbc)
    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L))

    it("extracts the variance estimate for sigma^2", {

        ## VanderWeele 2015. p470 states
        ## 2 * (sigma_hat^2)^2 / (n-p)
        ## Also see the comment in the body of the function.

        ## lm
        lm_fit <- lm(formula = alk.phos ~ trt + bili,
                     data    = pbc_cc)
        expect_equal(Sigma_sigma_sq_hat(lm_fit),
                     ## 2 * (sigma_hat^2)^2 / (n-p)
                     matrix(2 * ((sigma(lm_fit))^2)^2 / lm_fit$df.residual))
        ## Must be a matrix for later manipulation
        expect_equal(dim(Sigma_sigma_sq_hat(lm_fit)),
                     c(1,1))
    })
})


describe("Sigma_theta_hat", {
    data(pbc)
    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L),
               status = if_else(status == 0, 0L, 1L))
    ##
    describe("Sigma_theta_hat (NULL cvar)", {

        describe("Sigma_theta_hat (NULL cvar) for yreg linear", {
            it("extracts vcov correctly when there is no interaction", {
                yreg_fit0 <- fit_yreg(yreg = "linear",
                                      data = pbc_cc,
                                      yvar = "alk.phos",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- NULL
                vars <- c("(Intercept)","trt","bili","trt:bili")
                ref_vcov <- Matrix::bdiag(vcov(yreg_fit0)[vars1,vars1],
                                          matrix(0))
                dimnames(ref_vcov) <- list(vars,vars)
                expect_equal(Sigma_theta_hat(yreg = "linear",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = FALSE),
                             ref_vcov[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "linear",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = FALSE) %>% dim(),
                             dim(vcov(yreg_fit0)) + 1)
            })
            it("extracts vcov correctly when there is an interaction", {
                yreg_fit0 <- fit_yreg(yreg = "linear",
                                      data = pbc_cc,
                                      yvar = "alk.phos",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili")
                expect_equal(Sigma_theta_hat(yreg = "linear",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = TRUE),
                             vcov(yreg_fit0)[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "linear",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = TRUE) %>% dim(),
                             dim(vcov(yreg_fit0)) + 0)
            })
        })
        describe("Sigma_theta_hat (NULL cvar) for yreg logistic", {
            it("extracts vcov correctly when there is no interaction", {
                yreg_fit0 <- fit_yreg(yreg = "logistic",
                                      data = pbc_cc,
                                      yvar = "spiders",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- NULL
                vars <- c("(Intercept)","trt","bili","trt:bili")
                ref_vcov <- Matrix::bdiag(vcov(yreg_fit0)[vars1,vars1],
                                          matrix(0))
                dimnames(ref_vcov) <- list(vars,vars)
                expect_equal(Sigma_theta_hat(yreg = "logistic",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = FALSE),
                             ref_vcov[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "logistic",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = FALSE) %>% dim(),
                             dim(vcov(yreg_fit0)) + 1)
            })
            it("extracts vcov correctly when there is an interaction", {
                yreg_fit0 <- fit_yreg(yreg = "logistic",
                                      data = pbc_cc,
                                      yvar = "spiders",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili")
                expect_equal(Sigma_theta_hat(yreg = "logistic",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = TRUE),
                             vcov(yreg_fit0)[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "logistic",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = TRUE) %>% dim(),
                             dim(vcov(yreg_fit0)) + 0)
            })
        })
        describe("Sigma_theta_hat (NULL cvar) for yreg loglinear", {
            it("extracts vcov correctly when there is no interaction", {
                yreg_fit0 <- fit_yreg(yreg = "loglinear",
                                      data = pbc_cc,
                                      yvar = "spiders",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- NULL
                vars <- c("(Intercept)","trt","bili","trt:bili")
                ref_vcov <- Matrix::bdiag(vcov(yreg_fit0)[vars1,vars1],
                                          matrix(0))
                dimnames(ref_vcov) <- list(vars,vars)
                expect_equal(Sigma_theta_hat(yreg = "loglinear",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = FALSE),
                             ref_vcov[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "loglinear",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = FALSE) %>% dim(),
                             dim(vcov(yreg_fit0)) + 1)
            })
            it("extracts vcov correctly when there is an interaction", {
                yreg_fit0 <- fit_yreg(yreg = "loglinear",
                                      data = pbc_cc,
                                      yvar = "spiders",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili")
                expect_equal(Sigma_theta_hat(yreg = "poisson",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = TRUE),
                             vcov(yreg_fit0)[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "loglinear",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = TRUE) %>% dim(),
                             dim(vcov(yreg_fit0)) + 0)
            })
        })
        describe("Sigma_theta_hat (NULL cvar) for yreg poisson", {
            ## Use platelet as a fake count variable
            it("extracts vcov correctly when there is no interaction", {
                yreg_fit0 <- fit_yreg(yreg = "poisson",
                                      data = pbc_cc,
                                      yvar = "platelet",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- NULL
                vars <- c("(Intercept)","trt","bili","trt:bili")
                ref_vcov <- Matrix::bdiag(vcov(yreg_fit0)[vars1,vars1],
                                          matrix(0))
                dimnames(ref_vcov) <- list(vars,vars)
                expect_equal(Sigma_theta_hat(yreg = "poisson",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = FALSE),
                             ref_vcov[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "poisson",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = FALSE) %>% dim(),
                             dim(vcov(yreg_fit0)) + 1)

            })
            it("extracts vcov correctly when there is an interaction", {
                yreg_fit0 <- fit_yreg(yreg = "poisson",
                                      data = pbc_cc,
                                      yvar = "platelet",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili")
                expect_equal(Sigma_theta_hat(yreg = "poisson",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = TRUE),
                             vcov(yreg_fit0)[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "poisson",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = TRUE) %>% dim(),
                             dim(vcov(yreg_fit0)) + 0)
            })
        })
        describe("Sigma_theta_hat (NULL cvar) for yreg negbin", {
            ## Use platelet as a fake count variable
            it("extracts vcov correctly when there is no interaction", {
                yreg_fit0 <- fit_yreg(yreg = "negbin",
                                      data = pbc_cc,
                                      yvar = "platelet",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- NULL
                vars <- c("(Intercept)","trt","bili","trt:bili")
                ref_vcov <- Matrix::bdiag(vcov(yreg_fit0)[vars1,vars1],
                                          matrix(0))
                dimnames(ref_vcov) <- list(vars,vars)
                expect_equal(Sigma_theta_hat(yreg = "negbin",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = FALSE),
                             ref_vcov[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "negbin",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = FALSE) %>% dim(),
                             dim(vcov(yreg_fit0)) + 1)
            })
            it("extracts vcov correctly when there is an interaction", {
                yreg_fit0 <- fit_yreg(yreg = "negbin",
                                      data = pbc_cc,
                                      yvar = "platelet",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili")
                expect_equal(Sigma_theta_hat(yreg = "negbin",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = TRUE),
                             vcov(yreg_fit0)[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "negbin",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = TRUE) %>% dim(),
                             dim(vcov(yreg_fit0)) + 0)
            })
        })
        describe("Sigma_theta_hat (NULL cvar) for yreg survCox", {
            it("extracts vcov correctly when there is no interaction", {
                yreg_fit0 <- fit_yreg(yreg = "survCox",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      interaction = FALSE,
                                      eventvar = "status")
                vars1 <- c("trt","bili")
                vars2 <- NULL
                vars <- c("(Intercept)","trt","bili","trt:bili")
                ref_vcov <- Matrix::bdiag(matrix(0),
                                          vcov(yreg_fit0)[vars1,vars1],
                                          matrix(0))
                dimnames(ref_vcov) <- list(vars,vars)
                expect_equal(Sigma_theta_hat(yreg = "survCox",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = FALSE),
                             ref_vcov)
                expect_equal(Sigma_theta_hat(yreg = "survCox",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = FALSE) %>% dim(),
                             dim(vcov(yreg_fit0)) + 2)
            })
            it("extracts vcov correctly when there is an interaction", {
                yreg_fit0 <- fit_yreg(yreg = "survCox",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      interaction = TRUE,
                                      eventvar = "status")
                vars1 <- c("trt","bili","trt:bili")
                vars <- c("(Intercept)","trt","bili","trt:bili")
                ref_vcov <- Matrix::bdiag(matrix(0),
                                          vcov(yreg_fit0)[vars1,vars1])
                dimnames(ref_vcov) <- list(vars,vars)
                expect_equal(Sigma_theta_hat(yreg = "survCox",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = TRUE),
                             ref_vcov)
                expect_equal(Sigma_theta_hat(yreg = "survCox",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = TRUE) %>% dim(),
                             dim(vcov(yreg_fit0)) + 1)
            })
        })
        describe("Sigma_theta_hat (NULL cvar) for yreg survAFT_exp", {
            it("extracts vcov correctly when there is no interaction", {
                yreg_fit0 <- fit_yreg(yreg = "survAFT_exp",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      interaction = FALSE,
                                      eventvar = "status")
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- NULL
                vars <- c("(Intercept)","trt","bili","trt:bili")
                ref_vcov <- Matrix::bdiag(vcov(yreg_fit0)[vars1,vars1],
                                          matrix(0))
                dimnames(ref_vcov) <- list(vars,vars)
                expect_equal(Sigma_theta_hat(yreg = "survAFT_exp",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = FALSE),
                             ref_vcov[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "survAFT_exp",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = FALSE) %>% dim(),
                             dim(vcov(yreg_fit0)) + 1)
            })
            it("extracts vcov correctly when there is an interaction", {
                yreg_fit0 <- fit_yreg(yreg = "survAFT_exp",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      interaction = TRUE,
                                      eventvar = "status")
                vars <- c("(Intercept)","trt","bili","trt:bili")
                expect_equal(Sigma_theta_hat(yreg = "survAFT_exp",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = TRUE),
                             vcov(yreg_fit0)[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "logistic",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = TRUE) %>% dim(),
                             dim(vcov(yreg_fit0)) + 0)
            })
        })
        describe("Sigma_theta_hat (NULL cvar) for yreg survAFT_weibull", {
            it("extracts vcov correctly when there is no interaction", {
                yreg_fit0 <- fit_yreg(yreg = "survAFT_weibull",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      interaction = FALSE,
                                      eventvar = "status")
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- NULL
                vars <- c("(Intercept)","trt","bili","trt:bili")
                ref_vcov <- Matrix::bdiag(vcov(yreg_fit0)[vars1,vars1],
                                          matrix(0))
                dimnames(ref_vcov) <- list(vars,vars)
                expect_equal(Sigma_theta_hat(yreg = "survAFT_weibull",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = FALSE),
                             ref_vcov[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "survAFT_weibull",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = FALSE) %>% dim(),
                             ## -1 for Log(scale)
                             dim(vcov(yreg_fit0)) -1 + 1)
            })
            it("extracts vcov correctly when there is an interaction", {
                yreg_fit0 <- fit_yreg(yreg = "survAFT_weibull",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = NULL,
                                      interaction = TRUE,
                                      eventvar = "status")
                vars <- c("(Intercept)","trt","bili","trt:bili")
                expect_equal(Sigma_theta_hat(yreg = "survAFT_weibull",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = TRUE),
                             vcov(yreg_fit0)[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "survAFT_weibull",
                                             yreg_fit = yreg_fit0,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = NULL,
                                             interaction = TRUE) %>% dim(),
                             ## -1 for Log(scale)
                             dim(vcov(yreg_fit0)) -1 + 0)
            })
        })
    })
    ##
    describe("Sigma_theta_hat (1 cvar)", {

        describe("Sigma_theta_hat (1 cvar) for yreg linear", {
            it("extracts vcov correctly when there is no interaction", {
                yreg_fit1 <- fit_yreg(yreg = "linear",
                                      data = pbc_cc,
                                      yvar = "alk.phos",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age")
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                ref_vcov <- Matrix::bdiag(vcov(yreg_fit1)[vars1,vars1],
                                          matrix(0),
                                          vcov(yreg_fit1)[vars2,vars2])
                dimnames(ref_vcov) <- list(vars,vars)
                expect_equal(Sigma_theta_hat(yreg = "linear",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = FALSE),
                             ref_vcov[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "linear",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = FALSE) %>% dim(),
                             dim(vcov(yreg_fit1)) + 1)
            })
            it("extracts vcov correctly when there is an interaction", {
                yreg_fit1 <- fit_yreg(yreg = "linear",
                                      data = pbc_cc,
                                      yvar = "alk.phos",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                expect_equal(Sigma_theta_hat(yreg = "linear",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = TRUE),
                             vcov(yreg_fit1)[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "linear",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = TRUE) %>% dim(),
                             dim(vcov(yreg_fit1)) + 0)
            })
        })
        describe("Sigma_theta_hat (1 cvar) for yreg logistic", {
            it("extracts vcov correctly when there is no interaction", {
                yreg_fit1 <- fit_yreg(yreg = "logistic",
                                      data = pbc_cc,
                                      yvar = "spiders",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age")
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                ref_vcov <- Matrix::bdiag(vcov(yreg_fit1)[vars1,vars1],
                                          matrix(0),
                                          vcov(yreg_fit1)[vars2,vars2])
                dimnames(ref_vcov) <- list(vars,vars)
                expect_equal(Sigma_theta_hat(yreg = "logistic",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = FALSE),
                             ref_vcov[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "logistic",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = FALSE) %>% dim(),
                             dim(vcov(yreg_fit1)) + 1)
            })
            it("extracts vcov correctly when there is an interaction", {
                yreg_fit1 <- fit_yreg(yreg = "logistic",
                                      data = pbc_cc,
                                      yvar = "spiders",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                expect_equal(Sigma_theta_hat(yreg = "logistic",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = TRUE),
                             vcov(yreg_fit1)[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "logistic",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = TRUE) %>% dim(),
                             dim(vcov(yreg_fit1)) + 0)
            })
        })
        describe("Sigma_theta_hat (1 cvar) for yreg loglinear", {
            it("extracts vcov correctly when there is no interaction", {
                yreg_fit1 <- fit_yreg(yreg = "loglinear",
                                      data = pbc_cc,
                                      yvar = "spiders",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age")
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                ref_vcov <- Matrix::bdiag(vcov(yreg_fit1)[vars1,vars1],
                                          matrix(0),
                                          vcov(yreg_fit1)[vars2,vars2])
                dimnames(ref_vcov) <- list(vars,vars)
                expect_equal(Sigma_theta_hat(yreg = "loglinear",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = FALSE),
                             ref_vcov[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "loglinear",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = FALSE) %>% dim(),
                             dim(vcov(yreg_fit1)) + 1)
            })
            it("extracts vcov correctly when there is an interaction", {
                yreg_fit1 <- fit_yreg(yreg = "loglinear",
                                      data = pbc_cc,
                                      yvar = "spiders",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                expect_equal(Sigma_theta_hat(yreg = "poisson",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = TRUE),
                             vcov(yreg_fit1)[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "loglinear",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = TRUE) %>% dim(),
                             dim(vcov(yreg_fit1)) + 0)
            })
        })
        describe("Sigma_theta_hat (1 cvar) for yreg poisson", {
            ## Use platelet as a fake count variable
            it("extracts vcov correctly when there is no interaction", {
                yreg_fit1 <- fit_yreg(yreg = "poisson",
                                      data = pbc_cc,
                                      yvar = "platelet",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age")
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                ref_vcov <- Matrix::bdiag(vcov(yreg_fit1)[vars1,vars1],
                                          matrix(0),
                                          vcov(yreg_fit1)[vars2,vars2])
                dimnames(ref_vcov) <- list(vars,vars)
                expect_equal(Sigma_theta_hat(yreg = "poisson",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = FALSE),
                             ref_vcov[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "poisson",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = FALSE) %>% dim(),
                             dim(vcov(yreg_fit1)) + 1)

            })
            it("extracts vcov correctly when there is an interaction", {
                yreg_fit1 <- fit_yreg(yreg = "poisson",
                                      data = pbc_cc,
                                      yvar = "platelet",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                expect_equal(Sigma_theta_hat(yreg = "poisson",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = TRUE),
                             vcov(yreg_fit1)[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "poisson",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = TRUE) %>% dim(),
                             dim(vcov(yreg_fit1)) + 0)
            })
        })
        describe("Sigma_theta_hat (1 cvar) for yreg negbin", {
            ## Use platelet as a fake count variable
            it("extracts vcov correctly when there is no interaction", {
                yreg_fit1 <- fit_yreg(yreg = "negbin",
                                      data = pbc_cc,
                                      yvar = "platelet",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age")
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                ref_vcov <- Matrix::bdiag(vcov(yreg_fit1)[vars1,vars1],
                                          matrix(0),
                                          vcov(yreg_fit1)[vars2,vars2])
                dimnames(ref_vcov) <- list(vars,vars)
                expect_equal(Sigma_theta_hat(yreg = "negbin",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = FALSE),
                             ref_vcov[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "negbin",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = FALSE) %>% dim(),
                             dim(vcov(yreg_fit1)) + 1)
            })
            it("extracts vcov correctly when there is an interaction", {
                yreg_fit1 <- fit_yreg(yreg = "negbin",
                                      data = pbc_cc,
                                      yvar = "platelet",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                expect_equal(Sigma_theta_hat(yreg = "negbin",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = TRUE),
                             vcov(yreg_fit1)[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "negbin",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = TRUE) %>% dim(),
                             dim(vcov(yreg_fit1)) + 0)
            })
        })
        describe("Sigma_theta_hat (1 cvar) for yreg survCox", {
            it("extracts vcov correctly when there is no interaction", {
                yreg_fit1 <- fit_yreg(yreg = "survCox",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      interaction = FALSE,
                                      eventvar = "status")
                vars1 <- c("trt","bili")
                vars2 <- c("age")
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                ref_vcov <- Matrix::bdiag(matrix(0),
                                          vcov(yreg_fit1)[vars1,vars1],
                                          matrix(0),
                                          vcov(yreg_fit1)[vars2,vars2])
                dimnames(ref_vcov) <- list(vars,vars)
                expect_equal(Sigma_theta_hat(yreg = "survCox",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = FALSE),
                             ref_vcov[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "survCox",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = FALSE) %>% dim(),
                             dim(vcov(yreg_fit1)) + 2)
            })
            it("extracts vcov correctly when there is an interaction", {
                yreg_fit1 <- fit_yreg(yreg = "survCox",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      interaction = TRUE,
                                      eventvar = "status")
                vars1 <- c("trt","bili","trt:bili","age")
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                ref_vcov <- Matrix::bdiag(matrix(0),
                                          vcov(yreg_fit1)[vars1,vars1])
                dimnames(ref_vcov) <- list(vars,vars)
                expect_equal(Sigma_theta_hat(yreg = "survCox",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = TRUE),
                             ref_vcov)
                expect_equal(Sigma_theta_hat(yreg = "survCox",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = TRUE) %>% dim(),
                             dim(vcov(yreg_fit1)) + 1)
            })
        })
        describe("Sigma_theta_hat (1 cvar) for yreg survAFT_exp", {
            it("extracts vcov correctly when there is no interaction", {
                yreg_fit1 <- fit_yreg(yreg = "survAFT_exp",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      interaction = FALSE,
                                      eventvar = "status")
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age")
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                ref_vcov <- Matrix::bdiag(vcov(yreg_fit1)[vars1,vars1],
                                          matrix(0),
                                          vcov(yreg_fit1)[vars2,vars2])
                dimnames(ref_vcov) <- list(vars,vars)
                expect_equal(Sigma_theta_hat(yreg = "survAFT_exp",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = FALSE),
                             ref_vcov[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "survAFT_exp",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = FALSE) %>% dim(),
                             dim(vcov(yreg_fit1)) + 1)
            })
            it("extracts vcov correctly when there is an interaction", {
                yreg_fit1 <- fit_yreg(yreg = "survAFT_exp",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      interaction = TRUE,
                                      eventvar = "status")
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                expect_equal(Sigma_theta_hat(yreg = "survAFT_exp",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = TRUE),
                             vcov(yreg_fit1)[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "logistic",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = TRUE) %>% dim(),
                             dim(vcov(yreg_fit1)) + 0)
            })
        })
        describe("Sigma_theta_hat (1 cvar) for yreg survAFT_weibull", {
            it("extracts vcov correctly when there is no interaction", {
                yreg_fit1 <- fit_yreg(yreg = "survAFT_weibull",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      interaction = FALSE,
                                      eventvar = "status")
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age")
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                ref_vcov <- Matrix::bdiag(vcov(yreg_fit1)[vars1,vars1],
                                          matrix(0),
                                          vcov(yreg_fit1)[vars2,vars2])
                dimnames(ref_vcov) <- list(vars,vars)
                expect_equal(Sigma_theta_hat(yreg = "survAFT_weibull",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = FALSE),
                             ref_vcov[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "survAFT_weibull",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = FALSE) %>% dim(),
                             ## -1 for Log(scale)
                             dim(vcov(yreg_fit1)) -1 + 1)
            })
            it("extracts vcov correctly when there is an interaction", {
                yreg_fit1 <- fit_yreg(yreg = "survAFT_weibull",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age"),
                                      interaction = TRUE,
                                      eventvar = "status")
                vars <- c("(Intercept)","trt","bili","trt:bili","age")
                expect_equal(Sigma_theta_hat(yreg = "survAFT_weibull",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = TRUE),
                             vcov(yreg_fit1)[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "survAFT_weibull",
                                             yreg_fit = yreg_fit1,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age"),
                                             interaction = TRUE) %>% dim(),
                             ## -1 for Log(scale)
                             dim(vcov(yreg_fit1)) -1 + 0)
            })
        })
    })
    ##
    describe("Sigma_theta_hat (3 cvar)", {

        describe("Sigma_theta_hat (3 cvar) for yreg linear", {
            it("extracts vcov correctly when there is no interaction", {
                yreg_fit3 <- fit_yreg(yreg = "linear",
                                      data = pbc_cc,
                                      yvar = "alk.phos",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age","male","stage")
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                ref_vcov <- Matrix::bdiag(vcov(yreg_fit3)[vars1,vars1],
                                          matrix(0),
                                          vcov(yreg_fit3)[vars2,vars2])
                dimnames(ref_vcov) <- list(vars,vars)
                expect_equal(Sigma_theta_hat(yreg = "linear",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = FALSE),
                             ref_vcov[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "linear",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = FALSE) %>% dim(),
                             dim(vcov(yreg_fit3)) + 1)
            })
            it("extracts vcov correctly when there is an interaction", {
                yreg_fit3 <- fit_yreg(yreg = "linear",
                                      data = pbc_cc,
                                      yvar = "alk.phos",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                expect_equal(Sigma_theta_hat(yreg = "linear",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = TRUE),
                             vcov(yreg_fit3)[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "linear",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = TRUE) %>% dim(),
                             dim(vcov(yreg_fit3)) + 0)
            })
        })
        describe("Sigma_theta_hat (3 cvar) for yreg logistic", {
            it("extracts vcov correctly when there is no interaction", {
                yreg_fit3 <- fit_yreg(yreg = "logistic",
                                      data = pbc_cc,
                                      yvar = "spiders",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age","male","stage")
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                ref_vcov <- Matrix::bdiag(vcov(yreg_fit3)[vars1,vars1],
                                          matrix(0),
                                          vcov(yreg_fit3)[vars2,vars2])
                dimnames(ref_vcov) <- list(vars,vars)
                expect_equal(Sigma_theta_hat(yreg = "logistic",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = FALSE),
                             ref_vcov[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "logistic",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = FALSE) %>% dim(),
                             dim(vcov(yreg_fit3)) + 1)
            })
            it("extracts vcov correctly when there is an interaction", {
                yreg_fit3 <- fit_yreg(yreg = "logistic",
                                      data = pbc_cc,
                                      yvar = "spiders",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                expect_equal(Sigma_theta_hat(yreg = "logistic",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = TRUE),
                             vcov(yreg_fit3)[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "logistic",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = TRUE) %>% dim(),
                             dim(vcov(yreg_fit3)) + 0)
            })
        })
        describe("Sigma_theta_hat (3 cvar) for yreg loglinear", {
            it("extracts vcov correctly when there is no interaction", {
                yreg_fit3 <- fit_yreg(yreg = "loglinear",
                                      data = pbc_cc,
                                      yvar = "spiders",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age","male","stage")
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                ref_vcov <- Matrix::bdiag(vcov(yreg_fit3)[vars1,vars1],
                                          matrix(0),
                                          vcov(yreg_fit3)[vars2,vars2])
                dimnames(ref_vcov) <- list(vars,vars)
                expect_equal(Sigma_theta_hat(yreg = "loglinear",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = FALSE),
                             ref_vcov[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "loglinear",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = FALSE) %>% dim(),
                             dim(vcov(yreg_fit3)) + 1)
            })
            it("extracts vcov correctly when there is an interaction", {
                yreg_fit3 <- fit_yreg(yreg = "loglinear",
                                      data = pbc_cc,
                                      yvar = "spiders",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                expect_equal(Sigma_theta_hat(yreg = "poisson",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = TRUE),
                             vcov(yreg_fit3)[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "loglinear",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = TRUE) %>% dim(),
                             dim(vcov(yreg_fit3)) + 0)
            })
        })
        describe("Sigma_theta_hat (3 cvar) for yreg poisson", {
            ## Use platelet as a fake count variable
            it("extracts vcov correctly when there is no interaction", {
                yreg_fit3 <- fit_yreg(yreg = "poisson",
                                      data = pbc_cc,
                                      yvar = "platelet",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age","male","stage")
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                ref_vcov <- Matrix::bdiag(vcov(yreg_fit3)[vars1,vars1],
                                          matrix(0),
                                          vcov(yreg_fit3)[vars2,vars2])
                dimnames(ref_vcov) <- list(vars,vars)
                expect_equal(Sigma_theta_hat(yreg = "poisson",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = FALSE),
                             ref_vcov[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "poisson",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = FALSE) %>% dim(),
                             dim(vcov(yreg_fit3)) + 1)

            })
            it("extracts vcov correctly when there is an interaction", {
                yreg_fit3 <- fit_yreg(yreg = "poisson",
                                      data = pbc_cc,
                                      yvar = "platelet",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                expect_equal(Sigma_theta_hat(yreg = "poisson",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = TRUE),
                             vcov(yreg_fit3)[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "poisson",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = TRUE) %>% dim(),
                             dim(vcov(yreg_fit3)) + 0)
            })
        })
        describe("Sigma_theta_hat (3 cvar) for yreg negbin", {
            ## Use platelet as a fake count variable
            it("extracts vcov correctly when there is no interaction", {
                yreg_fit3 <- fit_yreg(yreg = "negbin",
                                      data = pbc_cc,
                                      yvar = "platelet",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      interaction = FALSE,
                                      eventvar = NULL)
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age","male","stage")
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                ref_vcov <- Matrix::bdiag(vcov(yreg_fit3)[vars1,vars1],
                                          matrix(0),
                                          vcov(yreg_fit3)[vars2,vars2])
                dimnames(ref_vcov) <- list(vars,vars)
                expect_equal(Sigma_theta_hat(yreg = "negbin",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = FALSE),
                             ref_vcov[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "negbin",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = FALSE) %>% dim(),
                             dim(vcov(yreg_fit3)) + 1)
            })
            it("extracts vcov correctly when there is an interaction", {
                yreg_fit3 <- fit_yreg(yreg = "negbin",
                                      data = pbc_cc,
                                      yvar = "platelet",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      interaction = TRUE,
                                      eventvar = NULL)
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                expect_equal(Sigma_theta_hat(yreg = "negbin",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = TRUE),
                             vcov(yreg_fit3)[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "negbin",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = TRUE) %>% dim(),
                             dim(vcov(yreg_fit3)) + 0)
            })
        })
        describe("Sigma_theta_hat (3 cvar) for yreg survCox", {
            it("extracts vcov correctly when there is no interaction", {
                yreg_fit3 <- fit_yreg(yreg = "survCox",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      interaction = FALSE,
                                      eventvar = "status")
                vars1 <- c("trt","bili")
                vars2 <- c("age","male","stage")
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                ref_vcov <- Matrix::bdiag(matrix(0),
                                          vcov(yreg_fit3)[vars1,vars1],
                                          matrix(0),
                                          vcov(yreg_fit3)[vars2,vars2])
                dimnames(ref_vcov) <- list(vars,vars)
                expect_equal(Sigma_theta_hat(yreg = "survCox",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = FALSE),
                             ref_vcov[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "survCox",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = FALSE) %>% dim(),
                             dim(vcov(yreg_fit3)) + 2)
            })
            it("extracts vcov correctly when there is an interaction", {
                yreg_fit3 <- fit_yreg(yreg = "survCox",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      interaction = TRUE,
                                      eventvar = "status")
                vars1 <- c("trt","bili","trt:bili","age","male","stage")
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                ref_vcov <- Matrix::bdiag(matrix(0),
                                          vcov(yreg_fit3)[vars1,vars1])
                dimnames(ref_vcov) <- list(vars,vars)
                expect_equal(Sigma_theta_hat(yreg = "survCox",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = TRUE),
                             ref_vcov)
                expect_equal(Sigma_theta_hat(yreg = "survCox",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = TRUE) %>% dim(),
                             dim(vcov(yreg_fit3)) + 1)
            })
        })
        describe("Sigma_theta_hat (3 cvar) for yreg survAFT_exp", {
            it("extracts vcov correctly when there is no interaction", {
                yreg_fit3 <- fit_yreg(yreg = "survAFT_exp",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      interaction = FALSE,
                                      eventvar = "status")
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age","male","stage")
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                ref_vcov <- Matrix::bdiag(vcov(yreg_fit3)[vars1,vars1],
                                          matrix(0),
                                          vcov(yreg_fit3)[vars2,vars2])
                dimnames(ref_vcov) <- list(vars,vars)
                expect_equal(Sigma_theta_hat(yreg = "survAFT_exp",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = FALSE),
                             ref_vcov[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "survAFT_exp",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = FALSE) %>% dim(),
                             dim(vcov(yreg_fit3)) + 1)
            })
            it("extracts vcov correctly when there is an interaction", {
                yreg_fit3 <- fit_yreg(yreg = "survAFT_exp",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      interaction = TRUE,
                                      eventvar = "status")
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                expect_equal(Sigma_theta_hat(yreg = "survAFT_exp",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = TRUE),
                             vcov(yreg_fit3)[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "logistic",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = TRUE) %>% dim(),
                             dim(vcov(yreg_fit3)) + 0)
            })
        })
        describe("Sigma_theta_hat (3 cvar) for yreg survAFT_weibull", {
            it("extracts vcov correctly when there is no interaction", {
                yreg_fit3 <- fit_yreg(yreg = "survAFT_weibull",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      interaction = FALSE,
                                      eventvar = "status")
                vars1 <- c("(Intercept)","trt","bili")
                vars2 <- c("age","male","stage")
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                ref_vcov <- Matrix::bdiag(vcov(yreg_fit3)[vars1,vars1],
                                          matrix(0),
                                          vcov(yreg_fit3)[vars2,vars2])
                dimnames(ref_vcov) <- list(vars,vars)
                expect_equal(Sigma_theta_hat(yreg = "survAFT_weibull",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = FALSE),
                             ref_vcov[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "survAFT_weibull",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = FALSE) %>% dim(),
                             ## -1 for Log(scale)
                             dim(vcov(yreg_fit3)) -1 + 1)
            })
            it("extracts vcov correctly when there is an interaction", {
                yreg_fit3 <- fit_yreg(yreg = "survAFT_weibull",
                                      data = pbc_cc,
                                      yvar = "time",
                                      avar = "trt",
                                      mvar = "bili",
                                      cvar = c("age","male","stage"),
                                      interaction = TRUE,
                                      eventvar = "status")
                vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
                expect_equal(Sigma_theta_hat(yreg = "survAFT_weibull",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = TRUE),
                             vcov(yreg_fit3)[vars,vars])
                expect_equal(Sigma_theta_hat(yreg = "survAFT_weibull",
                                             yreg_fit = yreg_fit3,
                                             avar = "trt",
                                             mvar = "bili",
                                             cvar = c("age","male","stage"),
                                             interaction = TRUE) %>% dim(),
                             ## -1 for Log(scale)
                             dim(vcov(yreg_fit3)) -1 + 0)
            })
        })
    })
})
