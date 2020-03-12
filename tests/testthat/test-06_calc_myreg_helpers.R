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

    describe("yreg linear", {
        it("extracts coef correctly when there is no interaction", {
            yreg_fit3 <- fit_yreg(yreg = "linear",
                                  data = pbc_cc,
                                  yvar = "alk.phos",
                                  avar = "trt",
                                  mvar = "bili",
                                  cvar = c("age","male","stage"),
                                  interaction = FALSE,
                                  eventvar = NULL)
            vars <- c("(Intercept)","trt","bili","age","male","stage")
            expect_equal(theta_hat(yreg = "linear",
                                   yreg_fit = yreg_fit3,
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   interaction = FALSE),
                         coef(yreg_fit3)[vars])
        })
        it("extracts coef correctly when there is an interaction", {
            yreg_fit3 <- fit_yreg(yreg = "linear",
                                  data = pbc_cc,
                                  yvar = "alk.phos",
                                  avar = "trt",
                                  mvar = "bili",
                                  cvar = c("age","male","stage"),
                                  interaction = TRUE,
                                  eventvar = NULL)
            vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
            expect_equal(theta_hat(yreg = "linear",
                                   yreg_fit = yreg_fit3,
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   interaction = TRUE),
                         coef(yreg_fit3)[vars])
        })
    })
    describe("yreg logistic", {
        it("extracts coef correctly when there is no interaction", {
            yreg_fit3 <- fit_yreg(yreg = "logistic",
                                  data = pbc_cc,
                                  yvar = "spiders",
                                  avar = "trt",
                                  mvar = "bili",
                                  cvar = c("age","male","stage"),
                                  interaction = FALSE,
                                  eventvar = NULL)
            vars <- c("(Intercept)","trt","bili","age","male","stage")
            expect_equal(theta_hat(yreg = "logistic",
                                   yreg_fit = yreg_fit3,
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   interaction = FALSE),
                         coef(yreg_fit3)[vars])
        })
        it("extracts coef correctly when there is an interaction", {
            yreg_fit3 <- fit_yreg(yreg = "logistic",
                                  data = pbc_cc,
                                  yvar = "spiders",
                                  avar = "trt",
                                  mvar = "bili",
                                  cvar = c("age","male","stage"),
                                  interaction = TRUE,
                                  eventvar = NULL)
            vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
            expect_equal(theta_hat(yreg = "logistic",
                                   yreg_fit = yreg_fit3,
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   interaction = TRUE),
                         coef(yreg_fit3)[vars])

        })
    })
    describe("yreg loglinear", {
        it("extracts coef correctly when there is no interaction", {
            yreg_fit3 <- fit_yreg(yreg = "loglinear",
                                  data = pbc_cc,
                                  yvar = "spiders",
                                  avar = "trt",
                                  mvar = "bili",
                                  cvar = c("age","male","stage"),
                                  interaction = FALSE,
                                  eventvar = NULL)
            vars <- c("(Intercept)","trt","bili","age","male","stage")
            expect_equal(theta_hat(yreg = "loglinear",
                                   yreg_fit = yreg_fit3,
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   interaction = FALSE),
                         coef(yreg_fit3)[vars])
        })
        it("extracts coef correctly when there is an interaction", {
            yreg_fit3 <- fit_yreg(yreg = "loglinear",
                                  data = pbc_cc,
                                  yvar = "spiders",
                                  avar = "trt",
                                  mvar = "bili",
                                  cvar = c("age","male","stage"),
                                  interaction = TRUE,
                                  eventvar = NULL)
            vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
            expect_equal(theta_hat(yreg = "poisson",
                                   yreg_fit = yreg_fit3,
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   interaction = TRUE),
                         coef(yreg_fit3)[vars])
        })
    })
    describe("yreg poisson", {
        ## Use platelet as a fake count variable
        it("extracts coef correctly when there is no interaction", {
            yreg_fit3 <- fit_yreg(yreg = "poisson",
                                  data = pbc_cc,
                                  yvar = "platelet",
                                  avar = "trt",
                                  mvar = "bili",
                                  cvar = c("age","male","stage"),
                                  interaction = FALSE,
                                  eventvar = NULL)
            vars <- c("(Intercept)","trt","bili","age","male","stage")
            expect_equal(theta_hat(yreg = "poisson",
                                   yreg_fit = yreg_fit3,
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   interaction = FALSE),
                         coef(yreg_fit3)[vars])
        })
        it("extracts coef correctly when there is an interaction", {
            yreg_fit3 <- fit_yreg(yreg = "poisson",
                                  data = pbc_cc,
                                  yvar = "platelet",
                                  avar = "trt",
                                  mvar = "bili",
                                  cvar = c("age","male","stage"),
                                  interaction = TRUE,
                                  eventvar = NULL)
            vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
            expect_equal(theta_hat(yreg = "poisson",
                                   yreg_fit = yreg_fit3,
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   interaction = TRUE),
                         coef(yreg_fit3)[vars])
        })
    })
    describe("yreg negbin", {
        ## Use platelet as a fake count variable
        it("extracts coef correctly when there is no interaction", {
            yreg_fit3 <- fit_yreg(yreg = "negbin",
                                  data = pbc_cc,
                                  yvar = "platelet",
                                  avar = "trt",
                                  mvar = "bili",
                                  cvar = c("age","male","stage"),
                                  interaction = FALSE,
                                  eventvar = NULL)
            vars <- c("(Intercept)","trt","bili","age","male","stage")
            expect_equal(theta_hat(yreg = "negbin",
                                   yreg_fit = yreg_fit3,
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   interaction = FALSE),
                         coef(yreg_fit3)[vars])
        })
        it("extracts coef correctly when there is an interaction", {
            yreg_fit3 <- fit_yreg(yreg = "negbin",
                                  data = pbc_cc,
                                  yvar = "platelet",
                                  avar = "trt",
                                  mvar = "bili",
                                  cvar = c("age","male","stage"),
                                  interaction = TRUE,
                                  eventvar = NULL)
            vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
            expect_equal(theta_hat(yreg = "negbin",
                                   yreg_fit = yreg_fit3,
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   interaction = TRUE),
                         coef(yreg_fit3)[vars])
        })
    })
    describe("yreg survCox", {
        it("extracts coef correctly when there is no interaction", {
            yreg_fit3 <- fit_yreg(yreg = "survCox",
                                  data = pbc_cc,
                                  yvar = "time",
                                  avar = "trt",
                                  mvar = "bili",
                                  cvar = c("age","male","stage"),
                                  interaction = FALSE,
                                  eventvar = "status")
            vars <- c("(Intercept)","trt","bili","age","male","stage")
            expect_equal(theta_hat(yreg = "survCox",
                                   yreg_fit = yreg_fit3,
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   interaction = FALSE),
                         coef(yreg_fit3)[vars])
        })
        it("extracts coef correctly when there is an interaction", {
            yreg_fit3 <- fit_yreg(yreg = "survCox",
                                  data = pbc_cc,
                                  yvar = "time",
                                  avar = "trt",
                                  mvar = "bili",
                                  cvar = c("age","male","stage"),
                                  interaction = TRUE,
                                  eventvar = "status")
            vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
            expect_equal(theta_hat(yreg = "survCox",
                                   yreg_fit = yreg_fit3,
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   interaction = TRUE),
                         coef(yreg_fit3)[vars])
        })
    })
    describe("yreg survAFT_exp", {
        it("extracts coef correctly when there is no interaction", {
            yreg_fit3 <- fit_yreg(yreg = "survAFT_exp",
                                  data = pbc_cc,
                                  yvar = "time",
                                  avar = "trt",
                                  mvar = "bili",
                                  cvar = c("age","male","stage"),
                                  interaction = FALSE,
                                  eventvar = "status")
            vars <- c("(Intercept)","trt","bili","age","male","stage")
            expect_equal(theta_hat(yreg = "survAFT_exp",
                                   yreg_fit = yreg_fit3,
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   interaction = FALSE),
                         coef(yreg_fit3)[vars])
        })
        it("extracts coef correctly when there is an interaction", {
            yreg_fit3 <- fit_yreg(yreg = "survAFT_exp",
                                  data = pbc_cc,
                                  yvar = "time",
                                  avar = "trt",
                                  mvar = "bili",
                                  cvar = c("age","male","stage"),
                                  interaction = TRUE,
                                  eventvar = "status")
            vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
            expect_equal(theta_hat(yreg = "survAFT_exp",
                                   yreg_fit = yreg_fit3,
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   interaction = TRUE),
                         coef(yreg_fit3)[vars])
        })
    })
    describe("yreg survAFT_weibull", {
        it("extracts coef correctly when there is no interaction", {
            yreg_fit3 <- fit_yreg(yreg = "survAFT_weibull",
                                  data = pbc_cc,
                                  yvar = "time",
                                  avar = "trt",
                                  mvar = "bili",
                                  cvar = c("age","male","stage"),
                                  interaction = FALSE,
                                  eventvar = "status")
            vars <- c("(Intercept)","trt","bili","age","male","stage")
            expect_equal(theta_hat(yreg = "survAFT_weibull",
                                   yreg_fit = yreg_fit3,
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   interaction = FALSE),
                         coef(yreg_fit3)[vars])
        })
        it("extracts coef correctly when there is an interaction", {
            yreg_fit3 <- fit_yreg(yreg = "survAFT_weibull",
                                  data = pbc_cc,
                                  yvar = "time",
                                  avar = "trt",
                                  mvar = "bili",
                                  cvar = c("age","male","stage"),
                                  interaction = TRUE,
                                  eventvar = "status")
            vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
            expect_equal(theta_hat(yreg = "survAFT_weibull",
                                   yreg_fit = yreg_fit3,
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   interaction = TRUE),
                         coef(yreg_fit3)[vars])
        })
    })

})


###
### Tests for vcov extractors
################################################################################

## BDD-style
## https://github.com/r-lib/testthat/issues/747
describe("Sigma_beta_hat", {

    data(pbc)
    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L))

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


describe("Sigma_sigma_sq_hat", {

    data(pbc)
    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L))

    it("extracts the variance estimate for sigma^2", {

        ## FIXME: This is tentative and suspicious.
        ## Chapter 2
        ## Quadratic Forms of Random Variables
        ## http://pages.stat.wisc.edu/~st849-1/lectures/Ch02.pdf
        ## Corollary 6. In a full-rank Gaussian model for M with p covariates X.
        ## ||M - X beta_hat||^2 / sigma^2 ~ (central chi-squared DF = (n - p))
        ## Var(central chi-squared DF = (n - p)) = (n - p)
        ## Var(sigma_hat^2) = Var(1/(n-p) * ||M - X beta_hat||^2)
        ##                  = 1/(n-p)^2 * Var(||M - X beta_hat||^2 * sigma^2/sigma^2)
        ##                  = 1/(n-p)^2 * (sigma^2)^2 * Var(||M - X beta_hat||^2 /sigma^2)
        ##                  = 1/(n-p)^2 * (sigma^2)^2 * (n-p)
        ##                  = (sigma^2)^2 / (n-p)

        ## lm
        lm_fit <- lm(formula = alk.phos ~ trt + bili,
                     data    = pbc_cc)
        ## Derivation above.
        expect_equal(Sigma_sigma_sq_hat(lm_fit),
                     ## (sigma_hat^2)^2 / (n-p)
                     matrix(((sigma(lm_fit))^2)^2 / lm_fit$df.residual))
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

    describe("yreg linear", {
        it("extracts vcov correctly when there is no interaction", {
            yreg_fit3 <- fit_yreg(yreg = "linear",
                                  data = pbc_cc,
                                  yvar = "alk.phos",
                                  avar = "trt",
                                  mvar = "bili",
                                  cvar = c("age","male","stage"),
                                  interaction = FALSE,
                                  eventvar = NULL)
            vars <- c("(Intercept)","trt","bili","age","male","stage")
            expect_equal(Sigma_theta_hat(yreg = "linear",
                                         yreg_fit = yreg_fit3,
                                         avar = "trt",
                                         mvar = "bili",
                                         cvar = c("age","male","stage"),
                                         interaction = FALSE),
                         vcov(yreg_fit3)[vars,vars])
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
        })
    })
    describe("yreg logistic", {
        it("extracts vcov correctly when there is no interaction", {
            yreg_fit3 <- fit_yreg(yreg = "logistic",
                                  data = pbc_cc,
                                  yvar = "spiders",
                                  avar = "trt",
                                  mvar = "bili",
                                  cvar = c("age","male","stage"),
                                  interaction = FALSE,
                                  eventvar = NULL)
            vars <- c("(Intercept)","trt","bili","age","male","stage")
            expect_equal(Sigma_theta_hat(yreg = "logistic",
                                         yreg_fit = yreg_fit3,
                                         avar = "trt",
                                         mvar = "bili",
                                         cvar = c("age","male","stage"),
                                         interaction = FALSE),
                         vcov(yreg_fit3)[vars,vars])
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

        })
    })
    describe("yreg loglinear", {
        it("extracts vcov correctly when there is no interaction", {
            yreg_fit3 <- fit_yreg(yreg = "loglinear",
                                  data = pbc_cc,
                                  yvar = "spiders",
                                  avar = "trt",
                                  mvar = "bili",
                                  cvar = c("age","male","stage"),
                                  interaction = FALSE,
                                  eventvar = NULL)
            vars <- c("(Intercept)","trt","bili","age","male","stage")
            expect_equal(Sigma_theta_hat(yreg = "loglinear",
                                         yreg_fit = yreg_fit3,
                                         avar = "trt",
                                         mvar = "bili",
                                         cvar = c("age","male","stage"),
                                         interaction = FALSE),
                         vcov(yreg_fit3)[vars,vars])
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
        })
    })
    describe("yreg poisson", {
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
            vars <- c("(Intercept)","trt","bili","age","male","stage")
            expect_equal(Sigma_theta_hat(yreg = "poisson",
                                         yreg_fit = yreg_fit3,
                                         avar = "trt",
                                         mvar = "bili",
                                         cvar = c("age","male","stage"),
                                         interaction = FALSE),
                         vcov(yreg_fit3)[vars,vars])
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
        })
    })
    describe("yreg negbin", {
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
            vars <- c("(Intercept)","trt","bili","age","male","stage")
            expect_equal(Sigma_theta_hat(yreg = "negbin",
                                         yreg_fit = yreg_fit3,
                                         avar = "trt",
                                         mvar = "bili",
                                         cvar = c("age","male","stage"),
                                         interaction = FALSE),
                         vcov(yreg_fit3)[vars,vars])
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
        })
    })
    describe("yreg survCox", {
        it("extracts vcov correctly when there is no interaction", {
            yreg_fit3 <- fit_yreg(yreg = "survCox",
                                  data = pbc_cc,
                                  yvar = "time",
                                  avar = "trt",
                                  mvar = "bili",
                                  cvar = c("age","male","stage"),
                                  interaction = FALSE,
                                  eventvar = "status")
            vars <- c("(Intercept)","trt","bili","age","male","stage")
            expect_equal(Sigma_theta_hat(yreg = "survCox",
                                         yreg_fit = yreg_fit3,
                                         avar = "trt",
                                         mvar = "bili",
                                         cvar = c("age","male","stage"),
                                         interaction = FALSE),
                         vcov(yreg_fit3)[vars,vars])
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
            vars <- c("(Intercept)","trt","bili","trt:bili","age","male","stage")
            expect_equal(Sigma_theta_hat(yreg = "survCox",
                                         yreg_fit = yreg_fit3,
                                         avar = "trt",
                                         mvar = "bili",
                                         cvar = c("age","male","stage"),
                                         interaction = TRUE),
                         vcov(yreg_fit3)[vars,vars])
        })
    })
    describe("yreg survAFT_exp", {
        it("extracts vcov correctly when there is no interaction", {
            yreg_fit3 <- fit_yreg(yreg = "survAFT_exp",
                                  data = pbc_cc,
                                  yvar = "time",
                                  avar = "trt",
                                  mvar = "bili",
                                  cvar = c("age","male","stage"),
                                  interaction = FALSE,
                                  eventvar = "status")
            vars <- c("(Intercept)","trt","bili","age","male","stage")
            expect_equal(Sigma_theta_hat(yreg = "survAFT_exp",
                                         yreg_fit = yreg_fit3,
                                         avar = "trt",
                                         mvar = "bili",
                                         cvar = c("age","male","stage"),
                                         interaction = FALSE),
                         vcov(yreg_fit3)[vars,vars])
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
        })
    })
    describe("yreg survAFT_weibull", {
        it("extracts vcov correctly when there is no interaction", {
            yreg_fit3 <- fit_yreg(yreg = "survAFT_weibull",
                                  data = pbc_cc,
                                  yvar = "time",
                                  avar = "trt",
                                  mvar = "bili",
                                  cvar = c("age","male","stage"),
                                  interaction = FALSE,
                                  eventvar = "status")
            vars <- c("(Intercept)","trt","bili","age","male","stage")
            expect_equal(Sigma_theta_hat(yreg = "survAFT_weibull",
                                         yreg_fit = yreg_fit3,
                                         avar = "trt",
                                         mvar = "bili",
                                         cvar = c("age","male","stage"),
                                         interaction = FALSE),
                         vcov(yreg_fit3)[vars,vars])
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
        })
    })

})
