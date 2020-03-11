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
### Helper functions
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

test_that("Sigma_sigma_hat_sq extracts the variance estimate for sigma^2", {

    data(pbc)
    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L))

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
    expect_equal(Sigma_sigma_hat_sq(lm_fit),
                 ## (sigma_hat^2)^2 / (n-p)
                 ((sigma(lm_fit))^2)^2 / lm_fit$df.residual)

})
