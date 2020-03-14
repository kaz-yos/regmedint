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
###
################################################################################

describe("calc_myreg_mreg_linear_yreg_logistic", {

    data(pbc)
    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L),
               ## Combine transplant and death for testing purpose
               status = if_else(status == 0, 0L, 1L))

    describe("calc_myreg_mreg_linear_yreg_logistic(cvar = NULL)", {
        mreg_fit <- fit_mreg(mreg = "linear",
                             data = pbc_cc,
                             avar = "trt",
                             mvar = "bili",
                             cvar = NULL)
        yreg_fit <- fit_yreg(yreg = "logistic",
                             data = pbc_cc,
                             yvar = "spiders",
                             avar = "trt",
                             mvar = "bili",
                             cvar = NULL,
                             interaction = FALSE,
                             eventvar = NULL)
        myreg_funs <-
            calc_myreg_mreg_linear_yreg_logistic(mreg = "linear",
                                                 mreg_fit = mreg_fit,
                                                 yreg = "logistic",
                                                 yreg_fit = yreg_fit,
                                                 avar = "trt",
                                                 mvar = "bili",
                                                 cvar = NULL,
                                                 interaction = FALSE)
        ##
        it("returns a list of two functions", {
            expect_equal(class(myreg_funs),
                         "list")
            expect_equal(length(myreg_funs),
                         2)
        })
        it("returns functions that take 4 arguments", {
            expect_equal(formals(myreg_funs[[1]]),
                         c("a0","a1","m_cde","c_cond"))
            expect_equal(formals(myreg_funs[[2]]),
                         c("a0","a1","m_cde","c_cond"))
        })
        it("returns functions that return named vector of effect estimates", {
            expect(names(myreg_funs[[1]](1,2,3,NULL)),
                   c("cde","pnde","tnie","tnde","pnie","te","pm"))
            expect(names(myreg_funs[[2]](1,2,3,NULL)),
                   c("se_cde","se_pnde","se_tnie","se_tnde","se_pnie","se_te","se_pm"))
        })
        it("returns functions that error on inconsistent c_cond", {
            expect_error(myreg_funs[[1]](1,2,3,4))
            expect_error(myreg_funs[[2]](1,2,3,4))
        })
    })
})


###
### Effect estimation function constructor
################################################################################

## FIXME: Most of these functionalities should be factored out to be shared
## among calc_myreg_mreg_*_yreg_*_est.
describe("calc_myreg_mreg_linear_yreg_logistic_est", {
    describe("calc_myreg_mreg_linear_yreg_logistic_est (error handling)", {
        it("errors given inconsistent beta2 and theta4", {
            expect_error(
                calc_myreg_mreg_linear_yreg_logistic_est(beta0 = 1,
                                                         beta1 = 2,
                                                         beta2 = 3,
                                                         theta1 = 4,
                                                         theta2 = 5,
                                                         theta3 = 6,
                                                         theta4 = NULL,
                                                         sigma_sq = 8))
            expect_error(
                calc_myreg_mreg_linear_yreg_logistic_est(beta0 = 1,
                                                         beta1 = 2,
                                                         beta2 = 1:2,
                                                         theta1 = 4,
                                                         theta2 = 5,
                                                         theta3 = 6,
                                                         theta4 = NULL,
                                                         sigma_sq = 8))
            expect_error(
                calc_myreg_mreg_linear_yreg_logistic_est(beta0 = 1,
                                                         beta1 = 2,
                                                         beta2 = NULL,
                                                         theta1 = 4,
                                                         theta2 = 5,
                                                         theta3 = 6,
                                                         theta4 = 7,
                                                         sigma_sq = 8))
            expect_error(
                calc_myreg_mreg_linear_yreg_logistic_est(beta0 = 1,
                                                         beta1 = 2,
                                                         beta2 = NULL,
                                                         theta1 = 4,
                                                         theta2 = 5,
                                                         theta3 = 6,
                                                         theta4 = 7:8,
                                                         sigma_sq = 8))
        })
        it("errors given vector inputs in arguments other than beta2 and theta4", {
            expect_error(
                calc_myreg_mreg_linear_yreg_logistic_est(beta0 = 1:2,
                                                         beta1 = 2,
                                                         beta2 = 3,
                                                         theta1 = 4,
                                                         theta2 = 5,
                                                         theta3 = 6,
                                                         theta4 = 7,
                                                         sigma_sq = 8))
            expect_error(
                calc_myreg_mreg_linear_yreg_logistic_est(beta0 = 1,
                                                         beta1 = 2:3,
                                                         beta2 = 3,
                                                         theta1 = 4,
                                                         theta2 = 5,
                                                         theta3 = 6,
                                                         theta4 = 7,
                                                         sigma_sq = 8))
            expect_error(
                calc_myreg_mreg_linear_yreg_logistic_est(beta0 = 1,
                                                         beta1 = 2,
                                                         beta2 = 3,
                                                         theta1 = 4:5,
                                                         theta2 = 5,
                                                         theta3 = 6,
                                                         theta4 = 7,
                                                         sigma_sq = 8))
            expect_error(
                calc_myreg_mreg_linear_yreg_logistic_est(beta0 = 1,
                                                         beta1 = 2,
                                                         beta2 = 3,
                                                         theta1 = 4,
                                                         theta2 = 5:6,
                                                         theta3 = 6,
                                                         theta4 = 7,
                                                         sigma_sq = 8))
            expect_error(
                calc_myreg_mreg_linear_yreg_logistic_est(beta0 = 1,
                                                         beta1 = 2,
                                                         beta2 = 3,
                                                         theta1 = 4,
                                                         theta2 = 5,
                                                         theta3 = 6:7,
                                                         theta4 = 7,
                                                         sigma_sq = 8))
            expect_error(
                calc_myreg_mreg_linear_yreg_logistic_est(beta0 = 1,
                                                         beta1 = 2,
                                                         beta2 = 3,
                                                         theta1 = 4,
                                                         theta2 = 5,
                                                         theta3 = 6:7,
                                                         theta4 = 7,
                                                         sigma_sq = 8:9))
        })
        it("errors given NULL inputs in arguments other than beta2 and theta4", {
            expect_error(
                calc_myreg_mreg_linear_yreg_logistic_est(beta0 = NULL,
                                                         beta1 = 2,
                                                         beta2 = 3,
                                                         theta1 = 4,
                                                         theta2 = 5,
                                                         theta3 = 6,
                                                         theta4 = 7,
                                                         sigma_sq = 8))
            expect_error(
                calc_myreg_mreg_linear_yreg_logistic_est(beta0 = 1,
                                                         beta1 = NULL,
                                                         beta2 = 3,
                                                         theta1 = 4,
                                                         theta2 = 5,
                                                         theta3 = 6,
                                                         theta4 = 7,
                                                         sigma_sq = 8))
            expect_error(
                calc_myreg_mreg_linear_yreg_logistic_est(beta0 = 1,
                                                         beta1 = 2,
                                                         beta2 = 3,
                                                         theta1 = NULL,
                                                         theta2 = 5,
                                                         theta3 = 6,
                                                         theta4 = 7,
                                                         sigma_sq = 8))
            expect_error(
                calc_myreg_mreg_linear_yreg_logistic_est(beta0 = 1,
                                                         beta1 = 2,
                                                         beta2 = 3,
                                                         theta1 = 4,
                                                         theta2 = NULL,
                                                         theta3 = 6,
                                                         theta4 = 7,
                                                         sigma_sq = 8))
            expect_error(
                calc_myreg_mreg_linear_yreg_logistic_est(beta0 = 1,
                                                         beta1 = 2,
                                                         beta2 = 3,
                                                         theta1 = 4,
                                                         theta2 = 5,
                                                         theta3 = NULL,
                                                         theta4 = 7,
                                                         sigma_sq = 8))
            expect_error(
                calc_myreg_mreg_linear_yreg_logistic_est(beta0 = 1,
                                                         beta1 = 2,
                                                         beta2 = 3,
                                                         theta1 = 4,
                                                         theta2 = 5,
                                                         theta3 = NULL,
                                                         theta4 = 7,
                                                         sigma_sq = 8:9))
        })
    })
    ## Note that this function does not require a model object and easy to test.
    describe("calc_myreg_mreg_linear_yreg_logistic_est (NULL cvar)", {
        est_fun <-
            calc_myreg_mreg_linear_yreg_logistic_est(beta0 = 1,
                                                     beta1 = 2,
                                                     beta2 = NULL,
                                                     theta1 = 4,
                                                     theta2 = 5,
                                                     theta3 = 6,
                                                     theta4 = NULL,
                                                     sigma_sq = 8)
        it("returns a function", {
            expect_equal(class(est_fun),
                         "function")
        })
        it("returns a function that takes a0, a1, m_cde, and c_cond", {
            expect_equal(names(formals(est_fun)),
                         c("a0","a1","m_cde","c_cond"))
        })
        it("returns a function with parameters in the enslosing environment", {
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "beta0"), 1)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "beta1"), 2)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "beta2"), NULL)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "theta1"), 4)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "theta2"), 5)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "theta3"), 6)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "theta4"), NULL)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "sigma_sq"), 8)
        })
        it("returns a function that errors given inconsistent c_cond", {
            expect_vector(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = NULL))
            expect_error(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1))
            expect_error(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:2))
            expect_error(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:3))
            expect_error(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:4))
        })
    })
    ##
    describe("calc_myreg_mreg_linear_yreg_logistic_est (one cvar)", {
        est_fun <-
            calc_myreg_mreg_linear_yreg_logistic_est(beta0 = 1,
                                                     beta1 = 2,
                                                     beta2 = 3,
                                                     theta1 = 4,
                                                     theta2 = 5,
                                                     theta3 = 6,
                                                     theta4 = 7,
                                                     sigma_sq = 8)
        it("returns a function", {
            expect_equal(class(est_fun),
                         "function")
        })
        it("returns a function that takes a0, a1, m_cde, and c_cond", {
            expect_equal(names(formals(est_fun)),
                         c("a0","a1","m_cde","c_cond"))
        })
        it("returns a function with parameters in the enslosing environment", {
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "beta0"), 1)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "beta1"), 2)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "beta2"), 3)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "theta1"), 4)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "theta2"), 5)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "theta3"), 6)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "theta4"), 7)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "sigma_sq"), 8)
        })
        it("returns a function that errors given inconsistent c_cond", {
            expect_error(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = NULL))
            expect_vector(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1))
            expect_error(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:2))
            expect_error(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:3))
            expect_error(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:4))
        })
    })
    ##
    describe("calc_myreg_mreg_linear_yreg_logistic_est (three cvar)", {
        est_fun <-
            calc_myreg_mreg_linear_yreg_logistic_est(beta0 = 1,
                                                     beta1 = 2,
                                                     beta2 = 3:5,
                                                     theta1 = 4,
                                                     theta2 = 5,
                                                     theta3 = 6,
                                                     theta4 = 7:9,
                                                     sigma_sq = 8)
        it("returns a function", {
            expect_equal(class(est_fun),
                         "function")
        })
        it("returns a function that takes a0, a1, m_cde, and c_cond", {
            expect_equal(names(formals(est_fun)),
                         c("a0","a1","m_cde","c_cond"))
        })
        it("returns a function with parameters in the enslosing environment", {
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "beta0"), 1)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "beta1"), 2)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "beta2"), 3:5)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "theta1"), 4)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "theta2"), 5)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "theta3"), 6)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "theta4"), 7:9)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "sigma_sq"), 8)
        })
        it("returns a function that errors given inconsistent c_cond", {
            expect_error(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = NULL))
            expect_error(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1))
            expect_error(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:2))
            expect_vector(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:3))
            expect_error(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:4))
        })
    })
})


###
### Standard error estimation function constructor
################################################################################

describe("calc_myreg_mreg_linear_yreg_logistic_se", {

})
