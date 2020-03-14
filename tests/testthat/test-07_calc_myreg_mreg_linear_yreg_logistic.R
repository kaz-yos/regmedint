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
        it("returns a function that return named vector of effect estimates", {
            expect(names(myreg_funs[[1]](1,2,3,4)),
                   c("cde","pnde","tnie","tnde","pnie","te","pm"))
            expect(names(myreg_funs[[2]](1,2,3,4)),
                   c("se_cde","se_pnde","se_tnie","se_tnde","se_pnie","se_te","se_pm"))
        })
    })
})


###
### Effect estimation function constructor
################################################################################

describe("calc_myreg_mreg_linear_yreg_logistic_est", {

    est_fun <- calc_myreg_mreg_linear_yreg_logistic_est()
    res <- est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = c(1,2,3))

    it("returns a numeric vector of length 7 with appropriate named elements", {
        expect_equal(length(res),
                     7)
        expect_equal(names(res),
                     c("cde","pnde","tnie","tnde","pnie","te","pm"))
    })
})


###
### Standard error estimation function constructor
################################################################################

describe("calc_myreg_mreg_linear_yreg_logistic_se", {

})
