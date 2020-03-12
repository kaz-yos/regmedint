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
