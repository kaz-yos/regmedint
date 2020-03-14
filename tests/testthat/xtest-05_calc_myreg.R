################################################################################
### Specifications for the main mediation analysis function
##
## Created on: 2020-03-12
## Author: Kazuki Yoshida
################################################################################

## Load testthat in case this is run in isolation.
library(testthat)
library(survival)
library(tidyverse)


###
### calc_myreg
################################################################################

describe("calc_myreg", {
    res_calc_myreg <- calc_myreg(

    )
    it("returns a list of two functions", {
        expect_equal(length(res_calc_myreg),
                     2)
        expect_equal(class(res_calc_myreg[[1]]),
                     "function")
        expect_equal(class(res_calc_myreg[[2]]),
                     "function")
    })
    it("returns a list of two functions with four argumens (a0, a1, m_cde, c_cond)", {
        expect_equal(formals(res_calc_myreg[[1]]),
                     c("a0","a1","m_cde","c_cond"))
        expect_equal(formals(res_calc_myreg[[2]]),
                     c("a0","a1","m_cde","c_cond"))
    })

})
