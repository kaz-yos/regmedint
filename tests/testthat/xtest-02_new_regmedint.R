################################################################################
###
##
## Created on: 2020-03-12
## Author: Kazuki Yoshida
################################################################################

## Load testthat in case this is run in isolation.
library(testthat)
library(survival)
library(tidyverse)


###
### Internal constructor
################################################################################

describe("new_regmedint", {
    it("returns an object with a class of regmedint", {
        res_new_regmedint <- new_regmedint(
            NULL
        )
        expect_equal(class(res_new_regmedint),
                     c("regmedint", "list"))
    })
})
