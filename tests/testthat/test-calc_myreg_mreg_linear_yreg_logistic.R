################################################################################
### Tests for regression-based causal mediation analysis
##
## Created on: 2020-03-11
## Author: Kazuki Yoshida
################################################################################

## Load testthat in case this is run in isolation.
library(testthat)
library(tidyverse)


###
### Effect estimation
################################################################################

test_that("calc_myreg_mreg_linear_yreg_logistic_est works given coef", {

})


###
### Standard error estimation
################################################################################

test_that("calc_myreg_mreg_linear_yreg_logistic_se works given coef and vcov", {

})


###
### Helper functions
################################################################################

test_that("variance estimates for sigma^2 is extracted", {

    data(pbc)

    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L))


    ## glm
    glm_fit <- glm(formula = alk.phos ~ trt + bili,
                   family  = gaussian(link = "identity"),
                   data    = pbc_cc)
    ## summary(glm_fit)$dispersion has no se.
    expect_equal(dim(vcov(glm_fit)),
                 ## 4x4 including the dispersion parameter.
                 c(4,4))

    ## lm
    lm_fit <- lm(formula = alk.phos ~ trt + bili,
                 data    = pbc_cc)
    expect_equal(dim(vcov(lm_fit)),
                 ## 4x4 including the dispersion parameter.
                 c(4,4))


})
