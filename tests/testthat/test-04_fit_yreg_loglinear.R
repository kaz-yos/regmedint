################################################################################
### Tests for internal functions
##
## Created on: 2020-03-09
## Author: Kazuki Yoshida
################################################################################

## Load testthat in case this is run in isolation.
library(sandwich)
library(geepack)
library(testthat)
library(survival)
library(tidyverse)

###
### Internal function for yreg model fitting (loglinear)
################################################################################

## Zou 2004
## Am J Epidemiol. 2004 Apr 1;159(7):702-6.
## A modified poisson regression approach to prospective studies with binary data.
## https://www.ncbi.nlm.nih.gov/pubmed/15033648
##
## Zeileis 2006
## J Stat Softw. 2006;16:1-16.
## Object-oriented Computation of Sandwich Estimators
## https://www.jstatsoft.org/article/view/v016i09


describe("Modified Poisson regression in R", {

    ## Examine whether glm/sandwich implementation
    ## is comparable to geepack implementation.
    ## https://rpubs.com/kaz_yos/epi204-lab4

    data(pbc)
    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L),
               id = factor(id))

    it("implemented in glm/sandwich matches geepack exchangeable (no covariates)", {
        glm_fit0 <- glm(formula = spiders ~ trt + bili + trt:bili,
                        family = poisson(link = "log"),
                        data = pbc_cc)

        geeglm_fit0 <- geeglm(formula   = spiders ~ trt + bili + trt:bili,
                              family    = poisson(link = "log"),
                              id        = id,
                              data      = pbc_cc,
                              corstr    = "exchangeable")

        expect_equal(coef(glm_fit0),
                     coef(geeglm_fit0),
                     ## Tolerance for absolute differences
                     tolerance = 0.00000001, scale = 1)

        expect_equal(sandwich::sandwich(glm_fit0),
                     vcov(geeglm_fit0),
                     ## Tolerance for absolute differences
                     tolerance = 0.000001, scale = 1)

        expect_equal(sqrt(diag(sandwich::sandwich(glm_fit0))),
                     sqrt(diag(vcov(geeglm_fit0))),
                     ## Tolerance for absolute differences
                     tolerance = 0.000001, scale = 1)

    })

    it("implemented in glm/sandwich  geepack exchangeable (three covariates)", {
        glm_fit3 <- glm(formula = spiders ~ trt + bili + trt:bili + age + male + stage,
                        family = poisson(link = "log"),
                        data = pbc_cc)

        geeglm_fit3 <- geeglm(formula   = spiders ~ trt + bili + trt:bili + age + male + stage,
                              family    = poisson(link = "log"),
                              id        = id,
                              data      = pbc_cc,
                              corstr    = "exchangeable")

        expect_equal(coef(glm_fit3),
                     coef(geeglm_fit3),
                     ## Tolerance for absolute differences
                     tolerance = 0.0000001, scale = 1)

        expect_equal(sandwich::sandwich(glm_fit3),
                     vcov(geeglm_fit3),
                     ## Tolerance for absolute differences
                     tolerance = 0.00001, scale = 1)

        expect_equal(sqrt(diag(sandwich::sandwich(glm_fit3))),
                     sqrt(diag(vcov(geeglm_fit3))),
                     ## Tolerance for absolute differences
                     tolerance = 0.00001, scale = 1)
    })
    
    
})


describe("fit_yreg loglinear as modified poisson (no interaction)", {

    data(pbc)
    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L))

    it("fits a correct model with no covariates", {
        ## No covariates
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
        ref_fit0 <- glm(formula = spiders ~ trt + bili,
                        family = poisson(link = "log"),
                        data = pbc_cc)
        ## Same classes
        expect_equal(class(yreg_fit0),
                     c("regmedint_mod_poisson", class(ref_fit0)))
        ## Same formula
        expect_equal(as.character(yreg_fit0$call$formula),
                     as.character(ref_fit0$call$formula))
        ## Same coef
        expect_equal(coef(yreg_fit0),
                     coef(ref_fit0))
        ## Robust vcov
        expect_equal(vcov(yreg_fit0),
                     sandwich::sandwich(ref_fit0))
        ## Summary should use robust vcov
        expect_equal(coef(summary(yreg_fit0))[,"Std. Error"],
                     sqrt(diag(sandwich::sandwich(ref_fit0))))
    })

    it("fits a correct model with one covariate", {
        ## One covariates
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
        ref_fit1 <- glm(formula = spiders ~ trt + bili + age,
                        family = poisson(link = "log"),
                        data = pbc_cc)
        ## Same classes
        expect_equal(class(yreg_fit1),
                     c("regmedint_mod_poisson", class(ref_fit1)))
        ## Same formula
        expect_equal(as.character(yreg_fit1$call$formula),
                     as.character(ref_fit1$call$formula))
        ## Same coef
        expect_equal(coef(yreg_fit1),
                     coef(ref_fit1))
        ## Robust vcov
        expect_equal(vcov(yreg_fit1),
                     sandwich::sandwich(ref_fit1))
        ## Summary should use robust vcov
        expect_equal(coef(summary(yreg_fit1))[,"Std. Error"],
                     sqrt(diag(sandwich::sandwich(ref_fit1))))
    })

    it("fits a correct model with three covariates", {
        ## Three covariates
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
        ref_fit3 <- glm(formula = spiders ~ trt + bili + age + male + stage,
                        family = poisson(link = "log"),
                        data = pbc_cc)
        ## Same classes
        expect_equal(class(yreg_fit3),
                     c("regmedint_mod_poisson", class(ref_fit3)))
        ## Same formula
        expect_equal(as.character(yreg_fit3$call$formula),
                     as.character(ref_fit3$call$formula))
        ## Same coef
        expect_equal(coef(yreg_fit3),
                     coef(ref_fit3))
        ## Robust vcov
        expect_equal(vcov(yreg_fit3),
                     sandwich::sandwich(ref_fit3))
        ## Summary should use robust vcov
        expect_equal(coef(summary(yreg_fit3))[,"Std. Error"],
                     sqrt(diag(sandwich::sandwich(ref_fit3))))
    })
    
    # only test when emm_ac_yreg and emm_mc_yreg are both not null:
    it("fits a correct model with three covariates, and non-null emm_ac_yreg and non-null EMM_M", {
        ## Three covariates
        yreg_fit6 <- fit_yreg(yreg = "loglinear",
                              data = pbc_cc,
                              yvar = "spiders",
                              avar = "trt",
                              mvar = "bili",
                              cvar = c("age","male","stage"),
                              emm_ac_yreg = c("age"),
                              emm_mc_yreg = c("male", "stage"),
                              interaction = FALSE,
                              eventvar = NULL)
        ref_fit6 <- glm(formula = spiders ~ trt + bili + age + male + stage + 
                            trt:age + bili:male + bili:stage,
                        family = poisson(link = "log"),
                        data = pbc_cc)
        ## Same classes
        expect_equal(class(yreg_fit6),
                     c("regmedint_mod_poisson", class(ref_fit6)))
        ## Same formula
        expect_equal(as.character(yreg_fit6$call$formula),
                     as.character(ref_fit6$call$formula))
        ## Same coef
        expect_equal(coef(yreg_fit6),
                     coef(ref_fit6))
        ## Robust vcov
        expect_equal(vcov(yreg_fit6),
                     sandwich::sandwich(ref_fit6))
        ## Summary should use robust vcov
        expect_equal(coef(summary(yreg_fit6))[,"Std. Error"],
                     sqrt(diag(sandwich::sandwich(ref_fit6))))
    })

})


describe("fit_yreg loglinear as modified poisson (interaction)", {

    data(pbc)

    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L))

    it("fits a correct model with no covariates", {
        ## No covariates
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
        ref_fit0 <- glm(formula = spiders ~ trt + bili + trt:bili,
                        family = poisson(link = "log"),
                        data = pbc_cc)
        ## Same classes
        expect_equal(class(yreg_fit0),
                     c("regmedint_mod_poisson", class(ref_fit0)))
        ## Same formula
        expect_equal(as.character(yreg_fit0$call$formula),
                     as.character(ref_fit0$call$formula))
        ## Same coef
        expect_equal(coef(yreg_fit0),
                     coef(ref_fit0))
        ## Robust vcov
        expect_equal(vcov(yreg_fit0),
                     sandwich::sandwich(ref_fit0))
        ## Summary should use robust vcov
        expect_equal(coef(summary(yreg_fit0))[,"Std. Error"],
                     sqrt(diag(sandwich::sandwich(ref_fit0))))
    })

    it("fits a correct model with one covariate", {
        ## One covariates
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
        ref_fit1 <- glm(formula = spiders ~ trt + bili + trt:bili + age,
                        family = poisson(link = "log"),
                        data = pbc_cc)
        ## Same classes
        expect_equal(class(yreg_fit1),
                     c("regmedint_mod_poisson", class(ref_fit1)))
        ## Same formula
        expect_equal(as.character(yreg_fit1$call$formula),
                     as.character(ref_fit1$call$formula))
        ## Same coef
        expect_equal(coef(yreg_fit1),
                     coef(ref_fit1))
        ## Robust vcov
        expect_equal(vcov(yreg_fit1),
                     sandwich::sandwich(ref_fit1))
        ## Summary should use robust vcov
        expect_equal(coef(summary(yreg_fit1))[,"Std. Error"],
                     sqrt(diag(sandwich::sandwich(ref_fit1))))
    })

    it("fits a correct model with three covariates", {
        ## Three covariates
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
        ref_fit3 <- glm(formula = spiders ~ trt + bili + trt:bili + age + male + stage,
                        family = poisson(link = "log"),
                        data = pbc_cc)
        ## Same classes
        expect_equal(class(yreg_fit3),
                     c("regmedint_mod_poisson", class(ref_fit3)))
        ## Same formula
        expect_equal(as.character(yreg_fit3$call$formula),
                     as.character(ref_fit3$call$formula))
        ## Same coef
        expect_equal(coef(yreg_fit3),
                     coef(ref_fit3))
        ## Robust vcov
        expect_equal(vcov(yreg_fit3),
                     sandwich::sandwich(ref_fit3))
        ## Summary should use robust vcov
        expect_equal(coef(summary(yreg_fit3))[,"Std. Error"],
                     sqrt(diag(sandwich::sandwich(ref_fit3))))
    })
    
    # only test when emm_ac_yreg and emm_mc_yreg are both not null:
    it("fits a correct model with three covariates, and non-null emm_ac_yreg and non-null emm_mc_yreg", {
        ## Three covariates
        yreg_fit6 <- fit_yreg(yreg = "loglinear",
                              data = pbc_cc,
                              yvar = "spiders",
                              avar = "trt",
                              mvar = "bili",
                              cvar = c("age","male","stage"),
                              emm_ac_yreg = c("age"),
                              emm_mc_yreg = c("male", "stage"),
                              interaction = TRUE,
                              eventvar = NULL)
        ref_fit6 <- glm(formula = spiders ~ trt + bili + trt:bili + age + male + stage + 
                            trt:age + bili:male + bili:stage,
                        family = poisson(link = "log"),
                        data = pbc_cc)
        ## Same classes
        expect_equal(class(yreg_fit6),
                     c("regmedint_mod_poisson", class(ref_fit6)))
        ## Same formula
        expect_equal(as.character(yreg_fit6$call$formula),
                     as.character(ref_fit6$call$formula))
        ## Same coef
        expect_equal(coef(yreg_fit6),
                     coef(ref_fit6))
        ## Robust vcov
        expect_equal(vcov(yreg_fit6),
                     sandwich::sandwich(ref_fit6))
        ## Summary should use robust vcov
        expect_equal(coef(summary(yreg_fit6))[,"Std. Error"],
                     sqrt(diag(sandwich::sandwich(ref_fit6))))
    })

})
