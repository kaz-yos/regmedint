################################################################################
### Tests for regmedint
##
## Created on: 2020-03-16
## Author: Kazuki Yoshida
################################################################################

library(testthat)
library(survival)
library(tidyverse)


###
### regmedint user interface
################################################################################

describe("regmedint", {

    data(pbc)
    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc %>%
        as_tibble() %>%
        ## Missing data should be warned in validate_args()
        drop_na() %>%
        mutate(male = if_else(sex == "m", 1L, 0L),
               ## Combine transplant and death for testing purpose
               status = if_else(status == 0, 0L, 1L),
               ##
               event = if_else(status == 1L, 1L, 0L),
               ## Binary mvar
               bili_bin = if_else(bili > median(bili), 1L, 0L),
               alk_phos = alk.phos)

    describe("validate_args (regmedint argument validation)", {
        it("rejects missing data in the variales of interest", {
            msg_missing <- "Missing is not allowed in variables of intrest! Consider multiple imputation."
            expect_error(pbc_cc %>%
                         mutate(alk_phos = NA) %>%
                         regmedint(data = .,
                                   yvar = "alk_phos",
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   a0 = 1,
                                   a1 = 2,
                                   m_cde = 0,
                                   c_cond = c(50,2,4),
                                   mreg = "linear",
                                   yreg = "linear",
                                   interaction = FALSE,
                                   casecontrol = FALSE,
                                   eventvar = NULL),
                         msg_missing)
            expect_error(pbc_cc %>%
                         mutate(trt = NA) %>%
                         regmedint(data = .,
                                   yvar = "alk_phos",
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   a0 = 1,
                                   a1 = 2,
                                   m_cde = 0,
                                   c_cond = c(50,2,4),
                                   mreg = "linear",
                                   yreg = "linear",
                                   interaction = FALSE,
                                   casecontrol = FALSE,
                                   eventvar = NULL),
                         msg_missing)
            expect_error(pbc_cc %>%
                         mutate(bili = NA) %>%
                         regmedint(data = .,
                                   yvar = "alk_phos",
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   a0 = 1,
                                   a1 = 2,
                                   m_cde = 0,
                                   c_cond = c(50,2,4),
                                   mreg = "linear",
                                   yreg = "linear",
                                   interaction = FALSE,
                                   casecontrol = FALSE,
                                   eventvar = NULL),
                         msg_missing)
            expect_error(pbc_cc %>%
                         mutate(age = NA) %>%
                         regmedint(data = .,
                                   yvar = "alk_phos",
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   a0 = 1,
                                   a1 = 2,
                                   m_cde = 0,
                                   c_cond = c(50,2,4),
                                   mreg = "linear",
                                   yreg = "linear",
                                   interaction = FALSE,
                                   casecontrol = FALSE,
                                   eventvar = NULL),
                         msg_missing)
        })
        ##
        it("allows missing data in unused variables", {
            expect_equal(pbc_cc %>%
                         mutate(sex = NA) %>%
                         regmedint(data = .,
                                   yvar = "alk.phos",
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = NULL,
                                   a0 = 1,
                                   a1 = 2,
                                   m_cde = 0,
                                   c_cond = NULL,
                                   mreg = "linear",
                                   yreg = "linear",
                                   interaction = FALSE,
                                   casecontrol = FALSE,
                                   eventvar = NULL),
                         regmedint(data = pbc_cc,
                                   yvar = "alk.phos",
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = NULL,
                                   a0 = 1,
                                   a1 = 2,
                                   m_cde = 0,
                                   c_cond = NULL,
                                   mreg = "linear",
                                   yreg = "linear",
                                   interaction = FALSE,
                                   casecontrol = FALSE,
                                   eventvar = NULL))
        })
        ##
        it("rejects factor variables", {
            msg_factor <- "Factors are not allowed! Use numeric variables only."
            expect_error(pbc_cc %>%
                         mutate(alk_phos = factor(alk_phos)) %>%
                         regmedint(data = ,
                                   yvar = "alk_phos",
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   a0 = 1,
                                   a1 = 2,
                                   m_cde = 0,
                                   c_cond = c(50,2,4),
                                   mreg = "linear",
                                   yreg = "linear",
                                   interaction = FALSE,
                                   casecontrol = FALSE,
                                   eventvar = NULL),
                         msg_factor)
            expect_error(pbc_cc %>%
                         mutate(trt = factor(trt)) %>%
                         regmedint(data = .,
                                   yvar = "alk_phos",
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   a0 = 1,
                                   a1 = 2,
                                   m_cde = 0,
                                   c_cond = c(50,2,4),
                                   mreg = "linear",
                                   yreg = "linear",
                                   interaction = FALSE,
                                   casecontrol = FALSE,
                                   eventvar = NULL),
                         msg_factor)
            expect_error(pbc_cc %>%
                         mutate(bili = factor(bili)) %>%
                         regmedint(data = .,
                                   yvar = "alk_phos",
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   a0 = 1,
                                   a1 = 2,
                                   m_cde = 0,
                                   c_cond = c(50,2,4),
                                   mreg = "linear",
                                   yreg = "linear",
                                   interaction = FALSE,
                                   casecontrol = FALSE,
                                   eventvar = NULL),
                         msg_factor)
            expect_error(pbc_cc %>%
                         mutate(status = factor(status)) %>%
                         regmedint(data = .,
                                   yvar = "alk_phos",
                                   avar = "trt",
                                   mvar = "bili",
                                   cvar = c("age","male","stage"),
                                   a0 = 1,
                                   a1 = 2,
                                   m_cde = 0,
                                   c_cond = c(50,2,4),
                                   mreg = "linear",
                                   yreg = "linear",
                                   interaction = FALSE,
                                   casecontrol = FALSE,
                                   eventvar = "status"),
                         msg_factor)
        })
    })

    describe("regmedint mreg linear yreg linear", {
        it("runs with zero cvar with no interaction", {
            fit_regmedint <- regmedint(data = pbc_cc,
                                       yvar = "alk.phos",
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       a0 = 1,
                                       a1 = 2,
                                       m_cde = 0,
                                       c_cond = NULL,
                                       mreg = "linear",
                                       yreg = "linear",
                                       interaction = FALSE,
                                       casecontrol = FALSE,
                                       eventvar = NULL)
        })
    })
    ##
    describe("regmedint mreg linear yreg logistic", {
        it("runs with zero cvar with no interaction", {
            fit_regmedint <- regmedint(data = pbc_cc,
                                       yvar = "spiders",
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       a0 = 1,
                                       a1 = 2,
                                       m_cde = 0,
                                       c_cond = NULL,
                                       mreg = "linear",
                                       yreg = "logistic",
                                       interaction = FALSE,
                                       casecontrol = FALSE,
                                       eventvar = NULL)
        })
    })
})
