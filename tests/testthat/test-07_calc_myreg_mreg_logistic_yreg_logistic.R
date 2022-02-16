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
### Tests for calc_myreg_mreg_logistic_yreg_logistic
################################################################################

data(pbc)
## Missing data should be warned in validate_args()
pbc_cc <- pbc[complete.cases(pbc),] %>%
    mutate(male = if_else(sex == "m", 1L, 0L),
           ## Combine transplant and death for testing purpose
           status = if_else(status == 0, 0L, 1L),
           bili_bin = if_else(bili > median(bili), 1L, 0L),
           alk_phos = alk.phos)

describe("calc_myreg_mreg_logistic_yreg_logistic logistic no interaction", {

    describe("calc_myreg_mreg_logistic_yreg_logistic logistic no interaction (NULL cvar)", {
        mreg_fit <- fit_mreg(mreg = "logistic",
                             data = pbc_cc,
                             avar = "trt",
                             mvar = "bili_bin",
                             cvar = NULL)
        yreg_fit <- fit_yreg(yreg = "logistic",
                             data = pbc_cc,
                             yvar = "spiders",
                             avar = "trt",
                             mvar = "bili_bin",
                             cvar = NULL,
                             interaction = FALSE,
                             eventvar = NULL)
        myreg_funs <-
            calc_myreg_mreg_logistic_yreg_logistic(mreg = "logistic",
                                                   mreg_fit = mreg_fit,
                                                   yreg = "logistic",
                                                   yreg_fit = yreg_fit,
                                                   avar = "trt",
                                                   mvar = "bili_bin",
                                                   cvar = NULL,
                                                   interaction = FALSE,
                                                   emm_ac_mreg = NULL,
                                                   emm_ac_yreg = NULL,
                                                   emm_mc_yreg = NULL)
        ##
        it("returns a list of two functions", {
            expect_equal(class(myreg_funs),
                         "list")
            expect_equal(length(myreg_funs),
                         2)
        })
        it("returns functions that take 4 arguments", {
            expect_equal(names(formals(myreg_funs[[1]])),
                         c("a0","a1","m_cde","c_cond"))
            expect_equal(names(formals(myreg_funs[[2]])),
                         c("a0","a1","m_cde","c_cond"))
        })
        it("returns functions that return named vector of effect estimates", {
            expect_equal(names(myreg_funs[[1]](1,2,3,NULL)),
                         c("cde",
                           "pnde","tnie",
                           "tnde","pnie",
                           "te",
                           "pm"))
            expect_equal(names(myreg_funs[[2]](1,2,3,NULL)),
                         c("se_cde",
                           "se_pnde","se_tnie",
                           "se_tnde","se_pnie",
                           "se_te",
                           "se_pm"))
        })
        it("returns functions that error on inconsistent c_cond", {
            expect_error(myreg_funs[[1]](1,2,3,4), regexp = "c_cond")
            expect_error(myreg_funs[[2]](1,2,3,4), regexp = "c_cond")
            expect_error(myreg_funs[[1]](1,2,3,c(4,5,6)), regexp = "c_cond")
            expect_error(myreg_funs[[2]](1,2,3,c(4,5,6)), regexp = "c_cond")
        })
    })
    ##
    describe("calc_myreg_mreg_logistic_yreg_logistic logistic no interaction (1 cvar)", {
        mreg_fit <- fit_mreg(mreg = "logistic",
                             data = pbc_cc,
                             avar = "trt",
                             mvar = "bili_bin",
                             cvar = c("age"))
        yreg_fit <- fit_yreg(yreg = "logistic",
                             data = pbc_cc,
                             yvar = "spiders",
                             avar = "trt",
                             mvar = "bili_bin",
                             cvar = c("age"),
                             interaction = FALSE,
                             eventvar = NULL)
        myreg_funs <-
            calc_myreg_mreg_logistic_yreg_logistic(mreg = "logistic",
                                                   mreg_fit = mreg_fit,
                                                   yreg = "logistic",
                                                   yreg_fit = yreg_fit,
                                                   avar = "trt",
                                                   mvar = "bili_bin",
                                                   cvar = c("age"),
                                                   interaction = FALSE,
                                                   emm_ac_mreg = NULL,
                                                   emm_ac_yreg = NULL,
                                                   emm_mc_yreg = NULL)
        
        # add EMM
        myreg_funs_EMM1 <-
          calc_myreg_mreg_logistic_yreg_logistic(mreg = "logistic",
                                                 mreg_fit = mreg_fit,
                                                 yreg = "logistic",
                                                 yreg_fit = yreg_fit,
                                                 avar = "trt",
                                                 mvar = "bili_bin",
                                                 cvar = c("age"),
                                                 interaction = FALSE,
                                                 emm_ac_mreg = c("age"),
                                                 emm_ac_yreg = NULL,
                                                 emm_mc_yreg = NULL)
        
        myreg_funs_EMM2 <-
          calc_myreg_mreg_logistic_yreg_logistic(mreg = "logistic",
                                                 mreg_fit = mreg_fit,
                                                 yreg = "logistic",
                                                 yreg_fit = yreg_fit,
                                                 avar = "trt",
                                                 mvar = "bili_bin",
                                                 cvar = c("age"),
                                                 interaction = FALSE,
                                                 emm_ac_mreg = NULL,
                                                 emm_ac_yreg = c("age"),
                                                 emm_mc_yreg = NULL)
        
        myreg_funs_EMM3 <-
          calc_myreg_mreg_logistic_yreg_logistic(mreg = "logistic",
                                                 mreg_fit = mreg_fit,
                                                 yreg = "logistic",
                                                 yreg_fit = yreg_fit,
                                                 avar = "trt",
                                                 mvar = "bili_bin",
                                                 cvar = c("age"),
                                                 interaction = FALSE,
                                                 emm_ac_mreg = NULL,
                                                 emm_ac_yreg = NULL,
                                                 emm_mc_yreg = c("age"))
        
        
        ##
        it("returns a list of two functions", {
            expect_equal(class(myreg_funs),
                         "list")
            expect_equal(length(myreg_funs),
                         2)
        })
        it("returns functions that take 4 arguments", {
            expect_equal(names(formals(myreg_funs[[1]])),
                         c("a0","a1","m_cde","c_cond"))
            expect_equal(names(formals(myreg_funs[[2]])),
                         c("a0","a1","m_cde","c_cond"))
        })
        it("returns functions that return named vector of effect estimates", {
            expect_equal(names(myreg_funs[[1]](1,2,3,4)),
                         c("cde",
                           "pnde","tnie",
                           "tnde","pnie",
                           "te",
                           "pm"))
            expect_equal(names(myreg_funs[[2]](1,2,3,4)),
                         c("se_cde",
                           "se_pnde","se_tnie",
                           "se_tnde","se_pnie",
                           "se_te",
                           "se_pm"))
        })
        it("returns functions that error on inconsistent c_cond", {
            expect_error(myreg_funs[[1]](1,2,3,NULL), regexp = "c_cond")
            expect_error(myreg_funs[[2]](1,2,3,NULL), regexp = "c_cond")
            expect_error(myreg_funs[[1]](1,2,3,c(4,5,6)), regexp = "c_cond")
            expect_error(myreg_funs[[2]](1,2,3,c(4,5,6)), regexp = "c_cond")
        })
        
        # check functions with EMM differ from no EMM 
        it("check functions with EMM differ from no EMM", {
          expect_false(isTRUE(all.equal(myreg_funs, myreg_funs_EMM1)))
          expect_false(isTRUE(all.equal(myreg_funs_EMM1, myreg_funs_EMM2)))
          expect_false(isTRUE(all.equal(myreg_funs_EMM2, myreg_funs_EMM3)))
        })
    })
    ##
    describe("calc_myreg_mreg_logistic_yreg_logistic logistic no interaction (3 cvar)", {
        mreg_fit <- fit_mreg(mreg = "logistic",
                             data = pbc_cc,
                             avar = "trt",
                             mvar = "bili_bin",
                             cvar = c("age","male","stage"))
        yreg_fit <- fit_yreg(yreg = "logistic",
                             data = pbc_cc,
                             yvar = "spiders",
                             avar = "trt",
                             mvar = "bili_bin",
                             cvar = c("age","male","stage"),
                             interaction = FALSE,
                             eventvar = NULL)
        myreg_funs <-
            calc_myreg_mreg_logistic_yreg_logistic(mreg = "logistic",
                                                   mreg_fit = mreg_fit,
                                                   yreg = "logistic",
                                                   yreg_fit = yreg_fit,
                                                   avar = "trt",
                                                   mvar = "bili_bin",
                                                   cvar = c("age","male","stage"),
                                                   interaction = FALSE,
                                                   emm_ac_mreg = NULL,
                                                   emm_ac_yreg = NULL,
                                                   emm_mc_yreg = NULL)
        ##
        it("returns a list of two functions", {
            expect_equal(class(myreg_funs),
                         "list")
            expect_equal(length(myreg_funs),
                         2)
        })
        it("returns functions that take 4 arguments", {
            expect_equal(names(formals(myreg_funs[[1]])),
                         c("a0","a1","m_cde","c_cond"))
            expect_equal(names(formals(myreg_funs[[2]])),
                         c("a0","a1","m_cde","c_cond"))
        })
        it("returns functions that return named vector of effect estimates", {
            expect_equal(names(myreg_funs[[1]](1,2,3,c(4,5,6))),
                         c("cde",
                           "pnde","tnie",
                           "tnde","pnie",
                           "te",
                           "pm"))
            expect_equal(names(myreg_funs[[2]](1,2,3,c(4,5,6))),
                         c("se_cde",
                           "se_pnde","se_tnie",
                           "se_tnde","se_pnie",
                           "se_te",
                           "se_pm"))
        })
        it("returns functions that error on inconsistent c_cond", {
            expect_error(myreg_funs[[1]](1,2,3,NULL), regexp = "c_cond")
            expect_error(myreg_funs[[2]](1,2,3,NULL), regexp = "c_cond")
            expect_error(myreg_funs[[1]](1,2,3,4), regexp = "c_cond")
            expect_error(myreg_funs[[2]](1,2,3,4), regexp = "c_cond")
        })
    })
    describe("calc_myreg_mreg_logistic_yreg_logistic logistic no interaction (methodological correctness)", {
        mreg_fit <- fit_mreg(mreg = "logistic",
                             data = pbc_cc,
                             avar = "trt",
                             mvar = "bili_bin",
                             cvar = c("age","male","stage"))
        yreg_fit <- fit_yreg(yreg = "logistic",
                             data = pbc_cc,
                             yvar = "spiders",
                             avar = "trt",
                             mvar = "bili_bin",
                             cvar = c("age","male","stage"),
                             interaction = FALSE,
                             eventvar = NULL)
        myreg_funs <-
            calc_myreg_mreg_logistic_yreg_logistic(mreg = "logistic",
                                                   mreg_fit = mreg_fit,
                                                   yreg = "logistic",
                                                   yreg_fit = yreg_fit,
                                                   avar = "trt",
                                                   mvar = "bili_bin",
                                                   cvar = c("age","male","stage"),
                                                   interaction = FALSE,
                                                   emm_ac_mreg = NULL,
                                                   emm_ac_yreg = NULL,
                                                   emm_mc_yreg = NULL)
        it("returns functions where cde does not depend on m_cde", {
            expect_equal(myreg_funs[[1]](1,2,-3,c(4,5,6))["cde"],
                         myreg_funs[[1]](1,2,+3,c(4,5,6))["cde"])
            expect_equal(myreg_funs[[2]](1,2,-3,c(4,5,6))["cde"],
                         myreg_funs[[2]](1,2,+3,c(4,5,6))["cde"])
        })
        it("returns functions where nde do no depend on c_cond", {
            expect_equal(myreg_funs[[1]](1,2,-3,-1 * c(4,5,6))[c("pnde","tnde")],
                         myreg_funs[[1]](1,2,+3,+2 * c(4,5,6))[c("pnde","tnde")])
            expect_equal(myreg_funs[[2]](1,2,-3,-1 * c(4,5,6))[c("pnde","tnde")],
                         myreg_funs[[2]](1,2,+3,+2 * c(4,5,6))[c("pnde","tnde")])
        })
        it("returns functions where direct effects match up", {
            expect_equal(unname(myreg_funs[[1]](1,2,3,c(4,5,6))["cde"]),
                         unname(myreg_funs[[1]](1,2,3,c(4,5,6))["pnde"]))
            expect_equal(unname(myreg_funs[[2]](1,2,3,c(4,5,6))["cde"]),
                         unname(myreg_funs[[2]](1,2,3,c(4,5,6))["tnde"]))
        })
        it("returns functions where indirect effects match up", {
            expect_equal(unname(myreg_funs[[1]](1,2,3,c(4,5,6))["tnie"]),
                         unname(myreg_funs[[1]](1,2,3,c(4,5,6))["pnie"]))
        })
        it("returns functions where total effect is nde+nie", {
            expect_equal(unname(myreg_funs[[1]](1,2,3,c(4,5,6))["te"]),
                         ## Pearl decomposition
                         unname(myreg_funs[[1]](1,2,3,c(4,5,6))["pnde"]) +
                         unname(myreg_funs[[1]](1,2,3,c(4,5,6))["tnie"]))
            expect_equal(unname(myreg_funs[[1]](1,2,3,c(4,5,6))["te"]),
                         ## The other decomposition
                         unname(myreg_funs[[1]](1,2,3,c(4,5,6))["tnde"]) +
                         unname(myreg_funs[[1]](1,2,3,c(4,5,6))["pnie"]))
            ##
            expect_equal(unname(myreg_funs[[2]](1,2,3,c(4,5,6))["te"]),
                         ## Pearl decomposition
                         unname(myreg_funs[[2]](1,2,3,c(4,5,6))["pnde"]) +
                         unname(myreg_funs[[2]](1,2,3,c(4,5,6))["tnie"]))
            expect_equal(unname(myreg_funs[[2]](1,2,3,c(4,5,6))["te"]),
                         ## The other decomposition
                         unname(myreg_funs[[2]](1,2,3,c(4,5,6))["tnde"]) +
                         unname(myreg_funs[[2]](1,2,3,c(4,5,6))["pnie"]))
        })
        it("returns functions where pm is calculated from natural effects correctly", {
            log_nde <- unname(myreg_funs[[1]](1,2,3,c(4,5,6))["pnde"])
            log_nie <- unname(myreg_funs[[1]](1,2,3,c(4,5,6))["tnie"])
            expect_equal(unname(myreg_funs[[1]](1,2,3,c(4,5,6))["pm"]),
                         ## VanderWeele 2015. p48.
            ((exp(log_nde) * (exp(log_nie) - 1)) /
             ((exp(log_nde) * exp(log_nie)) - 1)))
        })
    })
})
##
describe("calc_myreg_mreg_logistic_yreg_logistic logistic interaction", {

    describe("calc_myreg_mreg_logistic_yreg_logistic logistic interaction (NULL cvar)", {
        mreg_fit <- fit_mreg(mreg = "logistic",
                             data = pbc_cc,
                             avar = "trt",
                             mvar = "bili_bin",
                             cvar = NULL)
        yreg_fit <- fit_yreg(yreg = "logistic",
                             data = pbc_cc,
                             yvar = "spiders",
                             avar = "trt",
                             mvar = "bili_bin",
                             cvar = NULL,
                             interaction = TRUE,
                             eventvar = NULL)
        myreg_funs <-
            calc_myreg_mreg_logistic_yreg_logistic(mreg = "logistic",
                                                   mreg_fit = mreg_fit,
                                                   yreg = "logistic",
                                                   yreg_fit = yreg_fit,
                                                   avar = "trt",
                                                   mvar = "bili_bin",
                                                   cvar = NULL,
                                                   interaction = TRUE,
                                                   emm_ac_mreg = NULL,
                                                   emm_ac_yreg = NULL,
                                                   emm_mc_yreg = NULL)
        ##
        it("returns a list of two functions", {
            expect_equal(class(myreg_funs),
                         "list")
            expect_equal(length(myreg_funs),
                         2)
        })
        it("returns functions that take 4 arguments", {
            expect_equal(names(formals(myreg_funs[[1]])),
                         c("a0","a1","m_cde","c_cond"))
            expect_equal(names(formals(myreg_funs[[2]])),
                         c("a0","a1","m_cde","c_cond"))
        })
        it("returns functions that return named vector of effect estimates", {
            expect_equal(names(myreg_funs[[1]](1,2,3,NULL)),
                         c("cde",
                           "pnde","tnie",
                           "tnde","pnie",
                           "te",
                           "pm"))
            expect_equal(names(myreg_funs[[2]](1,2,3,NULL)),
                         c("se_cde",
                           "se_pnde","se_tnie",
                           "se_tnde","se_pnie",
                           "se_te",
                           "se_pm"))
        })
        it("returns functions that error on inconsistent c_cond", {
            expect_error(myreg_funs[[1]](1,2,3,4), regexp = "c_cond")
            expect_error(myreg_funs[[2]](1,2,3,4), regexp = "c_cond")
            expect_error(myreg_funs[[1]](1,2,3,c(4,5,6)), regexp = "c_cond")
            expect_error(myreg_funs[[2]](1,2,3,c(4,5,6)), regexp = "c_cond")
        })
    })
    ##
    describe("calc_myreg_mreg_logistic_yreg_logistic logistic interaction (1 cvar)", {
        mreg_fit <- fit_mreg(mreg = "logistic",
                             data = pbc_cc,
                             avar = "trt",
                             mvar = "bili_bin",
                             cvar = c("age"))
        yreg_fit <- fit_yreg(yreg = "logistic",
                             data = pbc_cc,
                             yvar = "spiders",
                             avar = "trt",
                             mvar = "bili_bin",
                             cvar = c("age"),
                             interaction = TRUE,
                             eventvar = NULL)
        myreg_funs <-
            calc_myreg_mreg_logistic_yreg_logistic(mreg = "logistic",
                                                   mreg_fit = mreg_fit,
                                                   yreg = "logistic",
                                                   yreg_fit = yreg_fit,
                                                   avar = "trt",
                                                   mvar = "bili_bin",
                                                   cvar = c("age"),
                                                   interaction = TRUE,
                                                   emm_ac_mreg = NULL,
                                                   emm_ac_yreg = NULL,
                                                   emm_mc_yreg = NULL)
        ##
        it("returns a list of two functions", {
            expect_equal(class(myreg_funs),
                         "list")
            expect_equal(length(myreg_funs),
                         2)
        })
        it("returns functions that take 4 arguments", {
            expect_equal(names(formals(myreg_funs[[1]])),
                         c("a0","a1","m_cde","c_cond"))
            expect_equal(names(formals(myreg_funs[[2]])),
                         c("a0","a1","m_cde","c_cond"))
        })
        it("returns functions that return named vector of effect estimates", {
            expect_equal(names(myreg_funs[[1]](1,2,3,4)),
                         c("cde",
                           "pnde","tnie",
                           "tnde","pnie",
                           "te",
                           "pm"))
            expect_equal(names(myreg_funs[[2]](1,2,3,4)),
                         c("se_cde",
                           "se_pnde","se_tnie",
                           "se_tnde","se_pnie",
                           "se_te",
                           "se_pm"))
        })
        it("returns functions that error on inconsistent c_cond", {
            expect_error(myreg_funs[[1]](1,2,3,NULL), regexp = "c_cond")
            expect_error(myreg_funs[[2]](1,2,3,NULL), regexp = "c_cond")
            expect_error(myreg_funs[[1]](1,2,3,c(4,5,6)), regexp = "c_cond")
            expect_error(myreg_funs[[2]](1,2,3,c(4,5,6)), regexp = "c_cond")
        })
    })
    ##
    describe("calc_myreg_mreg_logistic_yreg_logistic logistic interaction (3 cvar)", {
        mreg_fit <- fit_mreg(mreg = "logistic",
                             data = pbc_cc,
                             avar = "trt",
                             mvar = "bili_bin",
                             cvar = c("age","male","stage"))
        yreg_fit <- fit_yreg(yreg = "logistic",
                             data = pbc_cc,
                             yvar = "spiders",
                             avar = "trt",
                             mvar = "bili_bin",
                             cvar = c("age","male","stage"),
                             interaction = TRUE,
                             eventvar = NULL)
        myreg_funs <-
            calc_myreg_mreg_logistic_yreg_logistic(mreg = "logistic",
                                                   mreg_fit = mreg_fit,
                                                   yreg = "logistic",
                                                   yreg_fit = yreg_fit,
                                                   avar = "trt",
                                                   mvar = "bili_bin",
                                                   cvar = c("age","male","stage"),
                                                   interaction = TRUE,
                                                   emm_ac_mreg = NULL,
                                                   emm_ac_yreg = NULL,
                                                   emm_mc_yreg = NULL)
        ##
        it("returns a list of two functions", {
            expect_equal(class(myreg_funs),
                         "list")
            expect_equal(length(myreg_funs),
                         2)
        })
        it("returns functions that take 4 arguments", {
            expect_equal(names(formals(myreg_funs[[1]])),
                         c("a0","a1","m_cde","c_cond"))
            expect_equal(names(formals(myreg_funs[[2]])),
                         c("a0","a1","m_cde","c_cond"))
        })
        it("returns functions that return named vector of effect estimates", {
            expect_equal(names(myreg_funs[[1]](1,2,3,c(4,5,6))),
                         c("cde",
                           "pnde","tnie",
                           "tnde","pnie",
                           "te",
                           "pm"))
            expect_equal(names(myreg_funs[[2]](1,2,3,c(4,5,6))),
                         c("se_cde",
                           "se_pnde","se_tnie",
                           "se_tnde","se_pnie",
                           "se_te",
                           "se_pm"))
        })
        it("returns functions that error on inconsistent c_cond", {
            expect_error(myreg_funs[[1]](1,2,3,NULL), regexp = "c_cond")
            expect_error(myreg_funs[[2]](1,2,3,NULL), regexp = "c_cond")
            expect_error(myreg_funs[[1]](1,2,3,4), regexp = "c_cond")
            expect_error(myreg_funs[[2]](1,2,3,4), regexp = "c_cond")
        })
    })
    describe("calc_myreg_mreg_logistic_yreg_logistic logistic interaction (methodological correctness)", {
        mreg_fit <- fit_mreg(mreg = "logistic",
                             data = pbc_cc,
                             avar = "trt",
                             mvar = "bili_bin",
                             cvar = c("age","male","stage"))
        yreg_fit <- fit_yreg(yreg = "logistic",
                             data = pbc_cc,
                             yvar = "spiders",
                             avar = "trt",
                             mvar = "bili_bin",
                             cvar = c("age","male","stage"),
                             interaction = TRUE,
                             eventvar = NULL)
        ## Sign of the interaction coefficient is important.
        theta3 <- coef(yreg_fit)["trt:bili_bin"]
        beta1 <- coef(mreg_fit)[c("trt")]
        beta2 <- coef(mreg_fit)[c("age","male","stage")]
        myreg_funs <-
            calc_myreg_mreg_logistic_yreg_logistic(mreg = "logistic",
                                                   mreg_fit = mreg_fit,
                                                   yreg = "logistic",
                                                   yreg_fit = yreg_fit,
                                                   avar = "trt",
                                                   mvar = "bili_bin",
                                                   cvar = c("age","male","stage"),
                                                   interaction = TRUE,
                                                   emm_ac_mreg = NULL,
                                                   emm_ac_yreg = NULL,
                                                   emm_mc_yreg = NULL)
        it("returns functions where cde depends on m_cde", {
            ## Positive (a1 - a0)
            if (theta3 > 0) {
                ## Increasing in m_cde
                expect_gt(myreg_funs[[1]](1,2,+3,c(4,5,6))["cde"],
                          myreg_funs[[1]](1,2,-3,c(4,5,6))["cde"])
            } else if (theta3 < 0) {
                ## Decreasing in m_cde
                expect_lt(myreg_funs[[1]](1,2,+3,c(4,5,6))["cde"],
                          myreg_funs[[1]](1,2,-3,c(4,5,6))["cde"])
            }
        })
        it("returns functions where cde does not depend on c_cond", {
            expect_equal(myreg_funs[[1]](1,2,+3,-1 * c(4,5,6))["cde"],
                         myreg_funs[[1]](1,2,+3,+2 * c(4,5,6))["cde"])
            expect_equal(myreg_funs[[2]](1,2,+3,-1 * c(4,5,6))["cde"],
                         myreg_funs[[2]](1,2,+3,+2 * c(4,5,6))["cde"])
        })
        it("returns functions where natural effects depend on c_cond", {
            ## NDE is seems monotone. NIE may be conincidental.
            c_cond <- c(4,5,6)
            beta2_c1 <- sum(-1 * beta2 * c_cond)
            beta2_c2 <- sum(+2 * beta2 * c_cond)
            ## Positive (a1 - a0)
            if (theta3 > 0) {
                ## Increasing in beta2_c
                if (beta2_c1 - beta2_c2 > 0) {
                    expect_gt(myreg_funs[[1]](1,2,+3,-1 * c_cond)["pnde"],
                              myreg_funs[[1]](1,2,+3,+2 * c_cond)["pnde"])
                    expect_gt(myreg_funs[[1]](1,2,+3,-1 * c_cond)["tnde"],
                              myreg_funs[[1]](1,2,+3,+2 * c_cond)["tnde"])
                    expect_gt(myreg_funs[[1]](1,2,+3,-1 * c_cond)["pnie"],
                              myreg_funs[[1]](1,2,+3,+2 * c_cond)["pnie"])
                    expect_gt(myreg_funs[[1]](1,2,+3,-1 * c_cond)["tnie"],
                              myreg_funs[[1]](1,2,+3,+2 * c_cond)["tnie"])
                    expect_gt(myreg_funs[[1]](1,2,+3,-1 * c_cond)["te"],
                              myreg_funs[[1]](1,2,+3,+2 * c_cond)["te"])
                } else if (beta2_c1 - beta2_c2 < 0) {
                    expect_lt(myreg_funs[[1]](1,2,+3,-1 * c_cond)["pnde"],
                              myreg_funs[[1]](1,2,+3,+2 * c_cond)["pnde"])
                    expect_lt(myreg_funs[[1]](1,2,+3,-1 * c_cond)["tnde"],
                              myreg_funs[[1]](1,2,+3,+2 * c_cond)["tnde"])
                    expect_lt(myreg_funs[[1]](1,2,+3,-1 * c_cond)["pnie"],
                              myreg_funs[[1]](1,2,+3,+2 * c_cond)["pnie"])
                    expect_lt(myreg_funs[[1]](1,2,+3,-1 * c_cond)["tnie"],
                              myreg_funs[[1]](1,2,+3,+2 * c_cond)["tnie"])
                    expect_lt(myreg_funs[[1]](1,2,+3,-1 * c_cond)["te"],
                              myreg_funs[[1]](1,2,+3,+2 * c_cond)["te"])
                }
            } else if (theta3 < 0) {
                ## Decreasing in beta2_c
                if (beta2_c1 - beta2_c2 > 0) {
                    expect_lt(myreg_funs[[1]](1,2,+3,-1 * c_cond)["pnde"],
                              myreg_funs[[1]](1,2,+3,+2 * c_cond)["pnde"])
                    expect_lt(myreg_funs[[1]](1,2,+3,-1 * c_cond)["tnde"],
                              myreg_funs[[1]](1,2,+3,+2 * c_cond)["tnde"])
                    expect_lt(myreg_funs[[1]](1,2,+3,-1 * c_cond)["te"],
                              myreg_funs[[1]](1,2,+3,+2 * c_cond)["te"])
                } else if (beta2_c1 - beta2_c2 < 0) {
                    expect_gt(myreg_funs[[1]](1,2,+3,-1 * c_cond)["pnde"],
                              myreg_funs[[1]](1,2,+3,+2 * c_cond)["pnde"])
                    expect_gt(myreg_funs[[1]](1,2,+3,-1 * c_cond)["tnde"],
                              myreg_funs[[1]](1,2,+3,+2 * c_cond)["tnde"])
                    expect_gt(myreg_funs[[1]](1,2,+3,-1 * c_cond)["te"],
                              myreg_funs[[1]](1,2,+3,+2 * c_cond)["te"])
                }
            }
        })
        it("returns functions where natural effects relate correctly", {
            ## Positive (a1 - a0)
            if (theta3 * beta1 > 0) {
                ## Increasing in a0 (pure) -> a1 (total)
                expect_gt(unname(myreg_funs[[1]](1,2,3,c(4,5,6))["tnde"]),
                          unname(myreg_funs[[1]](1,2,3,c(4,5,6))["pnde"]))
                expect_gt(unname(myreg_funs[[1]](1,2,3,c(4,5,6))["tnie"]),
                          unname(myreg_funs[[1]](1,2,3,c(4,5,6))["pnie"]))
            } else if (theta3 * beta1 < 0) {
                ## Decreasing in a0 (pure) -> a1 (total)
                expect_lt(unname(myreg_funs[[1]](1,2,3,c(4,5,6))["tnde"]),
                          unname(myreg_funs[[1]](1,2,3,c(4,5,6))["pnde"]))
                expect_lt(unname(myreg_funs[[1]](1,2,3,c(4,5,6))["tnie"]),
                          unname(myreg_funs[[1]](1,2,3,c(4,5,6))["pnie"]))
            }
        })
        it("returns functions where total effect is nde+nie", {
            expect_equal(unname(myreg_funs[[1]](1,2,3,c(4,5,6))["te"]),
                         ## Pearl decomposition
                         unname(myreg_funs[[1]](1,2,3,c(4,5,6))["pnde"]) +
                         unname(myreg_funs[[1]](1,2,3,c(4,5,6))["tnie"]))
            expect_equal(unname(myreg_funs[[1]](1,2,3,c(4,5,6))["te"]),
                         ## The other decomposition
                         unname(myreg_funs[[1]](1,2,3,c(4,5,6))["tnde"]) +
                         unname(myreg_funs[[1]](1,2,3,c(4,5,6))["pnie"]))
            ##
            expect_equal(unname(myreg_funs[[2]](1,2,3,c(4,5,6))["te"]),
                         ## Pearl decomposition
                         unname(myreg_funs[[2]](1,2,3,c(4,5,6))["pnde"]) +
                         unname(myreg_funs[[2]](1,2,3,c(4,5,6))["tnie"]))
            expect_equal(unname(myreg_funs[[2]](1,2,3,c(4,5,6))["te"]),
                         ## The other decomposition
                         unname(myreg_funs[[2]](1,2,3,c(4,5,6))["tnde"]) +
                         unname(myreg_funs[[2]](1,2,3,c(4,5,6))["pnie"]))
        })
        it("returns functions where pm is calculated from natural effects correctly", {
            log_nde <- unname(myreg_funs[[1]](1,2,3,c(4,5,6))["pnde"])
            log_nie <- unname(myreg_funs[[1]](1,2,3,c(4,5,6))["tnie"])
            expect_equal(unname(myreg_funs[[1]](1,2,3,c(4,5,6))["pm"]),
                         ## VanderWeele 2015. p48.
            ((exp(log_nde) * (exp(log_nie) - 1)) /
             ((exp(log_nde) * exp(log_nie)) - 1)))
        })
    })
})


###
### Effect estimation function constructor
################################################################################
## No need to repeat this part for Poisson etc because *_est functions work
## on extracted parameters not on model objects.

## FIXME: Most of these functionalities should be factored out to be shared
## among calc_myreg_mreg_*_yreg_*_est.
describe("calc_myreg_mreg_logistic_yreg_logistic_est function factory", {
    describe("calc_myreg_mreg_logistic_yreg_logistic_est (error handling)", {
        it("errors given inconsistent beta2 and theta4", {
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_est(beta0 = 1,
                                                           beta1 = 2,
                                                           beta2 = 3,
                                                           beta3 = NULL,
                                                           theta0 = 0,
                                                           theta1 = 4,
                                                           theta2 = 5,
                                                           theta3 = 6,
                                                           theta4 = NULL,
                                                           theta5 = NULL,
                                                           theta6 = NULL))
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_est(beta0 = 1,
                                                           beta1 = 2,
                                                           beta2 = 1:2,
                                                           beta3 = NULL,
                                                           theta0 = 0,
                                                           theta1 = 4,
                                                           theta2 = 5,
                                                           theta3 = 6,
                                                           theta4 = NULL,
                                                           theta5 = NULL,
                                                           theta6 = NULL))
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_est(beta0 = 1,
                                                           beta1 = 2,
                                                           beta2 = NULL,
                                                           beta3 = NULL,
                                                           theta0 = 0,
                                                           theta1 = 4,
                                                           theta2 = 5,
                                                           theta3 = 6,
                                                           theta4 = 7,
                                                           theta5 = NULL,
                                                           theta6 = NULL))
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_est(beta0 = 1,
                                                           beta1 = 2,
                                                           beta2 = NULL,
                                                           beta3 = NULL,
                                                           theta0 = 0,
                                                           theta1 = 4,
                                                           theta2 = 5,
                                                           theta3 = 6,
                                                           theta4 = 7:8,
                                                           theta5 = NULL,
                                                           theta6 = NULL))
        })
        it("errors given vector inputs in arguments other than beta2 and theta4", {
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_est(beta0 = 1:2,
                                                           beta1 = 2,
                                                           beta2 = 3,
                                                           beta3 = NULL,
                                                           theta0 = 0,
                                                           theta1 = 4,
                                                           theta2 = 5,
                                                           theta3 = 6,
                                                           theta4 = 7,
                                                           theta5 = NULL,
                                                           theta6 = NULL))
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_est(beta0 = 1,
                                                           beta1 = 2:3,
                                                           beta2 = 3,
                                                           beta3 = NULL,
                                                           theta0 = 0,
                                                           theta1 = 4,
                                                           theta2 = 5,
                                                           theta3 = 6,
                                                           theta4 = 7,
                                                           theta5 = NULL,
                                                           theta6 = NULL))
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_est(beta0 = 1,
                                                           beta1 = 2,
                                                           beta2 = 3,
                                                           beta3 = NULL,
                                                           theta0 = 0,
                                                           theta1 = 4:5,
                                                           theta2 = 5,
                                                           theta3 = 6,
                                                           theta4 = 7,
                                                           theta5 = NULL,
                                                           theta6 = NULL))
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_est(beta0 = 1,
                                                           beta1 = 2,
                                                           beta2 = 3,
                                                           beta3 = NULL,
                                                           theta0 = 0,
                                                           theta1 = 4,
                                                           theta2 = 5:6,
                                                           theta3 = 6,
                                                           theta4 = 7,
                                                           theta5 = NULL,
                                                           theta6 = NULL))
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_est(beta0 = 1,
                                                           beta1 = 2,
                                                           beta2 = 3,
                                                           beta3 = NULL,
                                                           theta0 = 0,
                                                           theta1 = 4,
                                                           theta2 = 5,
                                                           theta3 = 6:7,
                                                           theta4 = 7,
                                                           theta5 = NULL,
                                                           theta6 = NULL))
        })
        it("errors given NULL inputs in arguments other than beta2 and theta4", {
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_est(beta0 = NULL,
                                                           beta1 = 2,
                                                           beta2 = 3,
                                                           beta3 = NULL,
                                                           theta0 = 0,
                                                           theta1 = 4,
                                                           theta2 = 5,
                                                           theta3 = 6,
                                                           theta4 = 7,
                                                           theta5 = NULL,
                                                           theta6 = NULL))
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_est(beta0 = 1,
                                                           beta1 = NULL,
                                                           beta2 = 3,
                                                           beta3 = NULL,
                                                           theta0 = 0,
                                                           theta1 = 4,
                                                           theta2 = 5,
                                                           theta3 = 6,
                                                           theta4 = 7,
                                                           theta5 = NULL,
                                                           theta6 = NULL))
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_est(beta0 = 1,
                                                           beta1 = 2,
                                                           beta2 = 3,
                                                           beta3 = NULL,
                                                           theta0 = 0,
                                                           theta1 = NULL,
                                                           theta2 = 5,
                                                           theta3 = 6,
                                                           theta4 = 7,
                                                           theta5 = NULL,
                                                           theta6 = NULL))
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_est(beta0 = 1,
                                                           beta1 = 2,
                                                           beta2 = 3,
                                                           beta3 = NULL,
                                                           theta0 = 0,
                                                           theta1 = 4,
                                                           theta2 = NULL,
                                                           theta3 = 6,
                                                           theta4 = 7,
                                                           theta5 = NULL,
                                                           theta6 = NULL))
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_est(beta0 = 1,
                                                           beta1 = 2,
                                                           beta2 = 3,
                                                           beta3 = NULL,
                                                           theta0 = 0,
                                                           theta1 = 4,
                                                           theta2 = 5,
                                                           theta3 = NULL,
                                                           theta4 = 7,
                                                           theta5 = NULL,
                                                           theta6 = NULL))
        })
    })
    ## Note that this function does not require a model object and easy to test.
    describe("calc_myreg_mreg_logistic_yreg_logistic_est (NULL cvar)", {
        est_fun <-
            calc_myreg_mreg_logistic_yreg_logistic_est(beta0 = 1,
                                                       beta1 = 2,
                                                       beta2 = NULL,
                                                       beta3 = NULL,
                                                       theta0 = 0,
                                                       theta1 = 4,
                                                       theta2 = 5,
                                                       theta3 = 6,
                                                       theta4 = NULL,
                                                       theta5 = NULL,
                                                       theta6 = NULL)
        it("returns a function", {
            expect_equal(class(est_fun),
                         "function")
        })
        it("returns a function that takes a0, a1, m_cde, c_cond", {
            expect_equal(names(formals(est_fun)),
                         c("a0","a1","m_cde","c_cond"))
        })
        it("returns a function with parameters in the enslosing environment", {
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "beta0"), 1)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "beta1"), 2)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "beta2"), NULL)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "theta0"), 0)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "theta1"), 4)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "theta2"), 5)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "theta3"), 6)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "theta4"), NULL)
        })
        it("returns a function that errors given inconsistent c_cond", {
            expect_vector(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = NULL))
            expect_error(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1))
            expect_error(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:2))
            expect_error(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:3))
            expect_error(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:4))
        })
        it("returns a function that gives a numeric vector without NA", {
            res <- est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = NULL)
            expect_true(is.vector(res))
            expect_true(is.numeric(res))
            expect_true(all(!is.na(res)))
        })
    })
    ##
    describe("calc_myreg_mreg_logistic_yreg_logistic_est (one cvar)", {
        est_fun <-
            calc_myreg_mreg_logistic_yreg_logistic_est(beta0 = 1,
                                                       beta1 = 2,
                                                       beta2 = 3,
                                                       beta3 = NULL,
                                                       theta0 = 0,
                                                       theta1 = 4,
                                                       theta2 = 5,
                                                       theta3 = 6,
                                                       theta4 = 7,
                                                       theta5 = NULL,
                                                       theta6 = NULL)
        it("returns a function", {
            expect_equal(class(est_fun),
                         "function")
        })
        it("returns a function that takes a0, a1, m_cde, c_cond", {
            expect_equal(names(formals(est_fun)),
                         c("a0","a1","m_cde","c_cond"))
        })
        it("returns a function with parameters in the enslosing environment", {
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "beta0"), 1)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "beta1"), 2)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "beta2"), 3)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "theta0"), 0)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "theta1"), 4)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "theta2"), 5)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "theta3"), 6)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "theta4"), 7)
        })
        it("returns a function that errors given inconsistent c_cond", {
            expect_error(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = NULL))
            expect_vector(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1))
            expect_error(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:2))
            expect_error(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:3))
            expect_error(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:4))
        })
        it("returns a function that gives a numeric vector without NA", {
            res <- est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1)
            expect_true(is.vector(res))
            expect_true(is.numeric(res))
            expect_true(all(!is.na(res)))
        })
    })
    ##
    describe("calc_myreg_mreg_logistic_yreg_logistic_est (three cvar)", {
        est_fun <-
            calc_myreg_mreg_logistic_yreg_logistic_est(beta0 = 1,
                                                       beta1 = 2,
                                                       beta2 = 3:5,
                                                       beta3 = NULL,
                                                       theta0 = 0,
                                                       theta1 = 4,
                                                       theta2 = 5,
                                                       theta3 = 6,
                                                       theta4 = 7:9,
                                                       theta5 = NULL,
                                                       theta6 = NULL)
        it("returns a function", {
            expect_equal(class(est_fun),
                         "function")
        })
        it("returns a function that takes a0, a1, m_cde, c_cond", {
            expect_equal(names(formals(est_fun)),
                         c("a0","a1","m_cde","c_cond"))
        })
        it("returns a function with parameters in the enslosing environment", {
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "beta0"), 1)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "beta1"), 2)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "beta2"), 3:5)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "theta0"), 0)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "theta1"), 4)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "theta2"), 5)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "theta3"), 6)
            expect_equal(rlang::env_get(rlang::fn_env(est_fun), nm = "theta4"), 7:9)
        })
        it("returns a function that errors given inconsistent c_cond", {
            expect_error(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = NULL))
            expect_error(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1))
            expect_error(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:2))
            expect_vector(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:3))
            expect_error(est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:4))
        })
        it("returns a function that gives a numeric vector without NA", {
            res <- est_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:3)
            expect_true(is.vector(res))
            expect_true(is.numeric(res))
            expect_true(all(!is.na(res)))
        })
    })
})


###
### Standard error estimation function constructor
################################################################################
## No need to repeat this part for Poisson etc because *_se functions work
## on extracted parameters not on model objects.

describe("calc_myreg_mreg_logistic_yreg_logistic_se function factory", {
    describe("calc_myreg_mreg_logistic_yreg_logistic_se (error handling)", {
        it("errors given inconsistent beta2 and theta4", {
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_se(beta0 = 1,
                                                          beta1 = 2,
                                                          beta2 = 3,
                                                          beta3 = NULL,
                                                          theta0 = 0,
                                                          theta1 = 4,
                                                          theta2 = 5,
                                                          theta3 = 6,
                                                          theta4 = NULL,
                                                          theta5 = NULL,
                                                          theta6 = NULL,
                                                          Sigma_beta = diag(1, nrow = 3, ncol = 3),
                                                          Sigma_theta = diag(2, nrow = 4, ncol = 4)))
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_se(beta0 = 1,
                                                          beta1 = 2,
                                                          beta2 = 1:2,
                                                          beta3 = NULL,
                                                          theta0 = 0,
                                                          theta1 = 4,
                                                          theta2 = 5,
                                                          theta3 = 6,
                                                          theta4 = NULL,
                                                          theta5 = NULL,
                                                          theta6 = NULL,
                                                          Sigma_beta = diag(1, nrow = 3, ncol = 3),
                                                          Sigma_theta = diag(2, nrow = 4, ncol = 4)))
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_se(beta0 = 1,
                                                          beta1 = 2,
                                                          beta2 = NULL,
                                                          beta3 = NULL,
                                                          theta0 = 0,
                                                          theta1 = 4,
                                                          theta2 = 5,
                                                          theta3 = 6,
                                                          theta4 = 7,
                                                          theta5 = NULL,
                                                          theta6 = NULL,
                                                          Sigma_beta = diag(1, nrow = 3, ncol = 3),
                                                          Sigma_theta = diag(2, nrow = 4, ncol = 4)))
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_se(beta0 = 1,
                                                          beta1 = 2,
                                                          beta2 = NULL,
                                                          beta3 = NULL,
                                                          theta0 = 0,
                                                          theta1 = 4,
                                                          theta2 = 5,
                                                          theta3 = 6,
                                                          theta4 = 7:8,
                                                          theta5 = NULL,
                                                          theta6 = NULL,
                                                          Sigma_beta = diag(1, nrow = 3, ncol = 3),
                                                          Sigma_theta = diag(2, nrow = 4, ncol = 4)))
        })
        it("errors given vector inputs in arguments other than beta2 and theta4", {
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_se(beta0 = 1:2,
                                                          beta1 = 2,
                                                          beta2 = 3,
                                                          beta3 = NULL,
                                                          theta0 = 0,
                                                          theta1 = 4,
                                                          theta2 = 5,
                                                          theta3 = 6,
                                                          theta4 = 7,
                                                          theta5 = NULL,
                                                          theta6 = NULL,
                                                          Sigma_beta = diag(1, nrow = 3, ncol = 3),
                                                          Sigma_theta = diag(2, nrow = 4, ncol = 4)))
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_se(beta0 = 1,
                                                          beta1 = 2:3,
                                                          beta2 = 3,
                                                          beta3 = NULL,
                                                          theta0 = 0,
                                                          theta1 = 4,
                                                          theta2 = 5,
                                                          theta3 = 6,
                                                          theta4 = 7,
                                                          theta5 = NULL,
                                                          theta6 = NULL,
                                                          Sigma_beta = diag(1, nrow = 3, ncol = 3),
                                                          Sigma_theta = diag(2, nrow = 4, ncol = 4)))
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_se(beta0 = 1,
                                                          beta1 = 2,
                                                          beta2 = 3,
                                                          beta3 = NULL,
                                                          theta0 = 0,
                                                          theta1 = 4:5,
                                                          theta2 = 5,
                                                          theta3 = 6,
                                                          theta4 = 7,
                                                          theta5 = NULL,
                                                          theta6 = NULL,
                                                          Sigma_beta = diag(1, nrow = 3, ncol = 3),
                                                          Sigma_theta = diag(2, nrow = 4, ncol = 4)))
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_se(beta0 = 1,
                                                          beta1 = 2,
                                                          beta2 = 3,
                                                          beta3 = NULL,
                                                          theta0 = 0,
                                                          theta1 = 4,
                                                          theta2 = 5:6,
                                                          theta3 = 6,
                                                          theta4 = 7,
                                                          theta5 = NULL,
                                                          theta6 = NULL,
                                                          Sigma_beta = diag(1, nrow = 3, ncol = 3),
                                                          Sigma_theta = diag(2, nrow = 4, ncol = 4)))
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_se(beta0 = 1,
                                                          beta1 = 2,
                                                          beta2 = 3,
                                                          beta3 = NULL,
                                                          theta0 = 0,
                                                          theta1 = 4,
                                                          theta2 = 5,
                                                          theta3 = 6:7,
                                                          theta4 = 7,
                                                          theta5 = NULL,
                                                          theta6 = NULL,
                                                          Sigma_beta = diag(1, nrow = 3, ncol = 3),
                                                          Sigma_theta = diag(2, nrow = 4, ncol = 4)))
        })
        it("errors given NULL inputs in arguments other than beta2 and theta4", {
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_se(beta0 = NULL,
                                                          beta1 = 2,
                                                          beta2 = 3,
                                                          beta3 = NULL,
                                                          theta0 = 0,
                                                          theta1 = 4,
                                                          theta2 = 5,
                                                          theta3 = 6,
                                                          theta4 = 7,
                                                          theta5 = NULL,
                                                          theta6 = NULL,
                                                          Sigma_beta = diag(1, nrow = 3, ncol = 3),
                                                          Sigma_theta = diag(2, nrow = 4, ncol = 4)))
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_se(beta0 = 1,
                                                          beta1 = NULL,
                                                          beta2 = 3,
                                                          beta3 = NULL,
                                                          theta0 = 0,
                                                          theta1 = 4,
                                                          theta2 = 5,
                                                          theta3 = 6,
                                                          theta4 = 7,
                                                          theta5 = NULL,
                                                          theta6 = NULL,
                                                          Sigma_beta = diag(1, nrow = 3, ncol = 3),
                                                          Sigma_theta = diag(2, nrow = 4, ncol = 4)))
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_se(beta0 = 1,
                                                          beta1 = 2,
                                                          beta2 = 3,
                                                          beta3 = NULL,
                                                          theta0 = 0,
                                                          theta1 = NULL,
                                                          theta2 = 5,
                                                          theta3 = 6,
                                                          theta4 = 7,
                                                          theta5 = NULL,
                                                          theta6 = NULL,
                                                          Sigma_beta = diag(1, nrow = 3, ncol = 3),
                                                          Sigma_theta = diag(2, nrow = 4, ncol = 4)))
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_se(beta0 = 1,
                                                          beta1 = 2,
                                                          beta2 = 3,
                                                          beta3 = NULL,
                                                          theta0 = 0,
                                                          theta1 = 4,
                                                          theta2 = NULL,
                                                          theta3 = 6,
                                                          theta4 = 7,
                                                          theta5 = NULL,
                                                          theta6 = NULL,
                                                          Sigma_beta = diag(1, nrow = 3, ncol = 3),
                                                          Sigma_theta = diag(2, nrow = 4, ncol = 4)))
            expect_error(
                calc_myreg_mreg_logistic_yreg_logistic_se(beta0 = 1,
                                                          beta1 = 2,
                                                          beta2 = 3,
                                                          beta3 = NULL,
                                                          theta0 = 0,
                                                          theta1 = 4,
                                                          theta2 = 5,
                                                          theta3 = NULL,
                                                          theta4 = 7,
                                                          theta5 = NULL,
                                                          theta6 = NULL,
                                                          Sigma_beta = diag(1, nrow = 3, ncol = 3),
                                                          Sigma_theta = diag(2, nrow = 4, ncol = 4)))
        })
    })
    ## Note that this function does not require a model object and easy to tse.
    describe("calc_myreg_mreg_logistic_yreg_logistic_se (NULL cvar)", {
        se_fun <-
            calc_myreg_mreg_logistic_yreg_logistic_se(beta0 = 0.1,
                                                      beta1 = 0.2,
                                                      beta2 = NULL,
                                                      beta3 = NULL,
                                                      theta0 = 0,
                                                      theta1 = 0.4,
                                                      theta2 = 0.5,
                                                      theta3 = 0.6,
                                                      theta4 = NULL,
                                                      theta5 = NULL,
                                                      theta6 = NULL,
                                                      Sigma_beta = diag(1, nrow = 2, ncol = 2),
                                                      Sigma_theta = diag(1.2, nrow = 4, ncol = 4))
        it("returns a function", {
            expect_equal(class(se_fun),
                         "function")
        })
        it("returns a function that takes a0, a1, m_cde, c_cond", {
            expect_equal(names(formals(se_fun)),
                         c("a0","a1","m_cde","c_cond"))
        })
        it("returns a function with parameters in the enslosing environment", {
            expect_equal(rlang::env_get(rlang::fn_env(se_fun), nm = "beta0"), 0.1)
            expect_equal(rlang::env_get(rlang::fn_env(se_fun), nm = "beta1"), 0.2)
            expect_equal(rlang::env_get(rlang::fn_env(se_fun), nm = "beta2"), NULL)
            expect_equal(rlang::env_get(rlang::fn_env(se_fun), nm = "theta0"), 0.0)
            expect_equal(rlang::env_get(rlang::fn_env(se_fun), nm = "theta1"), 0.4)
            expect_equal(rlang::env_get(rlang::fn_env(se_fun), nm = "theta2"), 0.5)
            expect_equal(rlang::env_get(rlang::fn_env(se_fun), nm = "theta3"), 0.6)
            expect_equal(rlang::env_get(rlang::fn_env(se_fun), nm = "theta4"), NULL)
        })
        it("returns a function that errors given inconsistent c_cond", {
            expect_vector(se_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = NULL))
            expect_error(se_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1))
            expect_error(se_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:2))
            expect_error(se_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:3))
            expect_error(se_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:4))
        })
        it("returns a function that gives a numeric vector without NA", {
            res <- se_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = NULL)
            expect_true(is.vector(res))
            expect_true(is.numeric(res))
            expect_true(all(!is.na(res)))
        })
    })
    ##
    describe("calc_myreg_mreg_logistic_yreg_logistic_se (one cvar)", {
        se_fun <-
            calc_myreg_mreg_logistic_yreg_logistic_se(beta0 = 0.1,
                                                      beta1 = 0.2,
                                                      beta2 = 0.3,
                                                      beta3 = NULL,
                                                      theta0 = 0,
                                                      theta1 = 0.4,
                                                      theta2 = 0.5,
                                                      theta3 = 0.6,
                                                      theta4 = 0.7,
                                                      theta5 = NULL,
                                                      theta6 = NULL,
                                                      Sigma_beta = diag(1, nrow = 3, ncol = 3),
                                                      Sigma_theta = diag(1.1, nrow = 5, ncol = 5))
        it("returns a function", {
            expect_equal(class(se_fun),
                         "function")
        })
        it("returns a function that takes a0, a1, m_cde, c_cond", {
            expect_equal(names(formals(se_fun)),
                         c("a0","a1","m_cde","c_cond"))
        })
        it("returns a function with parameters in the enslosing environment", {
            expect_equal(rlang::env_get(rlang::fn_env(se_fun), nm = "beta0"), 0.1)
            expect_equal(rlang::env_get(rlang::fn_env(se_fun), nm = "beta1"), 0.2)
            expect_equal(rlang::env_get(rlang::fn_env(se_fun), nm = "beta2"), 0.3)
            expect_equal(rlang::env_get(rlang::fn_env(se_fun), nm = "theta0"), 0)
            expect_equal(rlang::env_get(rlang::fn_env(se_fun), nm = "theta1"), 0.4)
            expect_equal(rlang::env_get(rlang::fn_env(se_fun), nm = "theta2"), 0.5)
            expect_equal(rlang::env_get(rlang::fn_env(se_fun), nm = "theta3"), 0.6)
            expect_equal(rlang::env_get(rlang::fn_env(se_fun), nm = "theta4"), 0.7)
        })
        it("returns a function that errors given inconsistent c_cond", {
            expect_error(se_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = NULL))
            expect_vector(se_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1))
            expect_error(se_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:2))
            expect_error(se_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:3))
            expect_error(se_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:4))
        })
        it("returns a function that gives a numeric vector without NA", {
            res <- se_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1)
            expect_true(is.vector(res))
            expect_true(is.numeric(res))
            expect_true(all(!is.na(res)))
        })
    })
    ##
    describe("calc_myreg_mreg_logistic_yreg_logistic_se (three cvar)", {
        se_fun <-
            calc_myreg_mreg_logistic_yreg_logistic_se(beta0 = 0.1,
                                                      beta1 = 0.2,
                                                      beta2 = 3:5/10,
                                                      beta3 = NULL,
                                                      theta0 = 0,
                                                      theta1 = 0.4,
                                                      theta2 = 0.5,
                                                      theta3 = 0.6,
                                                      theta4 = 7:9/10,
                                                      theta5 = NULL,
                                                      theta6 = NULL,
                                                      Sigma_beta = diag(1, nrow = 5, ncol = 5),
                                                      Sigma_theta = diag(1.1, nrow = 7, ncol = 7))
        it("returns a function", {
            expect_equal(class(se_fun),
                         "function")
        })
        it("returns a function that takes a0, a1, m_cde, c_cond", {
            expect_equal(names(formals(se_fun)),
                         c("a0","a1","m_cde","c_cond"))
        })
        it("returns a function with parameters in the enslosing environment", {
            expect_equal(rlang::env_get(rlang::fn_env(se_fun), nm = "beta0"), 0.1)
            expect_equal(rlang::env_get(rlang::fn_env(se_fun), nm = "beta1"), 0.2)
            expect_equal(rlang::env_get(rlang::fn_env(se_fun), nm = "beta2"), 3:5/10)
            expect_equal(rlang::env_get(rlang::fn_env(se_fun), nm = "theta0"), 0)
            expect_equal(rlang::env_get(rlang::fn_env(se_fun), nm = "theta1"), 0.4)
            expect_equal(rlang::env_get(rlang::fn_env(se_fun), nm = "theta2"), 0.5)
            expect_equal(rlang::env_get(rlang::fn_env(se_fun), nm = "theta3"), 0.6)
            expect_equal(rlang::env_get(rlang::fn_env(se_fun), nm = "theta4"), 7:9/10)
        })
        it("returns a function that errors given inconsistent c_cond", {
            expect_error(se_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = NULL))
            expect_error(se_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1))
            expect_error(se_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:2))
            expect_vector(se_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:3))
            expect_error(se_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:4))
        })
        it("returns a function that gives a numeric vector without NA", {
            res <- se_fun(a0 = 0, a1 = 1, m_cde = 0, c_cond = 1:3)
            expect_true(is.vector(res))
            expect_true(is.numeric(res))
            expect_true(all(!is.na(res)))
        })
    })
})
