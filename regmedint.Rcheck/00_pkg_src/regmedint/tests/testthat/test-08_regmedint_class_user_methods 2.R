################################################################################
### Tests for the user methods for the regmedint class
##
## Created on: 2020-03-16
## Author: Kazuki Yoshida
################################################################################

###
### regmedint class
################################################################################

describe("methods for regmedint", {

    data(pbc)
    ## Missing data should be warned in validate_args()
    pbc_cc <- pbc[complete.cases(pbc),] %>%
        mutate(male = if_else(sex == "m", 1L, 0L),
               ## Combine transplant and death for testing purpose
               status = if_else(status == 0, 0L, 1L),
               bili_bin = if_else(bili > median(bili), 1L, 0L),
               alk_phos = alk.phos)

    describe("methods for regmedint mreg linear yreg linear", {
        fit_regmedint <- regmedint(data = pbc_cc,
                                   yvar = "alk_phos",
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
        fit_regmedint_int <- regmedint(data = pbc_cc,
                                       yvar = "alk_phos",
                                       avar = "trt",
                                       mvar = "bili",
                                       cvar = NULL,
                                       a0 = 1,
                                       a1 = 2,
                                       m_cde = 0,
                                       c_cond = NULL,
                                       mreg = "linear",
                                       yreg = "linear",
                                       interaction = TRUE,
                                       casecontrol = FALSE,
                                       eventvar = NULL)
        ##
        describe("print.regmedint", {
            it("prints the mreg results", {
                expect_output(print(fit_regmedint),
                              deparse(fit_regmedint$mreg$call)[1],
                              fixed = TRUE)
            })
            it("prints the yreg results", {
                expect_output(print(fit_regmedint),
                              deparse(fit_regmedint$yreg$call)[1],
                              fixed = TRUE)
            })
            it("prints mediation analysis results with expected elements", {
                expect_output(print(fit_regmedint), "cde")
                expect_output(print(fit_regmedint), "pnde")
                expect_output(print(fit_regmedint), "tnie")
                expect_output(print(fit_regmedint), "tnde")
                expect_output(print(fit_regmedint), "pnie")
                expect_output(print(fit_regmedint), "te")
                expect_output(print(fit_regmedint), "pm")
            })
        })
        describe("summary.regmedint", {
            it("returns an object of summary_regmedint class", {
                expect_equal(class(summary(fit_regmedint))[[1]],
                             "summary_regmedint")
            })
            it("returns an object containing the mreg summary", {
                expect_equal(summary(fit_regmedint)$summary_mreg_fit,
                             summary(fit_regmedint$mreg_fit))
            })
            it("returns an object containing the yreg summary", {
                expect_equal(summary(fit_regmedint)$summary_yreg_fit,
                             summary(fit_regmedint$yreg_fit))
            })
            it("returns an object containing the myreg summary matrix", {
                expect_equal(class(summary(fit_regmedint)$summary_myreg)[[1]],
                             "matrix")
            })
            it("returns an object with appropriate columns", {
                expect_equal(colnames(summary(fit_regmedint)$summary_myreg),
                             c("est","se","Z","p","lower","upper"))
                ## expect_equivalent ignores attributes
                expect_equivalent(summary(fit_regmedint)$summary_myreg[,"est"],
                                  coef(fit_regmedint))
                expect_equivalent(summary(fit_regmedint)$summary_myreg[,"se"],
                                  sqrt(diag(vcov(fit_regmedint))))
                expect_equivalent(summary(fit_regmedint)$summary_myreg[,"lower"],
                                  confint(fit_regmedint)[,"lower"])
                expect_equivalent(summary(fit_regmedint)$summary_myreg[,"upper"],
                                  confint(fit_regmedint)[,"upper"])
            })
            it("returns an object, ignoring exponentiate = TRUE", {
                expect_equal(colnames(summary(fit_regmedint, exponentiate = TRUE)$summary_myreg),
                             c("est","se","Z","p","lower","upper"))
                ## expect_equivalent ignores attributes
                expect_equivalent(summary(fit_regmedint, exponentiate = TRUE)$summary_myreg[,"est"],
                                  coef(fit_regmedint))
                expect_equivalent(summary(fit_regmedint, exponentiate = TRUE)$summary_myreg[,"se"],
                                  sqrt(diag(vcov(fit_regmedint))))
                expect_equivalent(summary(fit_regmedint, exponentiate = TRUE)$summary_myreg[,"lower"],
                                  confint(fit_regmedint)[,"lower"])
                expect_equivalent(summary(fit_regmedint, exponentiate = TRUE)$summary_myreg[,"upper"],
                                  confint(fit_regmedint)[,"upper"])
            })
        })
        ##
        describe("methods for summary_regmedint", {
            describe("print.summary_regmedint", {
                ## Explicit printing within the function.
                ## No need to print the return value.
                it("prints the mreg results", {
                    expect_output(print(summary(fit_regmedint)),
                                  deparse(fit_regmedint$mreg$call)[1],
                                  fixed = TRUE)
                })
                it("prints the yreg results", {
                    expect_output(print(summary(fit_regmedint)),
                                  deparse(fit_regmedint$yreg$call)[1],
                                  fixed = TRUE)
                })
                it("prints mediation analysis results with expected elements", {
                    expect_output(print(summary(fit_regmedint)), "cde")
                    expect_output(print(summary(fit_regmedint)), "pnde")
                    expect_output(print(summary(fit_regmedint)), "tnie")
                    expect_output(print(summary(fit_regmedint)), "tnde")
                    expect_output(print(summary(fit_regmedint)), "pnie")
                    expect_output(print(summary(fit_regmedint)), "te")
                    expect_output(print(summary(fit_regmedint)), "pm")
                })
                it("prints mediation analysis results with expected elements exponentiated", {
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "cde")
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "pnde")
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "tnie")
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "tnde")
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "pnie")
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "te")
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "pm")
                })
                it("prints evaluation information", {
                    expect_output(print(summary(fit_regmedint)),
                                  "Evaluated at:",
                                  fixed = TRUE)
                    expect_output(print(summary(fit_regmedint)),
                                  "a1 (intervened value of avar) = ",
                                  fixed = TRUE)
                    expect_output(print(summary(fit_regmedint)),
                                  "a0 (reference value of avar)  = ",
                                  fixed = TRUE)
                    expect_output(print(summary(fit_regmedint)),
                                  "m_cde (intervend value of mvar for cde) = ",
                                  fixed = TRUE)
                    expect_output(print(summary(fit_regmedint)),
                                  "c_cond (covariate vector value) =",
                                  fixed = TRUE)
                    expect_output(print(summary(fit_regmedint)),
                                  "Note that effect estimates do not vary over m_cde and c_cond values when interaction = FALSE.",
                                  fixed = TRUE)
                    expect_output(print(summary(fit_regmedint_int)),
                                  "Note that effect estimates can vary over m_cde and c_cond values when interaction = TRUE.",
                                  fixed = TRUE)
                })
            })
            ##
            describe("coef.summary_regmedint", {
                it("extract the matrix object", {
                    expect_equal(coef(summary(fit_regmedint)),
                                 summary(fit_regmedint)$summary_myreg)
                    expect_equal(coef(summary(fit_regmedint_int)),
                                 summary(fit_regmedint_int)$summary_myreg)
                })
            })
        })
        ##
        describe("coef.regmedint", {
            it("creates a vector of estimates", {
                expect_equal(names(coef(fit_regmedint)),
                             c("cde","pnde","tnie","tnde","pnie","te","pm"))
            })
        })
        describe("vcov.regmedint", {
            it("creates a matrix with correct dimension names", {
                expect_equal(class(vcov(fit_regmedint))[[1]],
                             "matrix")
                expect_equal(colnames(vcov(fit_regmedint)),
                             c("cde","pnde","tnie","tnde","pnie","te","pm"))
                expect_equal(rownames(vcov(fit_regmedint)),
                             c("cde","pnde","tnie","tnde","pnie","te","pm"))
            })
        })
        ##
        describe("confint.regmedint", {
            it("creates a matrix of estimates", {
                expect_equal(colnames(confint(fit_regmedint)),
                             c("lower","upper"))
                expect_equal(rownames(confint(fit_regmedint)),
                             c("cde","pnde","tnie","tnde","pnie","te","pm"))
            })
            it("creates the lower column less than the upper column", {
                expect_true(all(confint(fit_regmedint)[,"lower"] <
                                confint(fit_regmedint)[,"upper"]))
            })

        })
    })
    ##
    describe("methods for regmedint mreg linear yreg logisitc", {
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
        fit_regmedint_int <- regmedint(data = pbc_cc,
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
                                       interaction = TRUE,
                                       casecontrol = FALSE,
                                       eventvar = NULL)
        ##
        describe("print.regmedint", {
            it("prints the mreg results", {
                expect_output(print(fit_regmedint),
                              deparse(fit_regmedint$mreg$call)[1],
                              fixed = TRUE)
            })
            it("prints the yreg results", {
                expect_output(print(fit_regmedint),
                              deparse(fit_regmedint$yreg$call)[1],
                              fixed = TRUE)
            })
            it("prints mediation analysis results with expected elements", {
                expect_output(print(fit_regmedint), "cde")
                expect_output(print(fit_regmedint), "pnde")
                expect_output(print(fit_regmedint), "tnie")
                expect_output(print(fit_regmedint), "tnde")
                expect_output(print(fit_regmedint), "pnie")
                expect_output(print(fit_regmedint), "te")
                expect_output(print(fit_regmedint), "pm")
            })
        })
        describe("summary.regmedint", {
            it("returns an object of summary_regmedint class", {
                expect_equal(class(summary(fit_regmedint))[[1]],
                             "summary_regmedint")
            })
            it("returns an object containing the mreg summary", {
                expect_equal(summary(fit_regmedint)$summary_mreg_fit,
                             summary(fit_regmedint$mreg_fit))
            })
            it("returns an object containing the yreg summary", {
                expect_equal(summary(fit_regmedint)$summary_yreg_fit,
                             summary(fit_regmedint$yreg_fit))
            })
            it("returns an object containing the myreg summary matrix", {
                expect_equal(class(summary(fit_regmedint)$summary_myreg)[[1]],
                             "matrix")
            })
            it("returns an object with appropriate columns", {
                expect_equal(colnames(summary(fit_regmedint)$summary_myreg),
                             c("est","se","Z","p","lower","upper"))
                ## expect_equivalent ignores attributes
                expect_equivalent(summary(fit_regmedint)$summary_myreg[,"est"],
                                  coef(fit_regmedint))
                expect_equivalent(summary(fit_regmedint)$summary_myreg[,"se"],
                                  sqrt(diag(vcov(fit_regmedint))))
                expect_equivalent(summary(fit_regmedint)$summary_myreg[,"lower"],
                                  confint(fit_regmedint)[,"lower"])
                expect_equivalent(summary(fit_regmedint)$summary_myreg[,"upper"],
                                  confint(fit_regmedint)[,"upper"])
            })
            it("returns an object with appropriate columns (exponentiated, exp(pm) is NA)", {
                expect_equal(colnames(summary(fit_regmedint, exponentiate = TRUE)$summary_myreg),
                             c("est","se","Z","p","lower","upper",
                               "exp(est)","exp(lower)","exp(upper)"))
                ## expect_equivalent ignores attributes
                expect_equivalent(summary(fit_regmedint, exponentiate = TRUE)$summary_myreg[,"est"],
                                  coef(fit_regmedint))
                expect_equivalent(summary(fit_regmedint, exponentiate = TRUE)$summary_myreg[,"se"],
                                  sqrt(diag(vcov(fit_regmedint))))
                expect_equivalent(summary(fit_regmedint, exponentiate = TRUE)$summary_myreg[,"lower"],
                                  confint(fit_regmedint)[,"lower"])
                expect_equivalent(summary(fit_regmedint, exponentiate = TRUE)$summary_myreg[,"upper"],
                                  confint(fit_regmedint)[,"upper"])
                exp_coefs <- exp(coef(fit_regmedint))
                exp_coefs["pm"] <- NA
                expect_equivalent(summary(fit_regmedint, exponentiate = TRUE)$summary_myreg[,"exp(est)"],
                                  exp_coefs)
                exp_lower <- exp(confint(fit_regmedint)[,"lower"])
                exp_lower["pm"] <- NA
                expect_equivalent(summary(fit_regmedint, exponentiate = TRUE)$summary_myreg[,"exp(lower)"],
                                  exp_lower)
                exp_upper <- exp(confint(fit_regmedint)[,"upper"])
                exp_upper["pm"] <- NA
                expect_equivalent(summary(fit_regmedint, exponentiate = TRUE)$summary_myreg[,"exp(upper)"],
                                  exp_upper)
            })
        })
        ##
        describe("methods for summary_regmedint", {
            describe("print.summary_regmedint", {
                ## Explicit printing within the function.
                ## No need to print the return value.
                it("prints the mreg results", {
                    expect_output(print(summary(fit_regmedint)),
                                  deparse(fit_regmedint$mreg$call)[1],
                                  fixed = TRUE)
                })
                it("prints the yreg results", {
                    expect_output(print(summary(fit_regmedint)),
                                  deparse(fit_regmedint$yreg$call)[1],
                                  fixed = TRUE)
                })
                it("prints mediation analysis results with expected elements", {
                    expect_output(print(summary(fit_regmedint)), "cde")
                    expect_output(print(summary(fit_regmedint)), "pnde")
                    expect_output(print(summary(fit_regmedint)), "tnie")
                    expect_output(print(summary(fit_regmedint)), "tnde")
                    expect_output(print(summary(fit_regmedint)), "pnie")
                    expect_output(print(summary(fit_regmedint)), "te")
                    expect_output(print(summary(fit_regmedint)), "pm")
                })
                it("prints mediation analysis results with expected elements exponentiated", {
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "cde")
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "pnde")
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "tnie")
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "tnde")
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "pnie")
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "te")
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "pm")
                })
                it("prints evaluation information", {
                    expect_output(print(summary(fit_regmedint)),
                                  "Evaluated at:",
                                  fixed = TRUE)
                    expect_output(print(summary(fit_regmedint)),
                                  "a1 (intervened value of avar) = ",
                                  fixed = TRUE)
                    expect_output(print(summary(fit_regmedint)),
                                  "a0 (reference value of avar)  = ",
                                  fixed = TRUE)
                    expect_output(print(summary(fit_regmedint)),
                                  "m_cde (intervend value of mvar for cde) = ",
                                  fixed = TRUE)
                    expect_output(print(summary(fit_regmedint)),
                                  "c_cond (covariate vector value) =",
                                  fixed = TRUE)
                    expect_output(print(summary(fit_regmedint)),
                                  "Note that effect estimates do not vary over m_cde and c_cond values when interaction = FALSE.",
                                  fixed = TRUE)
                    expect_output(print(summary(fit_regmedint_int)),
                                  "Note that effect estimates can vary over m_cde and c_cond values when interaction = TRUE.",
                                  fixed = TRUE)
                })
            })
            ##
            describe("coef.summary_regmedint", {
                it("extract the matrix object", {
                    expect_equal(coef(summary(fit_regmedint)),
                                 summary(fit_regmedint)$summary_myreg)
                    expect_equal(coef(summary(fit_regmedint_int)),
                                 summary(fit_regmedint_int)$summary_myreg)
                })
            })
        })
        ##
        describe("coef.regmedint", {
            it("creates a vector of estimates", {
                expect_equal(names(coef(fit_regmedint)),
                             c("cde","pnde","tnie","tnde","pnie","te","pm"))
            })
        })
        ##
        describe("vcov.regmedint", {
            it("creates a matrix with correct dimension names", {
                expect_equal(class(vcov(fit_regmedint))[[1]],
                             "matrix")
                expect_equal(colnames(vcov(fit_regmedint)),
                             c("cde","pnde","tnie","tnde","pnie","te","pm"))
                expect_equal(rownames(vcov(fit_regmedint)),
                             c("cde","pnde","tnie","tnde","pnie","te","pm"))
            })
        })
        ##
        describe("confint.regmedint", {
            it("creates a matrix of estimates", {
                expect_equal(colnames(confint(fit_regmedint)),
                             c("lower","upper"))
                expect_equal(rownames(confint(fit_regmedint)),
                             c("cde","pnde","tnie","tnde","pnie","te","pm"))
            })
            it("creates the lower column less than the upper column", {
                expect_true(all(confint(fit_regmedint)[,"lower"] <
                                confint(fit_regmedint)[,"upper"]))
            })

        })
    })
    ##
    describe("methods for regmedint mreg logistic yreg linear", {
        fit_regmedint <- regmedint(data = pbc_cc,
                                   yvar = "alk_phos",
                                   avar = "trt",
                                   mvar = "bili_bin",
                                   cvar = NULL,
                                   a0 = 1,
                                   a1 = 2,
                                   m_cde = 0,
                                   c_cond = NULL,
                                   mreg = "logistic",
                                   yreg = "linear",
                                   interaction = FALSE,
                                   casecontrol = FALSE,
                                   eventvar = NULL)
        fit_regmedint_int <- regmedint(data = pbc_cc,
                                       yvar = "alk_phos",
                                       avar = "trt",
                                       mvar = "bili_bin",
                                       cvar = NULL,
                                       a0 = 1,
                                       a1 = 2,
                                       m_cde = 0,
                                       c_cond = NULL,
                                       mreg = "logistic",
                                       yreg = "linear",
                                       interaction = TRUE,
                                       casecontrol = FALSE,
                                       eventvar = NULL)
        ##
        describe("print.regmedint", {
            it("prints the mreg results", {
                expect_output(print(fit_regmedint),
                              deparse(fit_regmedint$mreg$call)[1],
                              fixed = TRUE)
            })
            it("prints the yreg results", {
                expect_output(print(fit_regmedint),
                              deparse(fit_regmedint$yreg$call)[1],
                              fixed = TRUE)
            })
            it("prints mediation analysis results with expected elements", {
                expect_output(print(fit_regmedint), "cde")
                expect_output(print(fit_regmedint), "pnde")
                expect_output(print(fit_regmedint), "tnie")
                expect_output(print(fit_regmedint), "tnde")
                expect_output(print(fit_regmedint), "pnie")
                expect_output(print(fit_regmedint), "te")
                expect_output(print(fit_regmedint), "pm")
            })
        })
        describe("summary.regmedint", {
            it("returns an object of summary_regmedint class", {
                expect_equal(class(summary(fit_regmedint))[[1]],
                             "summary_regmedint")
            })
            it("returns an object containing the mreg summary", {
                expect_equal(summary(fit_regmedint)$summary_mreg_fit,
                             summary(fit_regmedint$mreg_fit))
            })
            it("returns an object containing the yreg summary", {
                expect_equal(summary(fit_regmedint)$summary_yreg_fit,
                             summary(fit_regmedint$yreg_fit))
            })
            it("returns an object containing the myreg summary matrix", {
                expect_equal(class(summary(fit_regmedint)$summary_myreg)[[1]],
                             "matrix")
            })
            it("returns an object with appropriate columns", {
                expect_equal(colnames(summary(fit_regmedint)$summary_myreg),
                             c("est","se","Z","p","lower","upper"))
                ## expect_equivalent ignores attributes
                expect_equivalent(summary(fit_regmedint)$summary_myreg[,"est"],
                                  coef(fit_regmedint))
                expect_equivalent(summary(fit_regmedint)$summary_myreg[,"se"],
                                  sqrt(diag(vcov(fit_regmedint))))
                expect_equivalent(summary(fit_regmedint)$summary_myreg[,"lower"],
                                  confint(fit_regmedint)[,"lower"])
                expect_equivalent(summary(fit_regmedint)$summary_myreg[,"upper"],
                                  confint(fit_regmedint)[,"upper"])
            })
            it("returns an object, ignoring exponentiate = TRUE", {
                expect_equal(colnames(summary(fit_regmedint, exponentiate = TRUE)$summary_myreg),
                             c("est","se","Z","p","lower","upper"))
                ## expect_equivalent ignores attributes
                expect_equivalent(summary(fit_regmedint, exponentiate = TRUE)$summary_myreg[,"est"],
                                  coef(fit_regmedint))
                expect_equivalent(summary(fit_regmedint, exponentiate = TRUE)$summary_myreg[,"se"],
                                  sqrt(diag(vcov(fit_regmedint))))
                expect_equivalent(summary(fit_regmedint, exponentiate = TRUE)$summary_myreg[,"lower"],
                                  confint(fit_regmedint)[,"lower"])
                expect_equivalent(summary(fit_regmedint, exponentiate = TRUE)$summary_myreg[,"upper"],
                                  confint(fit_regmedint)[,"upper"])
            })
        })
        ##
        describe("methods for summary_regmedint", {
            describe("print.summary_regmedint", {
                ## Explicit printing within the function.
                ## No need to print the return value.
                it("prints the mreg results", {
                    expect_output(print(summary(fit_regmedint)),
                                  deparse(fit_regmedint$mreg$call)[1],
                                  fixed = TRUE)
                })
                it("prints the yreg results", {
                    expect_output(print(summary(fit_regmedint)),
                                  deparse(fit_regmedint$yreg$call)[1],
                                  fixed = TRUE)
                })
                it("prints mediation analysis results with expected elements", {
                    expect_output(print(summary(fit_regmedint)), "cde")
                    expect_output(print(summary(fit_regmedint)), "pnde")
                    expect_output(print(summary(fit_regmedint)), "tnie")
                    expect_output(print(summary(fit_regmedint)), "tnde")
                    expect_output(print(summary(fit_regmedint)), "pnie")
                    expect_output(print(summary(fit_regmedint)), "te")
                    expect_output(print(summary(fit_regmedint)), "pm")
                })
                it("prints mediation analysis results with expected elements exponentiated", {
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "cde")
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "pnde")
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "tnie")
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "tnde")
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "pnie")
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "te")
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "pm")
                })
                it("prints evaluation information", {
                    expect_output(print(summary(fit_regmedint)),
                                  "Evaluated at:",
                                  fixed = TRUE)
                    expect_output(print(summary(fit_regmedint)),
                                  "a1 (intervened value of avar) = ",
                                  fixed = TRUE)
                    expect_output(print(summary(fit_regmedint)),
                                  "a0 (reference value of avar)  = ",
                                  fixed = TRUE)
                    expect_output(print(summary(fit_regmedint)),
                                  "m_cde (intervend value of mvar for cde) = ",
                                  fixed = TRUE)
                    expect_output(print(summary(fit_regmedint)),
                                  "c_cond (covariate vector value) =",
                                  fixed = TRUE)
                    expect_output(print(summary(fit_regmedint)),
                                  "Note that effect estimates do not vary over m_cde and c_cond values when interaction = FALSE.",
                                  fixed = TRUE)
                    expect_output(print(summary(fit_regmedint_int)),
                                  "Note that effect estimates can vary over m_cde and c_cond values when interaction = TRUE.",
                                  fixed = TRUE)
                })
            })
            ##
            describe("coef.summary_regmedint", {
                it("extract the matrix object", {
                    expect_equal(coef(summary(fit_regmedint)),
                                 summary(fit_regmedint)$summary_myreg)
                    expect_equal(coef(summary(fit_regmedint_int)),
                                 summary(fit_regmedint_int)$summary_myreg)
                })
            })
        })
        ##
        describe("coef.regmedint", {
            it("creates a vector of estimates", {
                expect_equal(names(coef(fit_regmedint)),
                             c("cde","pnde","tnie","tnde","pnie","te","pm"))
            })
        })
        ##
        describe("vcov.regmedint", {
            it("creates a matrix with correct dimension names", {
                expect_equal(class(vcov(fit_regmedint))[[1]],
                             "matrix")
                expect_equal(colnames(vcov(fit_regmedint)),
                             c("cde","pnde","tnie","tnde","pnie","te","pm"))
                expect_equal(rownames(vcov(fit_regmedint)),
                             c("cde","pnde","tnie","tnde","pnie","te","pm"))
            })
        })
        ##
        describe("confint.regmedint", {
            it("creates a matrix of estimates", {
                expect_equal(colnames(confint(fit_regmedint)),
                             c("lower","upper"))
                expect_equal(rownames(confint(fit_regmedint)),
                             c("cde","pnde","tnie","tnde","pnie","te","pm"))
            })
            it("creates the lower column less than the upper column", {
                expect_true(all(confint(fit_regmedint)[,"lower"] <
                                confint(fit_regmedint)[,"upper"]))
            })

        })
    })
    ##
    describe("methods for regmedint mreg logistic yreg logisitic", {
        fit_regmedint <- regmedint(data = pbc_cc,
                                   yvar = "spiders",
                                   avar = "trt",
                                   mvar = "bili_bin",
                                   cvar = NULL,
                                   a0 = 1,
                                   a1 = 2,
                                   m_cde = 0,
                                   c_cond = NULL,
                                   mreg = "logistic",
                                   yreg = "logistic",
                                   interaction = FALSE,
                                   casecontrol = FALSE,
                                   eventvar = NULL)
        fit_regmedint_int <- regmedint(data = pbc_cc,
                                       yvar = "spiders",
                                       avar = "trt",
                                       mvar = "bili_bin",
                                       cvar = NULL,
                                       a0 = 1,
                                       a1 = 2,
                                       m_cde = 0,
                                       c_cond = NULL,
                                       mreg = "logistic",
                                       yreg = "logistic",
                                       interaction = TRUE,
                                       casecontrol = FALSE,
                                       eventvar = NULL)
        ##
        describe("print.regmedint", {
            it("prints the mreg results", {
                expect_output(print(fit_regmedint),
                              deparse(fit_regmedint$mreg$call)[1],
                              fixed = TRUE)
            })
            it("prints the yreg results", {
                expect_output(print(fit_regmedint),
                              deparse(fit_regmedint$yreg$call)[1],
                              fixed = TRUE)
            })
            it("prints mediation analysis results with expected elements", {
                expect_output(print(fit_regmedint), "cde")
                expect_output(print(fit_regmedint), "pnde")
                expect_output(print(fit_regmedint), "tnie")
                expect_output(print(fit_regmedint), "tnde")
                expect_output(print(fit_regmedint), "pnie")
                expect_output(print(fit_regmedint), "te")
                expect_output(print(fit_regmedint), "pm")
            })
        })
        describe("summary.regmedint", {
            it("returns an object of summary_regmedint class", {
                expect_equal(class(summary(fit_regmedint))[[1]],
                             "summary_regmedint")
            })
            it("returns an object containing the mreg summary", {
                expect_equal(summary(fit_regmedint)$summary_mreg_fit,
                             summary(fit_regmedint$mreg_fit))
            })
            it("returns an object containing the yreg summary", {
                expect_equal(summary(fit_regmedint)$summary_yreg_fit,
                             summary(fit_regmedint$yreg_fit))
            })
            it("returns an object containing the myreg summary matrix", {
                expect_equal(class(summary(fit_regmedint)$summary_myreg)[[1]],
                             "matrix")
            })
            it("returns an object with appropriate columns", {
                expect_equal(colnames(summary(fit_regmedint)$summary_myreg),
                             c("est","se","Z","p","lower","upper"))
                ## expect_equivalent ignores attributes
                expect_equivalent(summary(fit_regmedint)$summary_myreg[,"est"],
                                  coef(fit_regmedint))
                expect_equivalent(summary(fit_regmedint)$summary_myreg[,"se"],
                                  sqrt(diag(vcov(fit_regmedint))))
                expect_equivalent(summary(fit_regmedint)$summary_myreg[,"lower"],
                                  confint(fit_regmedint)[,"lower"])
                expect_equivalent(summary(fit_regmedint)$summary_myreg[,"upper"],
                                  confint(fit_regmedint)[,"upper"])
            })
            it("returns an object with appropriate columns (exponentiated, exp(pm) is NA)", {
                expect_equal(colnames(summary(fit_regmedint, exponentiate = TRUE)$summary_myreg),
                             c("est","se","Z","p","lower","upper",
                               "exp(est)","exp(lower)","exp(upper)"))
                ## expect_equivalent ignores attributes
                expect_equivalent(summary(fit_regmedint, exponentiate = TRUE)$summary_myreg[,"est"],
                                  coef(fit_regmedint))
                expect_equivalent(summary(fit_regmedint, exponentiate = TRUE)$summary_myreg[,"se"],
                                  sqrt(diag(vcov(fit_regmedint))))
                expect_equivalent(summary(fit_regmedint, exponentiate = TRUE)$summary_myreg[,"lower"],
                                  confint(fit_regmedint)[,"lower"])
                expect_equivalent(summary(fit_regmedint, exponentiate = TRUE)$summary_myreg[,"upper"],
                                  confint(fit_regmedint)[,"upper"])
                exp_coefs <- exp(coef(fit_regmedint))
                exp_coefs["pm"] <- NA
                expect_equivalent(summary(fit_regmedint, exponentiate = TRUE)$summary_myreg[,"exp(est)"],
                                  exp_coefs)
                exp_lower <- exp(confint(fit_regmedint)[,"lower"])
                exp_lower["pm"] <- NA
                expect_equivalent(summary(fit_regmedint, exponentiate = TRUE)$summary_myreg[,"exp(lower)"],
                                  exp_lower)
                exp_upper <- exp(confint(fit_regmedint)[,"upper"])
                exp_upper["pm"] <- NA
                expect_equivalent(summary(fit_regmedint, exponentiate = TRUE)$summary_myreg[,"exp(upper)"],
                                  exp_upper)
            })
        })
        ##
        describe("methods for summary_regmedint", {
            describe("print.summary_regmedint", {
                ## Explicit printing within the function.
                ## No need to print the return value.
                it("prints the mreg results", {
                    expect_output(print(summary(fit_regmedint)),
                                  deparse(fit_regmedint$mreg$call)[1],
                                  fixed = TRUE)
                })
                it("prints the yreg results", {
                    expect_output(print(summary(fit_regmedint)),
                                  deparse(fit_regmedint$yreg$call)[1],
                                  fixed = TRUE)
                })
                it("prints mediation analysis results with expected elements", {
                    expect_output(print(summary(fit_regmedint)), "cde")
                    expect_output(print(summary(fit_regmedint)), "pnde")
                    expect_output(print(summary(fit_regmedint)), "tnie")
                    expect_output(print(summary(fit_regmedint)), "tnde")
                    expect_output(print(summary(fit_regmedint)), "pnie")
                    expect_output(print(summary(fit_regmedint)), "te")
                    expect_output(print(summary(fit_regmedint)), "pm")
                })
                it("prints mediation analysis results with expected elements exponentiated", {
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "cde")
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "pnde")
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "tnie")
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "tnde")
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "pnie")
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "te")
                    expect_output(print(summary(fit_regmedint, exponentiate = TRUE)), "pm")
                })
                it("prints evaluation information", {
                    expect_output(print(summary(fit_regmedint)),
                                  "Evaluated at:",
                                  fixed = TRUE)
                    expect_output(print(summary(fit_regmedint)),
                                  "a1 (intervened value of avar) = ",
                                  fixed = TRUE)
                    expect_output(print(summary(fit_regmedint)),
                                  "a0 (reference value of avar)  = ",
                                  fixed = TRUE)
                    expect_output(print(summary(fit_regmedint)),
                                  "m_cde (intervend value of mvar for cde) = ",
                                  fixed = TRUE)
                    expect_output(print(summary(fit_regmedint)),
                                  "c_cond (covariate vector value) =",
                                  fixed = TRUE)
                    expect_output(print(summary(fit_regmedint)),
                                  "Note that effect estimates do not vary over m_cde and c_cond values when interaction = FALSE.",
                                  fixed = TRUE)
                    expect_output(print(summary(fit_regmedint_int)),
                                  "Note that effect estimates can vary over m_cde and c_cond values when interaction = TRUE.",
                                  fixed = TRUE)
                })
            })
            ##
            describe("coef.summary_regmedint", {
                it("extract the matrix object", {
                    expect_equal(coef(summary(fit_regmedint)),
                                 summary(fit_regmedint)$summary_myreg)
                    expect_equal(coef(summary(fit_regmedint_int)),
                                 summary(fit_regmedint_int)$summary_myreg)
                })
            })
        })
        ##
        describe("coef.regmedint", {
            it("creates a vector of estimates", {
                expect_equal(names(coef(fit_regmedint)),
                             c("cde","pnde","tnie","tnde","pnie","te","pm"))
            })
        })
        ##
        describe("vcov.regmedint", {
            it("creates a matrix with correct dimension names", {
                expect_equal(class(vcov(fit_regmedint))[[1]],
                             "matrix")
                expect_equal(colnames(vcov(fit_regmedint)),
                             c("cde","pnde","tnie","tnde","pnie","te","pm"))
                expect_equal(rownames(vcov(fit_regmedint)),
                             c("cde","pnde","tnie","tnde","pnie","te","pm"))
            })
        })
        ##
        describe("confint.regmedint", {
            it("creates a matrix of estimates", {
                expect_equal(colnames(confint(fit_regmedint)),
                             c("lower","upper"))
                expect_equal(rownames(confint(fit_regmedint)),
                             c("cde","pnde","tnie","tnde","pnie","te","pm"))
            })
            it("creates the lower column less than the upper column", {
                expect_true(all(confint(fit_regmedint)[,"lower"] <
                                confint(fit_regmedint)[,"upper"]))
            })

        })
    })

})
