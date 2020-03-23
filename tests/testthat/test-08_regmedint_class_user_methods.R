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
               status = if_else(status == 0, 0L, 1L))

    describe("methods for regmedint mreg linear yreg linear", {
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
        ##
        describe("print.regmedint", {
            it("prints the mreg results", {
                expect_output(print(fit_regmedint),
                              deparse(fit_regmedint$mreg$call)[1])
            })
            it("prints the yreg results", {
                expect_output(print(fit_regmedint),
                              deparse(fit_regmedint$yreg$call)[1])
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
        ##
        describe("summary.regmedint", {
            it("prints the mreg results", {
                expect_output(print(summary(fit_regmedint)),
                              deparse(fit_regmedint$mreg$call)[1])
            })
            it("prints the yreg results", {
                expect_output(print(summary(fit_regmedint)),
                              deparse(fit_regmedint$yreg$call)[1])
            })
            it("prints mediation analysis results with expected elements", {
                expect_output(summary(fit_regmedint), "cde")
                expect_output(summary(fit_regmedint), "pnde")
                expect_output(summary(fit_regmedint), "tnie")
                expect_output(summary(fit_regmedint), "tnde")
                expect_output(summary(fit_regmedint), "pnie")
                expect_output(summary(fit_regmedint), "te")
                expect_output(summary(fit_regmedint), "pm")
            })
            it("prints mediation analysis results with expected elements exponentiated", {
                expect_output(summary(fit_regmedint, exponentiate = TRUE), "cde")
                expect_output(summary(fit_regmedint, exponentiate = TRUE), "pnde")
                expect_output(summary(fit_regmedint, exponentiate = TRUE), "tnie")
                expect_output(summary(fit_regmedint, exponentiate = TRUE), "tnde")
                expect_output(summary(fit_regmedint, exponentiate = TRUE), "pnie")
                expect_output(summary(fit_regmedint, exponentiate = TRUE), "te")
                expect_output(summary(fit_regmedint, exponentiate = TRUE), "pm")
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
        describe("confint.regmedint", {
            it("creates a matrix of estimates", {
                expect_equal(colnames(confint(fit_regmedint)),
                             c("lower","upper"))
                expect_equal(rownames(confint(fit_regmedint)),
                             c("cde","pnde","tnie","tnde","pnie","te","pm"))
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
        ##
                describe("print.regmedint", {
            it("prints the mreg results", {
                expect_output(print(fit_regmedint),
                              deparse(fit_regmedint$mreg$call)[1])
            })
            it("prints the yreg results", {
                expect_output(print(fit_regmedint),
                              deparse(fit_regmedint$yreg$call)[1])
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
        ##
        describe("summary.regmedint", {
            it("prints the mreg results", {
                expect_output(print(summary(fit_regmedint)),
                              deparse(fit_regmedint$mreg$call)[1])
            })
            it("prints the yreg results", {
                expect_output(print(summary(fit_regmedint)),
                              deparse(fit_regmedint$yreg$call)[1])
            })
            it("prints mediation analysis results with expected elements", {
                expect_output(summary(fit_regmedint), "cde")
                expect_output(summary(fit_regmedint), "pnde")
                expect_output(summary(fit_regmedint), "tnie")
                expect_output(summary(fit_regmedint), "tnde")
                expect_output(summary(fit_regmedint), "pnie")
                expect_output(summary(fit_regmedint), "te")
                expect_output(summary(fit_regmedint), "pm")
            })
            it("prints mediation analysis results with expected elements exponentiated", {
                expect_output(summary(fit_regmedint, exponentiate = TRUE), "cde")
                expect_output(summary(fit_regmedint, exponentiate = TRUE), "pnde")
                expect_output(summary(fit_regmedint, exponentiate = TRUE), "tnie")
                expect_output(summary(fit_regmedint, exponentiate = TRUE), "tnde")
                expect_output(summary(fit_regmedint, exponentiate = TRUE), "pnie")
                expect_output(summary(fit_regmedint, exponentiate = TRUE), "te")
                expect_output(summary(fit_regmedint, exponentiate = TRUE), "pm")
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
        describe("confint.regmedint", {
            it("creates a matrix of estimates", {
                expect_equal(colnames(confint(fit_regmedint)),
                             c("lower","upper"))
                expect_equal(rownames(confint(fit_regmedint)),
                             c("cde","pnde","tnie","tnde","pnie","te","pm"))
            })
        })
    })

})
