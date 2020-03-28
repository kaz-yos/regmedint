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
        ##
        describe("summary.regmedint", {
            ## Explicit printing within the function.
            ## No need to print the return value.
            it("prints the mreg results", {
                expect_output(summary(fit_regmedint),
                              deparse(fit_regmedint$mreg$call)[1],
                              fixed = TRUE)
            })
            it("prints the yreg results", {
                expect_output(summary(fit_regmedint),
                              deparse(fit_regmedint$yreg$call)[1],
                              fixed = TRUE)
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
            it("prints evaluation information", {
                expect_output(summary(fit_regmedint),
                              "Evaluated at:")
                expect_output(summary(fit_regmedint),
                              "a1 (intervened value of avar) = ")
                expect_output(summary(fit_regmedint),
                              "a0 (reference value of avar)  = ")
                expect_output(summary(fit_regmedint),
                              "m_cde (intervend value of mvar for cde) = ")
                expect_output(summary(fit_regmedint),
                              "c_cond (covariate vector value) =")
                expect_output(summary(fit_regmedint),
                              "Note that effect estimates do not vary over m_cde and c_cond values when interaction = FALSE.")
                expect_output(summary(fit_regmedint_int),
                              "Note that effect estimates can vary over m_cde and c_cond values when interaction = TRUE.")
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
            it("creates the lower column less than the upper column", {
                expect_true(all(confint(fit_regmedint)[,"lower"] <
                                confint(fit_regmedint)[,"upper"]))
            })

        })
    })

})
