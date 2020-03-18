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
            it("prints results with expected elements", {
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
            expect_output(summary(fit_regmedint), "cde")
            expect_output(summary(fit_regmedint), "pnde")
            expect_output(summary(fit_regmedint), "tnie")
            expect_output(summary(fit_regmedint), "tnde")
            expect_output(summary(fit_regmedint), "pnie")
            expect_output(summary(fit_regmedint), "te")
            expect_output(summary(fit_regmedint), "pm")
        })
        ##
        describe("coef.regmedint", {
            expect_output(coef(fit_regmedint), "cde")
            expect_output(coef(fit_regmedint), "pnde")
            expect_output(coef(fit_regmedint), "tnie")
            expect_output(coef(fit_regmedint), "tnde")
            expect_output(coef(fit_regmedint), "pnie")
            expect_output(coef(fit_regmedint), "te")
            expect_output(coef(fit_regmedint), "pm")
        })
        ##
        describe("confint.regmedint", {
            expect_output(confint(fit_regmedint), "cde")
            expect_output(confint(fit_regmedint), "pnde")
            expect_output(confint(fit_regmedint), "tnie")
            expect_output(confint(fit_regmedint), "tnde")
            expect_output(confint(fit_regmedint), "pnie")
            expect_output(confint(fit_regmedint), "te")
            expect_output(confint(fit_regmedint), "pm")
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
            it("prints results with expected elements", {
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
            expect_output(summary(fit_regmedint), "cde")
            expect_output(summary(fit_regmedint), "pnde")
            expect_output(summary(fit_regmedint), "tnie")
            expect_output(summary(fit_regmedint), "tnde")
            expect_output(summary(fit_regmedint), "pnie")
            expect_output(summary(fit_regmedint), "te")
            expect_output(summary(fit_regmedint), "pm")
        })
        ##
        describe("coef.regmedint", {
            expect_output(coef(fit_regmedint), "cde")
            expect_output(coef(fit_regmedint), "pnde")
            expect_output(coef(fit_regmedint), "tnie")
            expect_output(coef(fit_regmedint), "tnde")
            expect_output(coef(fit_regmedint), "pnie")
            expect_output(coef(fit_regmedint), "te")
            expect_output(coef(fit_regmedint), "pm")
        })
        ##
        describe("confint.regmedint", {
            expect_output(confint(fit_regmedint), "cde")
            expect_output(confint(fit_regmedint), "pnde")
            expect_output(confint(fit_regmedint), "tnie")
            expect_output(confint(fit_regmedint), "tnde")
            expect_output(confint(fit_regmedint), "pnie")
            expect_output(confint(fit_regmedint), "te")
            expect_output(confint(fit_regmedint), "pm")
        })
    })
})
