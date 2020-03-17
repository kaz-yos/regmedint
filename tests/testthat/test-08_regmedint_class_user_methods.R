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

    fit_mreg_lin_yreg_lin <- regmedint(data = pbc_cc,
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
            expect_output(print(fit_mreg_lin_yreg_lin), "cde")
            expect_output(print(fit_mreg_lin_yreg_lin), "pnde")
            expect_output(print(fit_mreg_lin_yreg_lin), "tnie")
            expect_output(print(fit_mreg_lin_yreg_lin), "tnde")
            expect_output(print(fit_mreg_lin_yreg_lin), "pnie")
            expect_output(print(fit_mreg_lin_yreg_lin), "te")
            expect_output(print(fit_mreg_lin_yreg_lin), "pm")
        })
    })
    ##
    describe("summary.regmedint", {
            expect_output(summary(fit_mreg_lin_yreg_lin), "cde")
            expect_output(summary(fit_mreg_lin_yreg_lin), "pnde")
            expect_output(summary(fit_mreg_lin_yreg_lin), "tnie")
            expect_output(summary(fit_mreg_lin_yreg_lin), "tnde")
            expect_output(summary(fit_mreg_lin_yreg_lin), "pnie")
            expect_output(summary(fit_mreg_lin_yreg_lin), "te")
            expect_output(summary(fit_mreg_lin_yreg_lin), "pm")
    })
    ##
    describe("coef.regmedint", {
            expect_output(coef(fit_mreg_lin_yreg_lin), "cde")
            expect_output(coef(fit_mreg_lin_yreg_lin), "pnde")
            expect_output(coef(fit_mreg_lin_yreg_lin), "tnie")
            expect_output(coef(fit_mreg_lin_yreg_lin), "tnde")
            expect_output(coef(fit_mreg_lin_yreg_lin), "pnie")
            expect_output(coef(fit_mreg_lin_yreg_lin), "te")
            expect_output(coef(fit_mreg_lin_yreg_lin), "pm")
    })
    ##
    describe("confint.regmedint", {
            expect_output(confint(fit_mreg_lin_yreg_lin), "cde")
            expect_output(confint(fit_mreg_lin_yreg_lin), "pnde")
            expect_output(confint(fit_mreg_lin_yreg_lin), "tnie")
            expect_output(confint(fit_mreg_lin_yreg_lin), "tnde")
            expect_output(confint(fit_mreg_lin_yreg_lin), "pnie")
            expect_output(confint(fit_mreg_lin_yreg_lin), "te")
            expect_output(confint(fit_mreg_lin_yreg_lin), "pm")
    })
})
