pkgname <- "regmedint"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "regmedint-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('regmedint')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("coef.regmedint")
### * coef.regmedint

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: coef.regmedint
### Title: Extract point estimates.
### Aliases: coef.regmedint

### ** Examples

library(regmedint)
data(vv2015)
regmedint_obj <- regmedint(data = vv2015,
                           ## Variables
                           yvar = "y",
                           avar = "x",
                           mvar = "m",
                           cvar = c("c"),
                           eventvar = "event",
                           ## Values at which effects are evaluated
                           a0 = 0,
                           a1 = 1,
                           m_cde = 1,
                           c_cond = 0.5,
                           ## Model types
                           mreg = "logistic",
                           yreg = "survAFT_weibull",
                           ## Additional specification
                           interaction = TRUE,
                           casecontrol = FALSE)
coef(regmedint_obj)
## Evaluate at different values
coef(regmedint_obj, m_cde = 0, c_cond = 1)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("coef.regmedint", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("coef.summary_regmedint")
### * coef.summary_regmedint

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: coef.summary_regmedint
### Title: Extract the result matrix from a summary_regmedint object.
### Aliases: coef.summary_regmedint

### ** Examples

library(regmedint)
data(vv2015)
regmedint_obj <- regmedint(data = vv2015,
                           ## Variables
                           yvar = "y",
                           avar = "x",
                           mvar = "m",
                           cvar = c("c"),
                           eventvar = "event",
                           ## Values at which effects are evaluated
                           a0 = 0,
                           a1 = 1,
                           m_cde = 1,
                           c_cond = 0.5,
                           ## Model types
                           mreg = "logistic",
                           yreg = "survAFT_weibull",
                           ## Additional specification
                           interaction = TRUE,
                           casecontrol = FALSE)
coef(summary(regmedint_obj))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("coef.summary_regmedint", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("confint.regmedint")
### * confint.regmedint

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: confint.regmedint
### Title: Confidence intervals for mediation prameter estimates.
### Aliases: confint.regmedint

### ** Examples

library(regmedint)
data(vv2015)
regmedint_obj <- regmedint(data = vv2015,
                           ## Variables
                           yvar = "y",
                           avar = "x",
                           mvar = "m",
                           cvar = c("c"),
                           eventvar = "event",
                           ## Values at which effects are evaluated
                           a0 = 0,
                           a1 = 1,
                           m_cde = 1,
                           c_cond = 0.5,
                           ## Model types
                           mreg = "logistic",
                           yreg = "survAFT_weibull",
                           ## Additional specification
                           interaction = TRUE,
                           casecontrol = FALSE)
confint(regmedint_obj)
## Evaluate at different values
confint(regmedint_obj, m_cde = 0, c_cond = 1)
## Change confidence level
confint(regmedint_obj, m_cde = 0, c_cond = 1, level = 0.99)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("confint.regmedint", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print.regmedint")
### * print.regmedint

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print.regmedint
### Title: print method for regmedint object
### Aliases: print.regmedint

### ** Examples

library(regmedint)
data(vv2015)
regmedint_obj <- regmedint(data = vv2015,
                           ## Variables
                           yvar = "y",
                           avar = "x",
                           mvar = "m",
                           cvar = c("c"),
                           eventvar = "event",
                           ## Values at which effects are evaluated
                           a0 = 0,
                           a1 = 1,
                           m_cde = 1,
                           c_cond = 0.5,
                           ## Model types
                           mreg = "logistic",
                           yreg = "survAFT_weibull",
                           ## Additional specification
                           interaction = TRUE,
                           casecontrol = FALSE)
## Implicit printing
regmedint_obj
## Explicit printing
print(regmedint_obj)
## Evaluate at different values
print(regmedint_obj, m_cde = 0, c_cond = 1)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print.regmedint", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print.summary_regmedint")
### * print.summary_regmedint

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print.summary_regmedint
### Title: Print method for summary objects from 'summary.regmedint'
### Aliases: print.summary_regmedint

### ** Examples

library(regmedint)
data(vv2015)
regmedint_obj <- regmedint(data = vv2015,
                           ## Variables
                           yvar = "y",
                           avar = "x",
                           mvar = "m",
                           cvar = c("c"),
                           eventvar = "event",
                           ## Values at which effects are evaluated
                           a0 = 0,
                           a1 = 1,
                           m_cde = 1,
                           c_cond = 0.5,
                           ## Model types
                           mreg = "logistic",
                           yreg = "survAFT_weibull",
                           ## Additional specification
                           interaction = TRUE,
                           casecontrol = FALSE)
## Implicit printing
summary(regmedint_obj)
## Explicit printing
print(summary(regmedint_obj))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print.summary_regmedint", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("regmedint")
### * regmedint

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: regmedint
### Title: regmedint: A package for regression-based causal mediation
###   analysis
### Aliases: regmedint

### ** Examples

library(regmedint)
data(vv2015)
regmedint_obj <- regmedint(data = vv2015,
                           ## Variables
                           yvar = "y",
                           avar = "x",
                           mvar = "m",
                           cvar = c("c"),
                           eventvar = "event",
                           ## Values at which effects are evaluated
                           a0 = 0,
                           a1 = 1,
                           m_cde = 1,
                           c_cond = 0.5,
                           ## Model types
                           mreg = "logistic",
                           yreg = "survAFT_weibull",
                           ## Additional specification
                           interaction = TRUE,
                           casecontrol = FALSE)
summary(regmedint_obj)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("regmedint", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary.regmedint")
### * summary.regmedint

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary.regmedint
### Title: summary method for regmedint object
### Aliases: summary.regmedint

### ** Examples

library(regmedint)
data(vv2015)
regmedint_obj <- regmedint(data = vv2015,
                           ## Variables
                           yvar = "y",
                           avar = "x",
                           mvar = "m",
                           cvar = c("c"),
                           eventvar = "event",
                           ## Values at which effects are evaluated
                           a0 = 0,
                           a1 = 1,
                           m_cde = 1,
                           c_cond = 0.5,
                           ## Model types
                           mreg = "logistic",
                           yreg = "survAFT_weibull",
                           ## Additional specification
                           interaction = TRUE,
                           casecontrol = FALSE)
## Detailed result with summary
summary(regmedint_obj)
## Add exponentiate results for non-linear outcome models
summary(regmedint_obj, exponentiate = TRUE)
## Evaluate at different values
summary(regmedint_obj, m_cde = 0, c_cond = 1)
## Change confidence level
summary(regmedint_obj, m_cde = 0, c_cond = 1, level = 0.99)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary.regmedint", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("vcov.regmedint")
### * vcov.regmedint

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: vcov.regmedint
### Title: Extract variance estimates in the vcov form.
### Aliases: vcov.regmedint

### ** Examples

library(regmedint)
data(vv2015)
regmedint_obj <- regmedint(data = vv2015,
                           ## Variables
                           yvar = "y",
                           avar = "x",
                           mvar = "m",
                           cvar = c("c"),
                           eventvar = "event",
                           ## Values at which effects are evaluated
                           a0 = 0,
                           a1 = 1,
                           m_cde = 1,
                           c_cond = 0.5,
                           ## Model types
                           mreg = "logistic",
                           yreg = "survAFT_weibull",
                           ## Additional specification
                           interaction = TRUE,
                           casecontrol = FALSE)
vcov(regmedint_obj)
## Evaluate at different values
vcov(regmedint_obj, m_cde = 0, c_cond = 1)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("vcov.regmedint", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
