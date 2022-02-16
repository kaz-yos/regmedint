################################################################################
### User interface functions
##
## Created on: 2020-03-09
## Author: Kazuki Yoshida
################################################################################


###
### Main user interface function
################################################################################

##
## https://adv-r.hadley.nz/s3.html#s3-classes
## Follow the recommendation for S3 classes in Advanced R.
## Provide a low-level constructor, validator, and user-friendly helper.
##
## https://adv-r.hadley.nz/s3.html#s3-constructor
## https://adv-r.hadley.nz/s3.html#validators
## https://adv-r.hadley.nz/s3.html#helpers


##' Conduct regression-based causal mediation analysis
##'
##' This is a user-interface for regression-based causal mediation analysis as described in Valeri & VanderWeele 2013 and Valeri & VanderWeele 2015.
##'
##' @param data Data frame containing the following relevant variables.
##' @param yvar A character vector of length 1. Outcome variable name. It should be the time variable for the survival outcome.
##' @param avar A character vector of length 1. Treatment variable name.
##' @param mvar A character vector of length 1. Mediator variable name.
##' @param cvar A character vector of length > 0. Covariate names. Use \code{NULL} if there is no covariate. However, this is a highly suspicious situation. Even if \code{avar} is randomized, \code{mvar} is not. Thus, there are usually some confounder(s) to account for the common cause structure (confounding) between \code{mvar} and \code{yvar}. 
##' @param emm_ac_mreg A character vector of length > 0. Effect modifiers names. The covariate vector in treatment-covariate product term in the mediator model.
##' @param emm_ac_yreg A character vector of length > 0. Effect modifiers names. The covariate vector in treatment-covariate product term in the outcome model. 
##' @param emm_mc_yreg A character vector of length > 0. Effect modifiers names. The covariate vector in mediator-covariate product term in outcome model. 
##' @param eventvar An character vector of length 1. Only required for survival outcome regression models. Note that the coding is 1 for event and 0 for censoring, following the R survival package convention.
##' @param a0 A numeric vector of length 1. The reference level of treatment variable that is considered "untreated" or "unexposed".
##' @param a1 A numeric vector of length 1.
##' @param m_cde A numeric vector of length 1. Mediator level at which controlled direct effect is evaluated at.
##' @param c_cond A numeric vector of the same length as \code{cvar}. Covariate levels at which natural direct and indirect effects are evaluated at. 
##' @param mreg A character vector of length 1. Mediator regression type: \code{"linear"} or \code{"logistic"}.
##' @param yreg A character vector of length 1. Outcome regression type: \code{"linear"}, \code{"logistic"}, \code{"loglinear"}, \code{"poisson"}, \code{"negbin"}, \code{"survCox"}, \code{"survAFT_exp"}, or \code{"survAFT_weibull"}.
##' @param interaction A logical vector of length 1. The presence of treatment-mediator interaction in the outcome model. Default to TRUE.
##' @param casecontrol A logical vector of length 1. Default to FALSE. Whether data comes from a case-control study.
##' @param na_omit A logical vector of length 1. Default to FALSE. Whether to remove NAs in the columns of interest before fitting the models.
##' 
##' @return regmedint object, which is a list containing the mediator regression object, the outcome regression object, and the regression-based mediation results.
##'
##' @examples
##' library(regmedint)
##' data(vv2015)
##' regmedint_obj1 <- regmedint(data = vv2015,
##'                             ## Variables
##'                             yvar = "y",
##'                             avar = "x",
##'                             mvar = "m",
##'                             cvar = c("c"),
##'                             eventvar = "event",
##'                             ## Values at which effects are evaluated
##'                             a0 = 0,
##'                             a1 = 1,
##'                             m_cde = 1,
##'                             c_cond = 3,
##'                             ## Model types
##'                             mreg = "logistic",
##'                             yreg = "survAFT_weibull",
##'                             ## Additional specification
##'                             interaction = TRUE,
##'                             casecontrol = FALSE)
##' summary(regmedint_obj1)
##' 
##' regmedint_obj2 <- regmedint(data = vv2015,
##'                             ## Variables
##'                             yvar = "y",
##'                             avar = "x",
##'                             mvar = "m",
##'                             cvar = c("c"),
##'                             emm_ac_mreg = c("c"), 
##'                             emm_ac_yreg = c("c"), 
##'                             emm_mc_yreg = c("c"), 
##'                             eventvar = "event",
##'                             ## Values at which effects are evaluated
##'                             a0 = 0,
##'                             a1 = 1,
##'                             m_cde = 1,
##'                             c_cond = 3,
##'                             ## Model types
##'                             mreg = "logistic",
##'                             yreg = "survAFT_weibull",
##'                             ## Additional specification
##'                             interaction = TRUE,
##'                             casecontrol = FALSE)
##' summary(regmedint_obj2)
##' 
##' 
##' 
##'
##' @export
regmedint <- function(data,
                      yvar,
                      avar,
                      mvar,
                      cvar,
                      emm_ac_mreg = NULL, 
                      emm_ac_yreg = NULL, 
                      emm_mc_yreg = NULL, 
                      eventvar = NULL,
                      a0,
                      a1,
                      m_cde,
                      c_cond,
                      mreg,
                      yreg,
                      interaction = TRUE,
                      casecontrol = FALSE,
                      na_omit = FALSE) {
    
    ## This is the user-friendly helper function with a name that is the class name.
    ## https://adv-r.hadley.nz/s3.html#helpers
    
    ## Handle missing value
    ## Select columns of interest only.
    data <- data[,c(yvar, avar, mvar, cvar, eventvar)]
    ## Report NA
    report_missing(data, yvar, avar, mvar, cvar, eventvar)
    ## Construct a complete case dataset if requested via na_omit
    if(any(is.na(data)) && na_omit) {
        data <- na.omit(data)
    }
    
    ## Check data contains yvar, avar, mvar, cvar, eventvar (if provided), emm_ac_mreg, emm_ac_yreg, emm_mc_yreg
    validate_args(data = data,
                  yvar = yvar,
                  avar = avar,
                  mvar = mvar,
                  cvar = cvar,
                  emm_ac_mreg = emm_ac_mreg, 
                  emm_ac_yreg = emm_ac_yreg, 
                  emm_mc_yreg = emm_mc_yreg, 
                  a0 = a0,
                  a1 = a1,
                  m_cde = m_cde,
                  c_cond = c_cond,
                  mreg = mreg,
                  yreg = yreg,
                  interaction = interaction,
                  casecontrol = casecontrol,
                  eventvar = eventvar)
    
    ## Construct a regmedint object after argument validation.
    ## This is the low-level constructor function.
    ## https://adv-r.hadley.nz/s3.html#s3-constructor
    res <- new_regmedint(data = data,
                         yvar = yvar,
                         avar = avar,
                         mvar = mvar,
                         cvar = cvar,
                         emm_ac_mreg = emm_ac_mreg, 
                         emm_ac_yreg = emm_ac_yreg, 
                         emm_mc_yreg = emm_mc_yreg, 
                         a0 = a0,
                         a1 = a1,
                         m_cde = m_cde,
                         c_cond = c_cond,
                         mreg = mreg,
                         yreg = yreg,
                         interaction = interaction,
                         casecontrol = casecontrol,
                         eventvar = eventvar)
    
    ## Check the resulting object for anomalies.
    ## This is the class validator function.
    ## https://adv-r.hadley.nz/s3.html#validators
    validate_regmedint(res)
    
    ## Return results
    res
}




###
#### Functions to handle missing data
################################################################################

##' Report variables with missing data
##'
##' Report the number of missing observations for each variables of interest relevant for the analysis
##'
##' @inheritParams regmedint
##'
##' @return No return value, called for side effects.
report_missing <- function(data, yvar, avar, mvar, cvar, eventvar){
    ## Only report missing in these variables.
    data <- data[c(yvar, avar, mvar, cvar, eventvar)]
    ## Add warning: NAs
    ## If the dataset contains NAs, print a general warning message
    if (any(is.na(data))) {
        message("Dataset contains NAs.")
    }
    ## Print the number of NAs by columns of interest
    for (i in 1:ncol(data)) {
        if (any(is.na(data[,i]))) {
            message(paste(colnames(data)[i],
                          "has",
                          sum(is.na(data[,i])),
                          "NAs."))
        }
    }
}





###
### Argument validation function
################################################################################

##' Validate arguments to regmedint before passing to other functions
##'
##' Internal functions (usually) do not validate arguments, thus, we need to make sure informative errors are raised when the arguments are not safe for subsequent computation.
##'
##' @inheritParams regmedint
##'
##' @return No return value, called for side effects.
validate_args <- function(data,
                          yvar,
                          avar,
                          mvar,
                          cvar,
                          emm_ac_mreg, 
                          emm_ac_yreg, 
                          emm_mc_yreg, 
                          eventvar,
                          a0,
                          a1,
                          m_cde,
                          c_cond,
                          mreg,
                          yreg,
                          interaction,
                          casecontrol) {
    
    ## Dataset
    assertthat::assert_that(is.data.frame(data))
    
    ##
    assertthat::assert_that(is.character(yvar))
    assertthat::assert_that(length(yvar) == 1)
    ##
    assertthat::assert_that(is.character(avar))
    assertthat::assert_that(length(avar) == 1)
    ##
    assertthat::assert_that(is.character(mvar))
    assertthat::assert_that(length(mvar) == 1)
    ##
    assertthat::assert_that(is.null(cvar) | is.character(cvar))
    ##
    assertthat::assert_that(is.null(emm_ac_mreg) | is.character(emm_ac_mreg))
    ##
    assertthat::assert_that(is.null(emm_ac_yreg) | is.character(emm_ac_yreg))
    ##
    assertthat::assert_that(is.null(emm_mc_yreg) | is.character(emm_mc_yreg))
    ##
    assertthat::assert_that(is.numeric(a0))
    assertthat::assert_that(length(a0) == 1)
    assertthat::assert_that(!is.na(a0))
    ##
    assertthat::assert_that(is.numeric(a1))
    assertthat::assert_that(length(a1) == 1)
    assertthat::assert_that(!is.na(a1))
    ##
    assertthat::assert_that(is.numeric(m_cde))
    assertthat::assert_that(length(m_cde) == 1)
    assertthat::assert_that(!is.na(m_cde))
    ##
    assertthat::assert_that(is.null(c_cond) | is.numeric(c_cond))
    assertthat::assert_that(length(c_cond) == length(cvar))
    assertthat::assert_that(all(!is.na(c_cond)))
    ##
    assertthat::assert_that(is.character(mreg))
    assertthat::assert_that(length(mreg) == 1)
    assertthat::assert_that(mreg %in% c("linear","logistic"))
    ##
    assertthat::assert_that(is.character(yreg))
    assertthat::assert_that(length(yreg) == 1)
    assertthat::assert_that(yreg %in% c("linear",
                                        "logistic","loglinear","poisson","negbin",
                                        "survCox","survAFT_exp","survAFT_weibull"))
    ##
    assertthat::assert_that(is.logical(interaction))
    assertthat::assert_that(length(interaction) == 1)
    ##
    assertthat::assert_that(is.logical(casecontrol))
    assertthat::assert_that(length(casecontrol) == 1)
    ##
    assertthat::assert_that(is.null(eventvar) | is.character(eventvar))
    assertthat::assert_that(length(eventvar) <= 1)
    
    ## Do not allow missing data in variables of interest
    ## because they may differ unexpected differences in sample sizes
    ## between yreg and mreg.
    vars_interest <- c(yvar, avar, mvar, cvar, eventvar)
    data_vars_interest <- data[, vars_interest, drop = FALSE]
    logical_vars_missing <- unlist(lapply(data_vars_interest, function(v) {
        any(is.na(v))
    }))
    ## Capture names of variables with missingness
    vars_missing <- vars_interest[logical_vars_missing]
    assertthat::assert_that(all(stats::complete.cases(data_vars_interest)),
                            msg = paste0("Missingness is not allowed in variables of interest!
For multiple imputation, see the multiple imputation vignette: vignette(\"vig_04_mi\")
To perform complete case analysis, use na_omit = TRUE.
Variables with missingness: ",
                                         paste0(vars_missing, collapse = ", ")))
    
    ## Do not allow factors as they can result in multiple
    ## dummy variables and coef results in different names
    ## from the original variable names.
    assertthat::assert_that(all(!vapply(X = data_vars_interest,
                                        FUN = is.factor,
                                        ## template value
                                        FUN.VALUE = TRUE)),
                            msg = "Factors are not allowed! Use numeric variables only. Create multiple dichotomous variables for multi-category variables.")
    
    NULL
}


###
### Argument validation function
################################################################################

##' Validate soundness of a regmedint object.
##'
##' Check the structure of a proposed regmedint object for soundness.
##'
##' @param x A \code{regmedint} object.
##'
##' @return No return value, called for side effects.
validate_regmedint <- function(x) {
    
    assertthat::assert_that(class(x)[[1]] == "regmedint")
    ##
    assertthat::assert_that("myreg" %in% names(x))
    assertthat::assert_that("args" %in% names(x))
    ##
    assertthat::assert_that(is.list(x$args))
    
    NULL
}



