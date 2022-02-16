################################################################################
### User interface methods
##
## Created on: 2020-03-14
## Author: Kazuki Yoshida
################################################################################

## Advanced R. 13.4 Generics and methods
## https://adv-r.hadley.nz/s3.html#s3-methods
## sloop::s3_dispatch(print(x))to check the method being called.
## sloop::s3_methods_class("regmedint") to check all methods for regmedint.


###
### Helper
################################################################################

validate_eval_at_values <- function(x, a0, a1, m_cde, c_cond) {
    assertthat::assert_that(is.null(a0) | (length(a0) == 1))
    assertthat::assert_that(is.null(a1) | (length(a1) == 1))
    assertthat::assert_that(is.null(m_cde) | (length(m_cde) == 1))
    assertthat::assert_that(is.null(c_cond) | (length(c_cond) == length(x$args$cvar)))
}

###
### print method
################################################################################

##' print method for regmedint object
##'
##' Print the \code{mreg_fit}, \code{yreg_fit}, and the mediation analysis effect estimates.
##'
##' @param x An object of the \code{\link{regmedint}} class.
##' @param a0 A numeric vector of length 1
##' @param a1 A numeric vector of length 1
##' @param m_cde A numeric vector of length 1 The mediator value at which the controlled direct effect (CDE) conditional on the adjustment covariates is evaluated. If not provided, the default value supplied to the call to \code{\link{regmedint}} will be used. Only the CDE is affected.
##' @param c_cond A numeric vector of the same length as \code{cvar}. A set of covariate values at which the conditional natural effects are evaluated.
##' @param args_mreg_fit A named list of argument to be passed to the method for the \code{mreg_fit} object.
##' @param args_yreg_fit A named list of argument to be passed to the method for the \code{mreg_fit} object.
##' @param ... For compatibility with the generic. Ignored.
##'
##' @return Invisibly return the \code{regmedint} class object as is.
##'
##' @examples
##' library(regmedint)
##' data(vv2015)
##' regmedint_obj <- regmedint(data = vv2015,
##'                            ## Variables
##'                            yvar = "y",
##'                            avar = "x",
##'                            mvar = "m",
##'                            cvar = c("c"),
##'                            eventvar = "event",
##'                            ## Values at which effects are evaluated
##'                            a0 = 0,
##'                            a1 = 1,
##'                            m_cde = 1,
##'                            c_cond = 0.5,
##'                            ## Model types
##'                            mreg = "logistic",
##'                            yreg = "survAFT_weibull",
##'                            ## Additional specification
##'                            interaction = TRUE,
##'                            casecontrol = FALSE)
##' ## Implicit printing
##' regmedint_obj
##' ## Explicit printing
##' print(regmedint_obj)
##' ## Evaluate at different values
##' print(regmedint_obj, m_cde = 0, c_cond = 1)
##'
##' @export
print.regmedint <- function(x,
                            a0 = NULL,
                            a1 = NULL,
                            m_cde = NULL,
                            c_cond = NULL,
                            args_mreg_fit = list(),
                            args_yreg_fit = list(),
                            ...) {

    ## This is a user function. Check arguments heavily.
    validate_eval_at_values(x = x,
                            a0 = a0,
                            a1 = a1,
                            m_cde = m_cde,
                            c_cond = c_cond)
    assertthat::assert_that(is.list(args_mreg_fit))
    assertthat::assert_that(is.list(args_yreg_fit))

    cat("### Mediator model\n")
    do.call(print, c(list(x$mreg),
                     args_mreg_fit))

    cat("### Outcome model\n")
    do.call(print, c(list(x$yreg),
                     args_yreg_fit))

    cat("### Mediation analysis \n")
    if (is.null(a0)) {
        a0 <- x$args$a0
    }
    if (is.null(a1)) {
        a1 <- x$args$a1
    }
    if (is.null(m_cde)) {
        m_cde <- x$args$m_cde
    }
    if (is.null(c_cond)) {
        c_cond <- x$args$c_cond
    }
    ## Compute point estimates
    res <- x$myreg$est_fun(a0 = a0,
                           a1 = a1,
                           m_cde = m_cde,
                           c_cond = c_cond)
    assertthat::assert_that(is.numeric(res),
                            is.vector(res))
    print(res)

    invisible(x)
}


###
### summary method
################################################################################

##' summary method for regmedint object
##'
##' Summarize the \code{mreg_fit}, \code{yreg_fit}, and the mediation analysis effect estimates.
##'
##' @inheritParams print.regmedint
##' @param object An object of the \code{\link{regmedint}} class.
##' @param exponentiate Whether to add exponentiated point and confidence limit estimates. When \code{yreg = "linear"}, it is ignored.
##' @param level Confidence level for the confidence intervals.
##'
##' @return A \code{summary_regmedint} object, which is a list containing the summary objects of the \code{mreg_fit} and the \code{yreg_fit} as well as the mediation analysis results.
##'
##' @examples
##' library(regmedint)
##' data(vv2015)
##' regmedint_obj <- regmedint(data = vv2015,
##'                            ## Variables
##'                            yvar = "y",
##'                            avar = "x",
##'                            mvar = "m",
##'                            cvar = c("c"),
##'                            eventvar = "event",
##'                            ## Values at which effects are evaluated
##'                            a0 = 0,
##'                            a1 = 1,
##'                            m_cde = 1,
##'                            c_cond = 0.5,
##'                            ## Model types
##'                            mreg = "logistic",
##'                            yreg = "survAFT_weibull",
##'                            ## Additional specification
##'                            interaction = TRUE,
##'                            casecontrol = FALSE)
##' ## Detailed result with summary
##' summary(regmedint_obj)
##' ## Add exponentiate results for non-linear outcome models
##' summary(regmedint_obj, exponentiate = TRUE)
##' ## Evaluate at different values
##' summary(regmedint_obj, m_cde = 0, c_cond = 1)
##' ## Change confidence level
##' summary(regmedint_obj, m_cde = 0, c_cond = 1, level = 0.99)
##'
##' @export
summary.regmedint <- function(object,
                              a0 = NULL,
                              a1 = NULL,
                              m_cde = NULL,
                              c_cond = NULL,
                              args_mreg_fit = list(),
                              args_yreg_fit = list(),
                              exponentiate = FALSE,
                              level = 0.95,
                              ...) {

    x <- object
    ## This is a user function. Check arguments heavily.
    validate_eval_at_values(x = x,
                            a0 = a0,
                            a1 = a1,
                            m_cde = m_cde,
                            c_cond = c_cond)
    assertthat::assert_that(is.list(args_mreg_fit))
    assertthat::assert_that(is.list(args_yreg_fit))

    ## Mediator model
    summary_mreg_fit <- do.call(summary, c(list(x$mreg),
                                           args_mreg_fit))

    ## Outcome model
    summary_yreg_fit <- do.call(summary, c(list(x$yreg),
                                           args_yreg_fit))

    ## Mediation analysis result matrix construction
    if (is.null(a0)) {
        a0 <- x$args$a0
    }
    if (is.null(a1)) {
        a1 <- x$args$a1
    }
    if (is.null(m_cde)) {
        m_cde <- x$args$m_cde
    }
    if (is.null(c_cond)) {
        c_cond <- x$args$c_cond
    }
    ## Compute point estimates
    res_est <- x$myreg$est_fun(a0 = a0,
                               a1 = a1,
                               m_cde = m_cde,
                               c_cond = c_cond)
    ## Compute SE estimates
    res_se <- x$myreg$se_fun(a0 = a0,
                             a1 = a1,
                             m_cde = m_cde,
                             c_cond = c_cond)

    assertthat::assert_that(length(res_est) == length(res_se))
    res_Z <- res_est / res_se
    res_p <- 2 * (1 - pnorm(q = abs(res_Z)))

    ## Compute CI estimates
    res_ci <- confint(object = x,
                      a0 = a0,
                      a1 = a1,
                      m_cde = m_cde,
                      c_cond = c_cond,
                      level = level)

    if (exponentiate & x$args$yreg != "linear") {
        res_mat <- cbind(est = res_est,
                         se = res_se,
                         Z = res_Z,
                         p = res_p,
                         lower = res_ci[,"lower"],
                         upper = res_ci[,"upper"],
                         `exp(est)` = exp(res_est),
                         `exp(lower)` = exp(res_ci[,"lower"]),
                         `exp(upper)` = exp(res_ci[,"upper"]))
        ## exp columns for pm should be NA.
        res_mat["pm",c("exp(est)","exp(lower)","exp(upper)")] <- NA
    } else {
        res_mat <- cbind(est = res_est,
                         se = res_se,
                         Z = res_Z,
                         p = res_p,
                         lower = res_ci[,"lower"],
                         upper = res_ci[,"upper"])
    }

    res <- list(summary_mreg_fit = summary_mreg_fit,
                summary_yreg_fit = summary_yreg_fit,
                summary_myreg = res_mat,
                eval_at = list(a0 = a0,
                               a1 = a1,
                               m_cde = m_cde,
                               c_cond = c_cond),
                args = x$args)

    class(res) <- c("summary_regmedint", class(res))
    return(res)
}

##' Print method for summary objects from \code{\link{summary.regmedint}}
##'
##' Print results contained in a \code{summary_regmedint} object with additional explanation regarding the evaluation settings.
##'
##' @param x An object of the class \code{summary_regmedint}.
##' @param ... For compatibility with the generic function.
##'
##' @return Invisibly return the first argument.
##'
##' @examples
##' library(regmedint)
##' data(vv2015)
##' regmedint_obj <- regmedint(data = vv2015,
##'                            ## Variables
##'                            yvar = "y",
##'                            avar = "x",
##'                            mvar = "m",
##'                            cvar = c("c"),
##'                            eventvar = "event",
##'                            ## Values at which effects are evaluated
##'                            a0 = 0,
##'                            a1 = 1,
##'                            m_cde = 1,
##'                            c_cond = 0.5,
##'                            ## Model types
##'                            mreg = "logistic",
##'                            yreg = "survAFT_weibull",
##'                            ## Additional specification
##'                            interaction = TRUE,
##'                            casecontrol = FALSE)
##' ## Implicit printing
##' summary(regmedint_obj)
##' ## Explicit printing
##' print(summary(regmedint_obj))
##'
##' @export
print.summary_regmedint <- function(x, ...) {


    cat("### Mediator model\n")
    print(x$summary_mreg_fit)

    cat("### Outcome model\n")
    print(x$summary_yreg_fit)

    cat("### Mediation analysis \n")
    res_mat <- x$summary_myreg

    ## Print before assignment of attributes for cleaness.
    print(res_mat, quote = FALSE)

    ## Print helpful information here
    print_eval_info_helper(a0 = x$eval_at$a0,
                           a1 = x$eval_at$a1,
                           m_cde = x$eval_at$m_cde,
                           c_cond = x$eval_at$c_cond,
                           ## Orignal arguments
                           avar = x$args$avar,
                           mvar = x$args$mvar,
                           cvar = x$args$cvar,
                           emm_ac_mreg =  x$args$emm_ac_mreg,
                           emm_ac_yreg = x$args$emm_ac_yreg,
                           emm_mc_yreg = x$args$emm_mc_yreg,
                           yreg = x$args$yreg,
                           mreg = x$args$mreg,
                           interaction = x$args$interaction,
                           casecontrol = x$args$casecontrol)

    invisible(x)
}


##' Extract the result matrix from a summary_regmedint object.
##'
##' Extract the result matrix from a summary_regmedint object.
##'
##' @param object An object with a class of \code{summary_regmedint}.
##' @param ... For compatibility with the generic.
##'
##' @return  A matrix populated with results.
##'
##' @examples
##' library(regmedint)
##' data(vv2015)
##' regmedint_obj <- regmedint(data = vv2015,
##'                            ## Variables
##'                            yvar = "y",
##'                            avar = "x",
##'                            mvar = "m",
##'                            cvar = c("c"),
##'                            eventvar = "event",
##'                            ## Values at which effects are evaluated
##'                            a0 = 0,
##'                            a1 = 1,
##'                            m_cde = 1,
##'                            c_cond = 0.5,
##'                            ## Model types
##'                            mreg = "logistic",
##'                            yreg = "survAFT_weibull",
##'                            ## Additional specification
##'                            interaction = TRUE,
##'                            casecontrol = FALSE)
##' coef(summary(regmedint_obj))
##'
##' @export
coef.summary_regmedint <- function(object, ...) {
    object$summary_myreg
}


print_eval_info_helper <- function(a0, a1, m_cde, c_cond,
                                   avar, mvar, cvar, emm_ac_mreg, emm_ac_yreg, emm_mc_yreg,
                                   yreg, mreg, interaction, casecontrol) {

    cat("\n")
    cat("Evaluated at:\n")

    ## avar
    cat("avar: ")
    cat(avar)
    cat("\n")

    cat(" a1 (intervened value of avar) = ")
    cat(a1)
    cat("\n")

    cat(" a0 (reference value of avar)  = ")
    cat(a0)
    cat("\n")

    ## mvar
    cat("mvar: ")
    cat(mvar)
    cat("\n")

    cat(" m_cde (intervend value of mvar for cde) = ")
    cat(m_cde)
    cat("\n")

    ## cvar
    cat("cvar: ")
    cat(paste(cvar, collapse = " "))
    cat("\n")

    cat(" c_cond (covariate vector value) = ")
    cat(c_cond)
    cat("\n")
    cat("\n")

    if (interaction) {
        cat("Note that effect estimates can vary over m_cde and c_cond values when interaction = TRUE.\n")
    } else {
        cat("Note that effect estimates do not vary over m_cde and c_cond values when interaction = FALSE.\n")
    }
}


###
### Others
################################################################################

##' Extract point estimates.
##'
##' Extract point estimates evaluated at \code{a0}, \code{a1}, \code{m_cde}, and \code{c_cond}.
##'
##' @inheritParams print.regmedint
##' @param object An object of the \code{\link{regmedint}} class.
##'
##' @return A numeric vector of point estimates.
##'
##' @examples
##' library(regmedint)
##' data(vv2015)
##' regmedint_obj <- regmedint(data = vv2015,
##'                            ## Variables
##'                            yvar = "y",
##'                            avar = "x",
##'                            mvar = "m",
##'                            cvar = c("c"),
##'                            eventvar = "event",
##'                            ## Values at which effects are evaluated
##'                            a0 = 0,
##'                            a1 = 1,
##'                            m_cde = 1,
##'                            c_cond = 0.5,
##'                            ## Model types
##'                            mreg = "logistic",
##'                            yreg = "survAFT_weibull",
##'                            ## Additional specification
##'                            interaction = TRUE,
##'                            casecontrol = FALSE)
##' coef(regmedint_obj)
##' ## Evaluate at different values
##' coef(regmedint_obj, m_cde = 0, c_cond = 1)
##'
##' @export
coef.regmedint <- function(object,
                           a0 = NULL,
                           a1 = NULL,
                           m_cde = NULL,
                           c_cond = NULL,
                           ...) {

    x <- object
    ## This is a user function. Check arguments heavily.
    validate_eval_at_values(x = x,
                            a0 = a0,
                            a1 = a1,
                            m_cde = m_cde,
                            c_cond = c_cond)

    if (is.null(a0)) {
        a0 <- x$args$a0
    }
    if (is.null(a1)) {
        a1 <- x$args$a1
    }
    if (is.null(m_cde)) {
        m_cde <- x$args$m_cde
    }
    if (is.null(c_cond)) {
        c_cond <- x$args$c_cond
    }

    res_est <- x$myreg$est_fun(a0 = a0,
                               a1 = a1,
                               m_cde = m_cde,
                               c_cond = c_cond)

    attr(res_est, which = "args") <- list(a0 = a0,
                                          a1 = a1,
                                          m_cde = m_cde,
                                          c_cond = c_cond)

    res_est
}

##' Extract variance estimates in the vcov form.
##'
##' Extract variance estimates evaluated at \code{a0}, \code{a1}, \code{m_cde}, and \code{c_cond}.
##'
##' @inheritParams print.regmedint
##' @param object An object of the \code{\link{regmedint}} class.
##'
##' @return A numeric matrix with the diagonals populated with variance estimates. Off-diagnonals are NA since these are not estimated.
##'
##' @examples
##' library(regmedint)
##' data(vv2015)
##' regmedint_obj <- regmedint(data = vv2015,
##'                            ## Variables
##'                            yvar = "y",
##'                            avar = "x",
##'                            mvar = "m",
##'                            cvar = c("c"),
##'                            eventvar = "event",
##'                            ## Values at which effects are evaluated
##'                            a0 = 0,
##'                            a1 = 1,
##'                            m_cde = 1,
##'                            c_cond = 0.5,
##'                            ## Model types
##'                            mreg = "logistic",
##'                            yreg = "survAFT_weibull",
##'                            ## Additional specification
##'                            interaction = TRUE,
##'                            casecontrol = FALSE)
##' vcov(regmedint_obj)
##' ## Evaluate at different values
##' vcov(regmedint_obj, m_cde = 0, c_cond = 1)
##'
##' @export
vcov.regmedint <- function(object,
                           a0 = NULL,
                           a1 = NULL,
                           m_cde = NULL,
                           c_cond = NULL,
                           ...) {

    x <- object
    ## This is a user function. Check arguments heavily.
    validate_eval_at_values(x = x,
                            a0 = a0,
                            a1 = a1,
                            m_cde = m_cde,
                            c_cond = c_cond)

    if (is.null(a0)) {
        a0 <- x$args$a0
    }
    if (is.null(a1)) {
        a1 <- x$args$a1
    }
    if (is.null(m_cde)) {
        m_cde <- x$args$m_cde
    }
    if (is.null(c_cond)) {
        c_cond <- x$args$c_cond
    }

    ## Compute SE estimates
    res_se <- x$myreg$se_fun(a0 = a0,
                             a1 = a1,
                             m_cde = m_cde,
                             c_cond = c_cond)

    ## Construct matrix with the diagonals as variances.
    res_vcov <- diag(res_se^2)

    ## NA out off diagonals as the are not zeros.
    res_vcov[upper.tri(res_vcov)] <- NA
    res_vcov[lower.tri(res_vcov)] <- NA

    ## Name
    colnames(res_vcov) <- names(coef(x))
    rownames(res_vcov) <- names(coef(x))

    attr(res_vcov, which = "args") <- list(a0 = a0,
                                           a1 = a1,
                                           m_cde = m_cde,
                                           c_cond = c_cond)
    res_vcov
}

##' Confidence intervals for mediation prameter estimates.
##'
##' Construct Wald approximate confidence intervals for the quantities of interest.
##'
##' @inheritParams print.regmedint
##' @param object An object of the \code{\link{regmedint}} class.
##' @param parm For compatibility with generic. Ignored.
##' @param level A numeric vector of length one. Requested confidence level. Defaults to 0.95.
##' @param ... For compatibility with generic.
##'
##' @return A numeric matrix of the lower limit and upper limit.
##'
##' @examples
##' library(regmedint)
##' data(vv2015)
##' regmedint_obj <- regmedint(data = vv2015,
##'                            ## Variables
##'                            yvar = "y",
##'                            avar = "x",
##'                            mvar = "m",
##'                            cvar = c("c"),
##'                            eventvar = "event",
##'                            ## Values at which effects are evaluated
##'                            a0 = 0,
##'                            a1 = 1,
##'                            m_cde = 1,
##'                            c_cond = 0.5,
##'                            ## Model types
##'                            mreg = "logistic",
##'                            yreg = "survAFT_weibull",
##'                            ## Additional specification
##'                            interaction = TRUE,
##'                            casecontrol = FALSE)
##' confint(regmedint_obj)
##' ## Evaluate at different values
##' confint(regmedint_obj, m_cde = 0, c_cond = 1)
##' ## Change confidence level
##' confint(regmedint_obj, m_cde = 0, c_cond = 1, level = 0.99)
##'
##' @export
confint.regmedint <- function(object,
                              parm = NULL,
                              level = 0.95,
                              a0 = NULL,
                              a1 = NULL,
                              m_cde = NULL,
                              c_cond = NULL,
                              ...) {

    x <- object
    ## This is a user function. Check arguments heavily.
    validate_eval_at_values(x = x,
                            a0 = a0,
                            a1 = a1,
                            m_cde = m_cde,
                            c_cond = c_cond)

    if (is.null(a0)) {
        a0 <- x$args$a0
    }
    if (is.null(a1)) {
        a1 <- x$args$a1
    }
    if (is.null(m_cde)) {
        m_cde <- x$args$m_cde
    }
    if (is.null(c_cond)) {
        c_cond <- x$args$c_cond
    }
    ## Compute point estimates
    res_est <- x$myreg$est_fun(a0 = a0,
                               a1 = a1,
                               m_cde = m_cde,
                               c_cond = c_cond)
    ## Compute SE estimates
    res_se <- x$myreg$se_fun(a0 = a0,
                             a1 = a1,
                             m_cde = m_cde,
                             c_cond = c_cond)

    alpha <- 1 - level
    res_mat <- cbind(lower = res_est - qnorm(p = (1 - (alpha / 2))) * res_se,
                     upper = res_est + qnorm(p = (1 - (alpha / 2))) * res_se)


    ## Set as attributes.
    attr(res_mat, which = "args") <- list(a0 = a0,
                                          a1 = a1,
                                          m_cde = m_cde,
                                          c_cond = c_cond)

    res_mat
}
