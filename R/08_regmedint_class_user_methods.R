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
### print method
################################################################################

##' print method for regmedint object
##'
##' Print the \code{mreg_fit}, \code{yreg_fit}, and the mediation analysis effect estimates.
##'
##' @param x An object of the \code{\link{regmedint}} class.
##' @param a0 A numeric vector of length one.
##' @param a1 A numeric vector of length one.
##' @param m_cde A numeric vector of length one. A mediator value at which the controlled direct effect (CDE) conditional on the adjustment covariates is evaluated. If not provided, the default value supplied to the call to \code{\link{regmedint}} will be used. Only the CDE is affected.
##' @param c_cond A numeric vector as long as the number of adjustment covariates. A set of covariate values at which the conditional natural effects are evaluated.
##' @param args_mreg_fit A named list of argument to be passed to the method for the \code{mreg_fit} object.
##' @param args_yreg_fit A named list of argument to be passed to the method for the \code{mreg_fit} object.
##' @param ... For compatibility with the generic. Ignored.
##'
##' @return Invisibly return the \code{regmedint} class object as is.
print.regmedint <- function(x,
                            a0 = NULL,
                            a1 = NULL,
                            m_cde = NULL,
                            c_cond = NULL,
                            args_mreg_fit = list(),
                            args_yreg_fit = list(),
                            ...) {

    ## This is a user function. Check arguments heavily.
    assertthat::assert_that(is.null(a0) | (length(a0) == 1))
    assertthat::assert_that(is.null(a1) | (length(a1) == 1))
    assertthat::assert_that(is.null(m_cde) | (length(m_cde) == 1))
    assertthat::assert_that(is.null(c_cond) | (length(c_cond) == length(x$args$cvar)))
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
##' @param exponentiate Whether to show exponentiated point estimates.
##'
##' @return A numeric matrix corresponding to what is displayed.
summary.regmedint <- function(x,
                              a0 = NULL,
                              a1 = NULL,
                              m_cde = NULL,
                              c_cond = NULL,
                              args_mreg_fit = list(),
                              args_yreg_fit = list(),
                              exponentiate = FALSE,
                              ...) {

    ## This is a user function. Check arguments heavily.
    assertthat::assert_that(is.null(a0) | (length(a0) == 1))
    assertthat::assert_that(is.null(a1) | (length(a1) == 1))
    assertthat::assert_that(is.null(m_cde) | (length(m_cde) == 1))
    assertthat::assert_that(is.null(c_cond) | (length(c_cond) == length(x$args$cvar)))
    assertthat::assert_that(is.list(args_mreg_fit))
    assertthat::assert_that(is.list(args_yreg_fit))

    cat("### Mediator model\n")
    do.call(summary, c(list(x$mreg),
                       args_mreg_fit))

    cat("### Outcome model\n")
    do.call(summary, c(list(x$yreg),
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
    res_p <- pnorm(q = res_Z)

    if (exponentiate & x$args$yreg != "linear") {
        ## Note only the point estimates are are exponentiated.
        res_mat <- cbind(`exp(est)` = exp(res_est),
                         `SE(est)` = res_se,
                         Z = res_Z,
                         p = res_p)
    } else {
        res_mat <- cbind(est = res_est,
                         `SE(est)` = res_se,
                         Z = res_Z,
                         p = res_p)
    }

    attr(res_mat, which = "args") <- list(a0 = a0,
                                          a1 = a1,
                                          m_cde = m_cde,
                                          c_cond = c_cond)

    print(res_mat, quote = FALSE)
    invisible(res_mat)
}


###
### Others
################################################################################

##' Extract point estimates.
##'
##' Extract point estimates evaluated at \code{a0}, \code{a1}, \code{m_cde}, and \code{c_cond}.
##'
##' @inheritParams print.regmedint
##'
##' @return A numeric vector of point estimates.
coef.regmedint <- function(x,
                           a0 = NULL,
                           a1 = NULL,
                           m_cde = NULL,
                           c_cond = NULL,
                           ...) {

    ## This is a user function. Check arguments heavily.
    assertthat::assert_that(is.null(a0) | (length(a0) == 1))
    assertthat::assert_that(is.null(a1) | (length(a1) == 1))
    assertthat::assert_that(is.null(m_cde) | (length(m_cde) == 1))
    assertthat::assert_that(is.null(c_cond) | (length(c_cond) == length(x$args$cvar)))

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

##' Confidence intervals for mediation prameter estimates.
##'
##' .. content for \details{} ..
##'
##' @param x
##' @param a0
##' @param a1
##' @param m_cde
##' @param c_cond
##' @param alpha
##' @param ...
##'
##' @return A numeric matrix
confint.regmedint <- function(x,
                              a0 = NULL,
                              a1 = NULL,
                              m_cde = NULL,
                              c_cond = NULL,
                              alpha = 0.05,
                              ...) {

    ## This is a user function. Check arguments heavily.
    assertthat::assert_that(is.null(a0) | (length(a0) == 1))
    assertthat::assert_that(is.null(a1) | (length(a1) == 1))
    assertthat::assert_that(is.null(m_cde) | (length(m_cde) == 1))
    assertthat::assert_that(is.null(c_cond) | (length(c_cond) == length(x$args$cvar)))

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

    res_mat <- cbind(lower = res_est - qnorm(p = (1 - (alpha / 2))) * res_se,
                     upper = res_est + qnorm(p = (1 - (alpha / 2))) * res_se)


    ## Set as attributes.
    attr(res_mat, which = "args") <- list(a0 = a0,
                                          a1 = a1,
                                          m_cde = m_cde,
                                          c_cond = c_cond)

    res_mat
}
