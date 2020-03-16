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
##' @param x An object of the \code{regmedint} class.
##' @param a0 A numeric vector of length one.
##' @param a1 A numeric vector of length one.
##' @param m_cde A numeric vector of length one. A mediator value at which the controlled direct effect (CDE) conditional on the adjustment covariates is evaluated. If not provided, the default value supplied to the call to \code{\link{regmedint}} will be used. Only the CDE is affected.
##' @param c_cond A numeric vector as long as the number of adjustment covariates. A set of covariate values at which the conditional natural effects are evaluated.
##' @param ...
##'
##' @return Invisibly return the \code{regmedint} class object as is.
print.regmedint <- function(x,
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

    cat("### Mediator model\n")
    print(x$mreg)

    cat("### Outcome model\n")
    print(x$yreg)

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
    print(x$myreg$est_fun(a0 = a0, a1 = a1, m_cde = m_cde, c_cond = c_cond))

    invisible(x)
}


###
### summary method
################################################################################

##' summary method for regmedint object
##'
##' Show summary of each model.
##'
##' @param x
##' @param a0
##' @param a1
##' @param m_cde
##' @param c_cond
##' @param ...
##'
##' @return
summary.regmedint <- function(x,
                              a0 = NULL,
                              a1 = NULL,
                              m_cde = NULL,
                              c_cond = NULL,
                              exponentiate = FALSE,
                              ...) {

    ## This is a user function. Check arguments heavily.
    assertthat::assert_that(is.null(a0) | (length(a0) == 1))
    assertthat::assert_that(is.null(a1) | (length(a1) == 1))
    assertthat::assert_that(is.null(m_cde) | (length(m_cde) == 1))
    assertthat::assert_that(is.null(c_cond) | (length(c_cond) == length(x$args$cvar)))

    cat("### Mediator model\n")
    summary(x$mreg)

    cat("### Outcome model\n")
    summary(x$yreg)

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

    print(res_mat, quote = FALSE)
    invisible(res_mat)
}


###
### Others
################################################################################

##' Extract coefficients
##'
##' .. content for \details{} ..
##'
##' @param x
##' @param a0
##' @param a1
##' @param m_cde
##' @param c_cond
##' @param ...
##' @return
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

    res_est <- x$myreg$est_fun(a0 = a0, a1 = a1, m_cde = m_cde, c_cond = c_cond)

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
##' @return
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

    res_est <- x$myreg$est_fun(a0 = a0, a1 = a1, m_cde = m_cde, c_cond = c_cond)
    res_se <- x$myreg$se_fun(a0 = a0, a1 = a1, m_cde = m_cde, c_cond = c_cond)

    res_mat <- cbind(lower = res_est - qnorm(p = (1 - alpha / 2)) * res_se,
                     upper = res_est + qnorm(p = (1 - alpha / 2)) * res_se)


    ## Set as attributes.
    attrr(res_mat, args = list(a0 = a0, a1 = a1, m_cde = m_cde, c_cond = c_cond))

    res_mat
}
