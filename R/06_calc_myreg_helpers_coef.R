################################################################################
### Helper functions for regression-based causal mediation analysis
##
## Created on: 2020-03-11
## Author: Kazuki Yoshida
################################################################################


###
### coef extractors
################################################################################

##' Create a vector of coefficients from the mediator model (mreg)
##'
##' This function extracts \code{\link{coef}} from \code{mreg_fit} and pads with zeros appropriately to create a named vector consistently having the following elements:
##' \code{(Intercept)},
##' \code{avar},
##' \code{cvar} (this part is eliminated when \code{cvar = NULL}),
##' \code{emm_ac_mreg} (this part is eliminated when \code{emm_ac_mreg = NULL}).
##'
##' @inheritParams regmedint
##' @param mreg_fit Model fit object for mreg (mediator model).
##'
##' @return A named numeric vector of coefficients.
beta_hat <- function(mreg, mreg_fit, avar, cvar, emm_ac_mreg = NULL) {
    if (is.null(cvar)) {
        vars <- c("(Intercept)", avar)
    } else if (!is.null(cvar) & is.null(emm_ac_mreg)){
        vars <- c("(Intercept)", avar, cvar)
    } else if (!is.null(cvar) & !is.null(emm_ac_mreg)){
        vars <- c("(Intercept)", avar, cvar, paste0(avar, ":", emm_ac_mreg))
    }
    coef(mreg_fit)[vars]
}



beta_hat_helper <- function(mreg, mreg_fit, avar, cvar, emm_ac_mreg = NULL) {
    beta_hat <- beta_hat(mreg = mreg,
                         mreg_fit = mreg_fit,
                         avar = avar,
                         cvar = cvar,
                         emm_ac_mreg = emm_ac_mreg)
    beta0 <- beta_hat["(Intercept)"]
    beta1 <- beta_hat[avar]
    if (is.null(cvar)) {
        ## beta_hat does not contain the cvar part in this case. 
        beta2 <- NULL
    } else {
        beta2 <- beta_hat[cvar]
    }
    if (is.null(emm_ac_mreg)){
        beta3 <- NULL
    } else {
        ## initialize to 0
        beta3 <- rep(0, length(cvar))
        names(beta3) <- paste0(avar, ":", cvar)
        ## plug in non-zeros to corresponding elements
        for(var_name in emm_ac_mreg){
          beta3[paste0(avar, ":", var_name)] <- beta_hat[paste0(avar, ":", var_name)] 
        }
    }
    list(beta0 = beta0,
         beta1 = beta1,
         beta2 = beta2,
         beta3 = beta3)
    
}

##' Create a vector of coefficients from the outcome model (yreg)
##'
##' This function extracts \code{\link{coef}} from \code{yreg_fit} and 3s with zeros appropriately to create a named vector consistently having the following elements:
##' \code{(Intercept)} (a zero element is added for \code{yreg = "survCox"} for which no intercept is estimated (the baseline hazard is left unspecified)),
##' \code{avar},
##' \code{mvar},
##' \code{avar:mvar} (a zero element is added when \code{interaction = FALSE}).
##' \code{cvar} (this part is eliminated when \code{cvar = NULL}),
##' \code{emm_ac_yreg} (this part is eliminated when \code{emm_ac_yreg = NULL}),
##' \code{emm_mc_yreg} (this part is eliminated when \code{emm_mc_yreg = NULL}).
##'
##' @inheritParams regmedint
##' @param yreg_fit Model fit object for yreg (outcome model).
##'
##' @return A named numeric vector of coefficients.
theta_hat <- function(yreg, yreg_fit, avar, mvar, cvar, 
                      emm_ac_yreg = NULL, emm_mc_yreg = NULL, 
                      interaction) {

    coef_raw <- coef(yreg_fit)

    ## Handle the absence of (Intercept) for Cox regression
    if (yreg == "survCox") {

        ## Pad with a zero for the missing (Intercept)
        coef_ready <- c(0, coef_raw)
        names(coef_ready) <- c("(Intercept)", names(coef_raw))

    } else {

        coef_ready <- coef_raw

    }

    ## Construct vars to extract from coef_ready
    ## Make sure theta3 for avar:mvar always stays
    if (interaction) {

        ## Interaction case
        ## No data manipulation is necessary.
        ## Technically, the first case can be used in both because NULL
        ## drops out in c(..., NULL). Here it is made explicit.
        if (is.null(cvar)) {
          coef_ready <- c(coef_ready[c("(Intercept)", avar, mvar, paste0(avar, ":", mvar))])
          names(coef_ready) <- c("(Intercept)",
                                 avar, mvar,
                                 paste0(avar, ":", mvar))
          vars <- c("(Intercept)", avar, mvar, paste0(avar, ":", mvar))
        } else if (!is.null(cvar) & is.null(emm_ac_yreg) & is.null(emm_mc_yreg)) {
          coef_ready <- c(coef_ready[c("(Intercept)", avar, mvar, paste0(avar, ":", mvar), cvar)])
          names(coef_ready) <- c("(Intercept)",
                                   avar, mvar,
                                   paste0(avar, ":", mvar),
                                   cvar)
          vars <- c("(Intercept)", avar, mvar, paste0(avar, ":", mvar), cvar)
        } else if (!is.null(cvar) & !is.null(emm_ac_yreg) & is.null(emm_mc_yreg)) {
          coef_ready <- c(coef_ready[c("(Intercept)", avar, mvar, paste0(avar, ":", mvar), cvar,
                                         paste0(avar, ":", emm_ac_yreg))])
          names(coef_ready) <- c("(Intercept)",
                                 avar, mvar,
                                 paste0(avar, ":", mvar),
                                 cvar,
                                 paste0(avar, ":", emm_ac_yreg))
          vars <- c("(Intercept)", avar, mvar, paste0(avar, ":", mvar), cvar, 
                      paste0(avar, ":", emm_ac_yreg))
        } else if (!is.null(cvar) & is.null(emm_ac_yreg) & !is.null(emm_mc_yreg)) {
          coef_ready <- c(coef_ready[c("(Intercept)", avar, mvar, paste0(avar, ":", mvar), cvar,
                                       paste0(mvar, ":", emm_mc_yreg))])
          names(coef_ready) <- c("(Intercept)",
                                 avar, mvar,
                                 paste0(avar, ":", mvar),
                                 cvar,
                                 paste0(mvar, ":", emm_mc_yreg))
          vars <- c("(Intercept)", avar, mvar, paste0(avar, ":", mvar), cvar, 
                      paste0(mvar, ":", emm_mc_yreg))
        } else if (!is.null(cvar) & !is.null(emm_ac_yreg) & !is.null(emm_mc_yreg)) {
          coef_ready <- c(coef_ready[c("(Intercept)", avar, mvar, paste0(avar, ":", mvar), cvar,
                                       paste0(avar, ":", emm_ac_yreg),
                                       paste0(mvar, ":", emm_mc_yreg))])
          names(coef_ready) <- c("(Intercept)",
                                 avar, mvar,
                                 paste0(avar, ":", mvar),
                                 cvar,
                                 paste0(avar, ":", emm_ac_yreg),
                                 paste0(mvar, ":", emm_mc_yreg))
          vars <- c("(Intercept)", avar, mvar, paste0(avar, ":", mvar), cvar, 
                      paste0(avar, ":", emm_ac_yreg), 
                      paste0(mvar, ":", emm_mc_yreg))
        }
      
      coef_ready[vars]
      
    } else {

        ## No interaction case
        if (is.null(cvar)) {
          coef_ready <- c(coef_ready[c("(Intercept)", avar, mvar)],
                          ## Add zero for AxM term
                          0)
          names(coef_ready) <- c("(Intercept)",
                                 avar, mvar,
                                 ## Name for AxM term
                                 paste0(avar, ":", mvar))
          vars <- c("(Intercept)", avar, mvar, paste0(avar, ":", mvar))
        } else if (!is.null(cvar) & is.null(emm_ac_yreg) & is.null(emm_mc_yreg)) {
            coef_ready <- c(coef_ready[c("(Intercept)", avar, mvar)],
                            ## Add zero for AxM term
                            0,
                            coef_ready[cvar])
            names(coef_ready) <- c("(Intercept)",
                                   avar, mvar,
                                   ## Name for AxM term
                                   paste0(avar, ":", mvar),
                                   cvar)
            vars <- c("(Intercept)", avar, mvar, paste0(avar, ":", mvar), cvar)
        } else if (!is.null(cvar) & !is.null(emm_ac_yreg) & is.null(emm_mc_yreg)){
          coef_ready <- c(coef_ready[c("(Intercept)", avar, mvar)],
                          ## Add zero for AxM term
                          0,
                          coef_ready[cvar],
                          coef_ready[paste0(avar, ":", emm_ac_yreg)])
          names(coef_ready) <- c("(Intercept)",
                                 avar, mvar,
                                 ## Name for AxM term
                                 paste0(avar, ":", mvar),
                                 cvar,
                                 ## Name for AxC term
                                 paste0(avar, ":", emm_ac_yreg))
          vars <- c("(Intercept)", avar, mvar, paste0(avar, ":", mvar), cvar, paste0(avar, ":", emm_ac_yreg))
        } else if (!is.null(cvar) & is.null(emm_ac_yreg) & !is.null(emm_mc_yreg)){
            coef_ready <- c(coef_ready[c("(Intercept)", avar, mvar)],
                            ## Add zero for AxM term
                            0,
                            coef_ready[cvar],
                            coef_ready[paste0(mvar, ":", emm_mc_yreg)])
            names(coef_ready) <- c("(Intercept)",
                                   avar, mvar,
                                   ## Name for AxM term
                                   paste0(avar, ":", mvar),
                                   cvar,
                                   ## Name for MxC term
                                   paste0(mvar, ":", emm_mc_yreg))
            vars <- c("(Intercept)", avar, mvar, paste0(avar, ":", mvar), cvar, paste0(mvar, ":", emm_mc_yreg))
        } else if (!is.null(cvar) & !is.null(emm_ac_yreg) & !is.null(emm_mc_yreg)){
          coef_ready <- c(coef_ready[c("(Intercept)", avar, mvar)],
                          ## Add zero for AxM term
                          0,
                          coef_ready[cvar],
                          coef_ready[paste0(avar, ":", emm_ac_yreg)],
                          coef_ready[paste0(mvar, ":", emm_mc_yreg)])
          names(coef_ready) <- c("(Intercept)",
                                 avar, mvar,
                                 ## Name for AxM term
                                 paste0(avar, ":", mvar),
                                 cvar,
                                 ## Name for AxC term
                                 paste0(avar, ":", emm_ac_yreg),
                                 ## Name for MxC term
                                 paste0(mvar, ":", emm_mc_yreg))
          vars <- c("(Intercept)", avar, mvar, paste0(avar, ":", mvar), cvar, paste0(avar, ":", emm_ac_yreg), paste0(mvar, ":", emm_mc_yreg))
        } 
      
      coef_ready[vars]
    }
    
}

theta_hat_helper <- function(yreg, yreg_fit, avar, mvar, cvar, emm_ac_yreg = NULL, emm_mc_yreg = NULL, interaction) {
    theta_hat <- theta_hat(yreg = yreg,
                           yreg_fit = yreg_fit,
                           avar = avar,
                           mvar = mvar,
                           cvar = cvar,
                           emm_ac_yreg = emm_ac_yreg,
                           emm_mc_yreg = emm_mc_yreg,
                           interaction = interaction)
    theta0 <- theta_hat["(Intercept)"]
    theta1 <- theta_hat[avar]
    theta2 <- theta_hat[mvar]
    theta3 <- theta_hat[paste0(avar,":", mvar)]
    if (is.null(cvar)) {
        theta4 <- NULL
    } else {
        theta4 <- theta_hat[cvar]
    }
    
    if (is.null(emm_ac_yreg)){
      theta5 <- NULL
    } else {
      ## initialize to 0
      theta5 <- rep(0, length(cvar))
      names(theta5) <- paste0(avar, ":", cvar)
      ## plug in non-zeros to corresponding elements
      for(var_name in emm_ac_yreg){
        theta5[paste0(avar, ":", var_name)] <- theta_hat[paste0(avar, ":", var_name)] 
      }
    }
    
    if (is.null(emm_mc_yreg)){
      theta6 <- NULL
    } else {
      ## initialize to 0
      theta6 <- rep(0, length(cvar))
      names(theta6) <- paste0(mvar, ":", cvar)
      ## plug in non-zeros to corresponding elements
      for(var_name in emm_mc_yreg){
        theta6[paste0(mvar, ":", var_name)] <- theta_hat[paste0(mvar, ":", var_name)] 
      }
    }
    
    list(theta0 = theta0,
         theta1 = theta1,
         theta2 = theta2,
         theta3 = theta3,
         theta4 = theta4,
         theta5 = theta5,
         theta6 = theta6)
    
}

sigma_hat_sq <- function(mreg_fit) {
    sigma(mreg_fit)^2
}

validate_myreg_coefs <- function(beta0,
                                 beta1,
                                 beta2,
                                 beta3,
                                 theta0,
                                 theta1,
                                 theta2,
                                 theta3,
                                 theta4,
                                 theta5,
                                 theta6,
                                 sigma_sq = NULL) {
    assertthat::assert_that(length(beta0) == 1)
    assertthat::assert_that(length(beta1) == 1)
    assertthat::assert_that(length(beta2) == length(theta4))
    # assertthat::assert_that(length(beta2) == length(beta3))
    assertthat::assert_that(length(theta0) == 1)
    assertthat::assert_that(length(theta1) == 1)
    assertthat::assert_that(length(theta2) == 1)
    assertthat::assert_that(length(theta3) == 1)
    # assertthat::assert_that(length(theta4) == length(theta5))
    # assertthat::assert_that(length(theta4) == length(theta6))
    if (!is.null(sigma_sq)) {
        assertthat::assert_that(length(sigma_sq) == 1)
    }
}

validate_myreg_vcovs <- function(beta0,
                                 beta1,
                                 beta2,
                                 beta3,
                                 theta0,
                                 theta1,
                                 theta2,
                                 theta3,
                                 theta4,
                                 theta5,
                                 theta6,
                                 sigma_sq = NULL,
                                 Sigma_beta,
                                 Sigma_theta,
                                 Sigma_sigma_sq = NULL) {

    Sigma_beta_size <- sum(1, # beta0 (Intercept)
                           1, # beta1 for avar
                           ## The following can be 0 (when cvar = NULL, length(NULL) = 0) 
                           length(beta2), 
                           length(beta3)) 
    # Can accommodate both scenarios:
    # if is.null(beta2), then Sigma_beta_size = sum(1,1); 
    # else, Sigma_beta_size = sum(1,1,length(cvar),length(cvar))
    assertthat::assert_that(dim(Sigma_beta)[1] == Sigma_beta_size)
    assertthat::assert_that(dim(Sigma_beta)[2] == Sigma_beta_size)

    Sigma_theta_size <- sum(1, # theta0 (Intercept). Never used so not in args.
                            1, # theta1 for avar
                            1, # theta2 for mvar
                            1, # theta3 for avar:mvar
                            ## The following can be 0 (when cvar = NULL, length(NULL) = 0) 
                            length(theta4),
                            length(theta5),
                            length(theta6)) 
    assertthat::assert_that(dim(Sigma_theta)[1] == Sigma_theta_size)
    assertthat::assert_that(dim(Sigma_theta)[2] == Sigma_theta_size)

    if (!is.null(Sigma_sigma_sq)) {
        assertthat::assert_that(dim(Sigma_sigma_sq)[1] == 1)
        assertthat::assert_that(dim(Sigma_sigma_sq)[2] == 1)
    }
}


###
### Proportion mediated helpers
################################################################################

##' Calculate the proportion mediated for yreg linear.
##'
##' Calculate the proportion mediated on the mean difference scale.
##'
##' @param pnde Pure natural direct effect.
##' @param tnie Total natural indirect effect.
##'
##' @return Proportion mediated value.
prop_med_yreg_linear <- function(pnde, tnie) {
    tnie / (pnde + tnie)
}

## Corresponding gradient: R^2 -> R^2
##' Calculate the gradient of the proportion mediated for yreg linear.
##'
##' Calculate the gradient of the proportion mediated for yreg linear case.
##'
##' @param pnde A numeric vector of length one. Pure natural direct effect.
##' @param tnie A numeric vector of length one. Total natural indirect effect.
##'
##' @return A numeric vector of length two. Gradient of the proportion mediated with respect to pnde and tnie.
grad_prop_med_yreg_linear <- Deriv::Deriv(prop_med_yreg_linear)
## function (pnde, tnie)
## {
##     .e1 <- pnde + tnie
##     c(pnde = -(tnie/.e1^2), tnie = (1 - tnie/.e1)/.e1)
## }

##' Calculate the proportion mediated for yreg logistic.
##'
##' Calculate the approximate proportion mediated on the risk difference scale.
##'
##' @param pnde Pure natural direct effect on the log scale.
##' @param tnie Total natural indirect effect on the log scale.
##'
##' @return Proportion mediated value.
prop_med_yreg_logistic <- function(pnde, tnie) {
    ## VanderWeele 2015. p48.
    (exp(pnde) * (exp(tnie) - 1)) / (exp(pnde) * exp(tnie) - 1)
}

## Corresponding gradient: R^2 -> R^2
##' Calculate the gradient of the proportion mediated for yreg logistic.
##'
##' Calculate the gradient of the proportion mediated for yreg logistic case.
##'
##' @param pnde A numeric vector of length one. Pure natural direct effect.
##' @param tnie A numeric vector of length one. Total natural indirect effect.
##'
##' @return A numeric vector of length two. Gradient of the proportion mediated with respect to pnde and tnie.
grad_prop_med_yreg_logistic <- Deriv::Deriv(prop_med_yreg_logistic)
## function (pnde, tnie)
## {
##     .e1 <- exp(pnde)
##     .e2 <- exp(tnie)
##     .e3 <- .e1 * .e2
##     .e4 <- .e3 - 1
##     .e5 <- .e2 - 1
##     c(pnde = (1 - .e3/.e4) * .e1 * .e5/.e4,
##       tnie = (1 - .e1 * .e5/.e4) * .e1 * .e2/.e4)
## }
## Use this vector to linearly combine Gamma_pnde and Gamma_tnie.
## Expanded
## c(pnde = (1 - (exp(pnde) * exp(tnie)) / (exp(pnde) * exp(tnie) - 1)) * exp(pnde) * (exp(tnie) - 1) / (exp(pnde) * exp(tnie) - 1),
##   tnie = (1 - exp(pnde) * (exp(tnie) - 1) / (exp(pnde) * exp(tnie) - 1)) * exp(pnde) * exp(tnie) / (exp(pnde) * exp(tnie) - 1))

## Symbolic differentiation by Deriv::Deriv
## https://github.com/sgsokol/Deriv
