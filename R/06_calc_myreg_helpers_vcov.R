################################################################################
### Helper functions for regression-based causal mediation analysis
##
## Created on: 2020-03-11
## Author: Kazuki Yoshida
################################################################################


###
### vcov extractors
################################################################################

Sigma_beta_hat <- function(mreg, mreg_fit, avar, cvar, emm_ac_mreg) {
    ## Assign row and column names, stratified by is.null(cvar), because paste0(avar, ":", cvar) will have avar:  .
    if(!is.null(cvar)){
      if(!is.null(emm_ac_mreg)){
        vcov_beta <- matrix(0, 2+2*length(cvar), 2+2*length(cvar))
        rownames(vcov_beta) <- colnames(vcov_beta) <- 
          c("(Intercept)", avar, cvar, paste0(avar, ":", cvar))
      } else {
        vcov_beta <- matrix(0, 2+length(cvar), 2+length(cvar))
        rownames(vcov_beta) <- colnames(vcov_beta) <- 
          c("(Intercept)", avar, cvar)
      }
    } else {
      vcov_beta <- matrix(0, 2, 2)
      rownames(vcov_beta) <- colnames(vcov_beta) <- 
        c("(Intercept)", avar)
    }
    
    # plug in non-zeros to corresponding elements:
    for(row_name in names(coef(mreg_fit))){
      for(col_name in names(coef(mreg_fit))){
        vcov_beta[row_name, col_name] <- vcov(mreg_fit)[row_name, col_name, drop = FALSE]
      }
    }
    
    return(vcov_beta)
    
}



# Note1: always has AxM row/column, even though there is interaction = FALSE
Sigma_theta_hat <- function(yreg, yreg_fit, avar, mvar, cvar, emm_ac_yreg, emm_mc_yreg, interaction) {
    ## Assign row and column names, stratified by is.null(cvar), because paste0(avar, ":", cvar) will generate string 'avar:'.
    if(!is.null(cvar)){
      if(!is.null(emm_ac_yreg) & !is.null(emm_mc_yreg)){
        vcov_theta <- matrix(0, 4+3*length(cvar), 4+3*length(cvar))
        rownames(vcov_theta) <- colnames(vcov_theta) <- c("(Intercept)", avar, mvar, paste0(avar,":", mvar), cvar, 
                                                          paste0(avar, ":", cvar), # AxC in Y model
                                                          paste0(mvar, ":", cvar)) # MxC
      } else if (!is.null(emm_ac_yreg) & is.null(emm_mc_yreg)){
        vcov_theta <- matrix(0, 4+2*length(cvar), 4+2*length(cvar))
        rownames(vcov_theta) <- colnames(vcov_theta) <- c("(Intercept)", avar, mvar, paste0(avar,":", mvar), cvar, 
                                                          paste0(avar, ":", cvar)) # AxC in Y model
      } else if (is.null(emm_ac_yreg) & !is.null(emm_mc_yreg)){
        vcov_theta <- matrix(0, 4+2*length(cvar), 4+2*length(cvar))
        rownames(vcov_theta) <- colnames(vcov_theta) <- c("(Intercept)", avar, mvar, paste0(avar,":", mvar), cvar, 
                                                          paste0(mvar, ":", cvar)) # MxC
      } else{
        vcov_theta <- matrix(0, 4+length(cvar), 4+length(cvar))
        rownames(vcov_theta) <- colnames(vcov_theta) <- c("(Intercept)", avar, mvar, paste0(avar,":", mvar), cvar) 
      }
     
    } else {
      vcov_theta <- matrix(0, 4, 4)
      rownames(vcov_theta) <- colnames(vcov_theta) <- c("(Intercept)", avar, mvar, paste0(avar, ":", mvar))
    }
    
    # plug in non-zeros to corresponding elements
    for(row_name in names(coef(yreg_fit))){
      for(col_name in names(coef(yreg_fit))){
        vcov_theta[row_name, col_name] <- vcov(yreg_fit)[row_name, col_name, drop = FALSE]
      }
    }
    
    return(vcov_theta)

    
}



Sigma_sigma_sq_hat <- function(mreg_fit) {

    ## VanderWeele 2015. p470 states
    ## 2 * (sigma_hat^2)^2 / (n-p)
    ##
    ## Chapter 2 Quadratic Forms of Random Variables
    ## http://pages.stat.wisc.edu/~st849-1/lectures/Ch02.pdf
    ## Corollary 6. In a full-rank Gaussian model for M with p covariates X.
    ## (||M - X beta_hat||^2 / sigma^2) ~ (central chi-squared DF = (n - p))
    ##
    ## https://en.wikipedia.org/wiki/Chi-squared_distribution
    ## Var(central chi-squared with DF = (n - p)) = 2 * (n - p)
    ##
    ## Var(sigma_hat^2) = Var(1/(n-p) * ||M - X beta_hat||^2)
    ##                  = 1/(n-p)^2 * Var(||M - X beta_hat||^2 * sigma^2/sigma^2)
    ##                  = 1/(n-p)^2 * (sigma^2)^2 * Var(||M - X beta_hat||^2 /sigma^2)
    ##                  = 1/(n-p)^2 * (sigma^2)^2 * (2 * (n-p))
    ##                  = 2 * (sigma^2)^2 / (n-p)
    ##
    matrix(2 * ((sigma(mreg_fit))^2)^2 / mreg_fit$df.residual)
}
