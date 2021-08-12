################################################################################
### Helper functions for regression-based causal mediation analysis
##
## Created on: 2020-03-11
## Author: Kazuki Yoshida
################################################################################


###
### vcov extractors
################################################################################

Sigma_beta_hat <- function(mreg, mreg_fit, avar, cvar, EMM_AC_Mmodel) {
  ## Initialize with 0: when computing SE later, t(Gamma) %*% Sigma %*% Gamma won't matter if some C or AC are missing
  ## Assign row and column names, stratified by is.null(cvar), because paste0(avar, ":", cvar) will have avar:  .
  ## Be careful to stratify is.null(EMM_AC_Mmodel), becaue vcov_beta won't have AxC columns if is.null(EMM_AC_Mmodel)
  if(!is.null(cvar)){
    if(!is.null(EMM_AC_Mmodel)){
      vcov_beta <- matrix(0, 2+2*length(cvar), 2+2*length(cvar))
      rownames(vcov_beta) <- colnames(vcov_beta) <- 
        c("(Intercept)", avar, cvar, paste0(avar, ":", cvar))
    } else{
      vcov_beta <- matrix(0, 2+length(cvar), 2+length(cvar))
      rownames(vcov_beta) <- colnames(vcov_beta) <- 
        c("(Intercept)", avar, cvar)
    }
  } else {
    rownames(vcov_beta) <- colnames(vcov_beta) <- 
      c("(Intercept)", avar)
  }
  
  # plug in non-zeros to corresponding elements:
  for(row_name in names(coef(mreg_fit))){
    for(col_name in names(coef(mreg_fit))){
      vcov_beta[row_name, col_name] <- vcov(mreg_fit)[row_name, col_name]
    }
  }
   return(vcov_beta) # as long as beta2 or beta3 is not null, don't delete any columns
  
  # Trick: drop zero rows and columns
  # vcov_beta[rowSums(!as.matrix(vcov_beta)) < ncol(vcov_beta), ]
  # return(vcov_beta[rowSums(!as.matrix(vcov_beta)) < ncol(vcov_beta), rowSums(!as.matrix(vcov_beta)) < ncol(vcov_beta)])
  
}




Sigma_theta_hat <- function(yreg, yreg_fit, avar, mvar, cvar, EMM_AC_Ymodel, EMM_MC, interaction) {
  
  # vcov_raw <- vcov(yreg_fit)
  
  ## Older versions of survival:::vcov.survreg() did not give dimension names.
  ## https://github.com/therneau/survival/commit/ed1c71b3817d4bfced43ed374e5e598e5f229bb8
  
  # if (yreg == "survAFT_exp") {
  # 
  #     vcov_ready <- vcov_raw
  #     if (is.null(dimnames(vcov_ready))) {
  #         ## Older vcov.survreg gives an unnamed vcov matrix
  #         dimnames(vcov_ready) <- list(names(coef(yreg_fit)),
  #                                      names(coef(yreg_fit)))
  #     }
  # 
  # } else if (yreg == "survAFT_weibull") {
  # 
  #     ## vcov.survreg(weibull_fit) has an extra row and column corresponding
  #     ## to the log(scale) parameter (See the above commit).
  #     coef_ind <- seq_along(coef(yreg_fit))
  #     vcov_ready <- vcov_raw[coef_ind, coef_ind]
  #     if (is.null(dimnames(vcov_ready))) {
  #         ## Older vcov.survreg gives an unnamed vcov matrix
  #         dimnames(vcov_ready) <- list(names(coef(yreg_fit)),
  #                                      names(coef(yreg_fit)))
  #     }
  # 
  # } else if (yreg == "survCox") {
  # 
  #     ## Pad the left and upper edges with zeros by creating a block diagonal.
  #     vcov_ready <- Matrix::bdiag(matrix(0),
  #                                 vcov_raw)
  #     vars_cox <- c("(Intercept)", names(coef(yreg_fit)))
  #     dimnames(vcov_ready) <- list(vars_cox,
  #                                  vars_cox)
  # 
  # } else {
  # 
  #     vcov_ready <- vcov_raw
  # 
  # }
  
  
  ## Assign row and column names, stratified by is.null(cvar), because paste0(avar, ":", cvar) will have the extra colon after avar.
  if(!is.null(cvar)){
    if(!is.null(EMM_AC_Ymodel) & !is.null(EMM_MC)){
      vcov_theta <- matrix(0, 4+3*length(cvar), 4+3*length(cvar))
      rownames(vcov_theta) <- colnames(vcov_theta) <- c("(Intercept)", avar, mvar, paste0(avar,":", mvar), cvar, 
                                                        paste0(avar, ":", cvar), # AxC in Y model
                                                        paste0(mvar, ":", cvar)) # MxC
    } else if (is.null(EMM_AC_Ymodel) & is.null(EMM_MC)){
      vcov_theta <- matrix(0, 4+length(cvar), 4+length(cvar))
      rownames(vcov_theta) <- colnames(vcov_theta) <- c("(Intercept)", avar, mvar, paste0(avar,":", mvar), cvar) 
    } else if (!is.null(EMM_AC_Ymodel) & is.null(EMM_MC)){
      vcov_theta <- matrix(0, 4+2*length(cvar), 4+2*length(cvar))
      rownames(vcov_theta) <- colnames(vcov_theta) <- c("(Intercept)", avar, mvar, paste0(avar,":", mvar), cvar,
                                                        paste0(avar, ":", cvar)) 
    } else if (is.null(EMM_AC_Ymodel) & !is.null(EMM_MC)){
      vcov_theta <- matrix(0, 4+2*length(cvar), 4+2*length(cvar))
      rownames(vcov_theta) <- colnames(vcov_theta) <- c("(Intercept)", avar, mvar, paste0(avar,":", mvar), cvar,
                                                        paste0(mvar, ":", cvar)) 
      }
    } else {
    rownames(vcov_theta) <- colnames(vcov_theta) <- c("(Intercept)", avar, mvar, paste0(avar,":", mvar))
  }
  
  # plug in non-zeros to corresponding elements:
  for(row_name in names(coef(yreg_fit))){
    for(col_name in names(coef(yreg_fit))){
      vcov_theta[row_name, col_name] <- vcov(yreg_fit)[row_name, col_name]
    }
  }
  
  return(vcov_theta) # as long as theta4-6 are not null, don't delete any columns
  # return(vcov_theta[rowSums(!as.matrix(vcov_theta)) < ncol(vcov_theta), rowSums(!as.matrix(vcov_theta)) < ncol(vcov_theta)])
  
  
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
