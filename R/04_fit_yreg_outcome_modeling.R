################################################################################
### Internal helper functions
##
## Created on: 2020-03-09
## Author: Kazuki Yoshida
################################################################################


###
### Model fitters for yreg
################################################################################

## The third and subsequent paragraphs go into details.
## http://r-pkgs.had.co.nz/man.html#roxygen-comments

##' Fit a model for the outcome given the treatment, mediator, and covariates.
##'
##' The outcome model type \code{yreg} can be one of the following \code{"linear"}, \code{"logistic"}, \code{"loglinear"} (implemented as modified Poisson), \code{"poisson"}, \code{"negbin"}, \code{"survCox"}, \code{"survAFT_exp"}, or \code{"survAFT_weibull"}.
##'
##' The outcome regression functions to be called are the following:
##' \itemize{
##'   \item \code{"linear"} \code{\link{lm}}
##'   \item \code{"logistic"} \code{\link{glm}}
##'   \item \code{"loglinear"} \code{\link{glm}} (modified Poisson)
##'   \item \code{"poisson"} \code{\link{glm}}
##'   \item \code{"negbin"} \code{\link[MASS]{glm.nb}}
##'   \item \code{"survCox"} \code{\link[survival]{coxph}}
##'   \item \code{"survAFT_exp"} \code{\link[survival]{survreg}}
##'   \item \code{"survAFT_weibull"} \code{\link[survival]{survreg}}
##' }
##'
##' @inheritParams regmedint
##'
##' @return Model fit object from on of the above regression functions.
fit_yreg <- function(yreg,
                     data,
                     yvar,
                     avar,
                     mvar,
                     cvar,
                     emm_ac_yreg = NULL,
                     emm_mc_yreg = NULL,
                     eventvar,
                     interaction) {

    ## Create a string representation of the formula
    string_formula <- string_yreg_formula(yvar,
                                          avar,
                                          mvar,
                                          cvar,
                                          emm_ac_yreg,
                                          emm_mc_yreg,
                                          interaction,
                                          eventvar)

    ## Quasi-quoting to make the formula readable.
    ## bquote suppresses evaluation except within .(...).
    ## Evaluate restart the evaluation with the .() part
    ## already expanded.
    if (yreg == "linear") {

        eval(
            bquote(
                lm(formula = .(as.formula(string_formula)),
                   data = data)
            )
        )

    } else if (yreg == "logistic") {

        eval(
            bquote(
                glm(formula = .(as.formula(string_formula)),
                    family = binomial(link = "logit"),
                    data = data)
            )
        )

    } else if (yreg == "loglinear") {

        message("loglinear is implemented as modified Poisson (Zou 2004).")
        yreg_fit <-
            eval(
                bquote(
                    glm(formula = .(as.formula(string_formula)),
                        family = poisson(link = "log"),
                        data = data)
                )
            )
        class(yreg_fit) <- c("regmedint_mod_poisson", class(yreg_fit))
        return(yreg_fit)

    } else if (yreg == "poisson") {

        eval(
            bquote(
                glm(formula = .(as.formula(string_formula)),
                    family = poisson(link = "log"),
                    data = data)
            )
        )

    } else if (yreg == "negbin") {

        eval(
            bquote(
                MASS::glm.nb(formula = .(as.formula(string_formula)),
                             data = data)
            )
        )

    } else if (yreg == "survCox") {

        eval(
            bquote(
                survival::coxph(formula = .(as.formula(string_formula)),
                                data = data,
                                ties = "efron")
            )
        )

    } else if (yreg == "survAFT_exp") {

        eval(
            bquote(
                survival::survreg(formula = .(as.formula(string_formula)),
                                  data = data,
                                  dist = "exponential")
            )
        )

    } else if (yreg == "survAFT_weibull") {

        eval(
            bquote(
                survival::survreg(formula = .(as.formula(string_formula)),
                                  data = data,
                                  dist = "weibull")
            )
        )

    } else {

        stop("Unsupported model type in yreg")

    }
}


###
### Formula string creators
################################################################################

string_yreg_formula <- function(yvar,
                                avar,
                                mvar,
                                cvar,
                                emm_ac_yreg,
                                emm_mc_yreg,
                                interaction,
                                eventvar) {

    assertthat::assert_that(!is.null(mvar))
    assertthat::assert_that(!is.null(avar))
    assertthat::assert_that(!is.null(yvar))

    ## Create A * M or A + M depending on interaction.
    if (interaction) {
      amvar_string <- paste(avar, mvar, paste0(avar, ":", mvar), sep = " + ")
    } else {
      amvar_string <- paste(avar, mvar, sep = " + ")
    }

    ## Add covariates if they exist.
    if (is.null(cvar)) {
      amcvar_string <- amvar_string
    } else {
      if(is.null(emm_ac_yreg) & is.null(emm_mc_yreg)){
        cvar_string <- paste0(cvar, collapse = " + ")
        amcvar_string <- sprintf("%s + %s", amvar_string, cvar_string)
      }
      else if(!is.null(emm_ac_yreg) & is.null(emm_mc_yreg)){
        cvar_string <- paste0(cvar, collapse = " + ")
        ## Add avar * emm_ac_yreg product terms:
        temp <- NULL
        for(i in 1:length(emm_ac_yreg)){
          temp <- paste0(c(temp, paste0(c(avar, emm_ac_yreg[i]), collapse = " : ")), collapse = " + ")
        }
        amcvar_string <- paste(amvar_string, cvar_string, temp, sep = " + ") 
      }
      else if(is.null(emm_ac_yreg) & !is.null(emm_mc_yreg)){
        cvar_string <- paste0(cvar, collapse = " + ")
        ## Add mvar * emm_mc_yreg product terms:
        temp <- NULL
        for(i in 1:length(emm_mc_yreg)){
          temp <- paste0(c(temp, paste0(c(mvar, emm_mc_yreg[i]), collapse = " : ")), collapse = " + ")
        }
        amcvar_string <- paste(amvar_string, cvar_string, temp, sep = " + ")  
      }
      else if(!is.null(emm_ac_yreg) & !is.null(emm_mc_yreg)){
        cvar_string <- paste0(cvar, collapse = " + ")
        ## Add avar * emm_ac_mreg product terms:
        temp1 <- temp2 <- NULL
        for(i in 1:length(emm_ac_yreg)){
          temp1 <- paste0(c(temp1, paste0(c(avar, emm_ac_yreg[i]), collapse = " : ")), collapse = " + ")
        }
        ## Add mvar * emm_mc_yreg product terms:
        for(j in 1:length(emm_mc_yreg)){
          temp2 <- paste0(c(temp2, paste0(c(mvar, emm_mc_yreg[j]), collapse = " : ")), collapse = " + ")
        }
        amcvar_string <- paste(amvar_string, cvar_string, temp1, temp2, sep = " + ") 
      }
    }
    
    ## eventvar must be NULL for a non-survival outcome model.
    if (is.null(eventvar)) {
      return(sprintf("%s ~ %s", yvar, amcvar_string))
      
    } else {
      ## Survival outcome
      surv_string <- sprintf("Surv(%s, %s)", yvar, eventvar)
      return(sprintf("%s ~ %s", surv_string, amcvar_string))
      
    }
}


##' Robust sandwich variance estimator for modified Poisson
##'
##' Provide robust sandwich variance-covariance estimate using \code{\link[sandwich]{sandwich}}.
##'
##' @param object A model object of the class \code{regmedint_mod_poisson}
##' @param ... For compatibility with the generic.
##'
##' @return A variance-covariance matrix using the \code{\link[sandwich]{sandwich}}.
vcov.regmedint_mod_poisson <- function(object, ...) {
    ## Drop the regmedint_mod_poisson class from the object because
    ## sandwich::sandwich internally invoke a summary method.
    ## We want the regular summary.glm to be invoked, not our version
    ## summary.regmedint_mod_poisson.
    class(object) <- setdiff(class(object), "regmedint_mod_poisson")
    sandwich::sandwich(object, ...)
}


##' Summary with robust sandwich variance estimator for modified Poisson
##'
##' This is a version of \code{\link{summary.glm}} modified to use the robust variance estimator \code{\link[sandwich]{sandwich}}.
##'
##' @param object A model object of the class \code{regmedint_mod_poisson}
##' @param ... For compatibility with the generic.
##'
##' @return An object of the class \code{summary.glm}
summary.regmedint_mod_poisson <- function(object, ...) {
    ## Set items that are irrelevant and not included in the arguments.
    dispersion <- NULL
    correlation <- FALSE
    symbolic.cor <- FALSE
    ## https://github.com/wch/r-source/blob/trunk/src/library/stats/R/lm.R
    ## using qr(<lm>)  as interface to  <lm>$qr :
    qr.lm <- function(x, ...) {
        if(is.null(r <- x$qr))
            stop("lm object does not have a proper 'qr' component.
 Rank zero or should not have used lm(.., qr=FALSE).")
        r
    }

    ## https://github.com/wch/r-source/blob/trunk/src/library/stats/R/glm.R
    est.disp <- FALSE
    df.r <- object$df.residual
    if(is.null(dispersion))	# calculate dispersion if needed
	dispersion <-
	    if(object$family$family %in% c("poisson", "binomial"))  1
	    else if(df.r > 0) {
                est.disp <- TRUE
		if(any(object$weights==0))
		    warning("observations with zero weight not used for calculating dispersion")
		sum((object$weights*object$residuals^2)[object$weights > 0])/ df.r
	    } else {
                est.disp <- TRUE
                NaN
            }

    ## calculate scaled and unscaled covariance matrix

    aliased <- is.na(coef(object))  # used in print method
    p <- object$rank
    if (p > 0) {
        p1 <- 1L:p
	Qr <- qr.lm(object)
        ## WATCHIT! doesn't this rely on pivoting not permuting 1L:p? -- that's guaranteed
        coef.p <- object$coefficients[Qr$pivot[p1]]
        ## Changed from below.
        ## covmat.unscaled <- chol2inv(Qr$qr[p1,p1,drop=FALSE])
        ## The dispersion is fixed at 1 for poisson, so vcov is ok.
        ## This should dispatch a robust variance estimator.
        covmat.unscaled <- vcov(object)
        dimnames(covmat.unscaled) <- list(names(coef.p),names(coef.p))
        covmat <- dispersion*covmat.unscaled
        var.cf <- diag(covmat)

        ## calculate coef table

        s.err <- sqrt(var.cf)
        tvalue <- coef.p/s.err

        dn <- c("Estimate", "Std. Error")
        if(!est.disp) { # known dispersion
            pvalue <- 2*pnorm(-abs(tvalue))
            coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
            dimnames(coef.table) <- list(names(coef.p),
                                         c(dn, "z value","Pr(>|z|)"))
        } else if(df.r > 0) {
            pvalue <- 2*pt(-abs(tvalue), df.r)
            coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
            dimnames(coef.table) <- list(names(coef.p),
                                         c(dn, "t value","Pr(>|t|)"))
        } else { # df.r == 0
            coef.table <- cbind(coef.p, NaN, NaN, NaN)
            dimnames(coef.table) <- list(names(coef.p),
                                         c(dn, "t value","Pr(>|t|)"))
        }
        df.f <- NCOL(Qr$qr)
    } else {
        coef.table <- matrix(, 0L, 4L)
        dimnames(coef.table) <-
            list(NULL, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
        covmat.unscaled <- covmat <- matrix(, 0L, 0L)
        df.f <- length(aliased)
    }
    ## return answer

    ## these need not all exist, e.g. na.action.
    keep <- match(c("call","terms","family","deviance", "aic",
                    "contrasts", "df.residual","null.deviance","df.null",
                    "iter", "na.action"), names(object), 0L)
    ans <- c(object[keep],
	     list(deviance.resid = residuals(object, type = "deviance"),
		  coefficients = coef.table,
                  aliased = aliased,
		  dispersion = dispersion,
		  df = c(object$rank, df.r, df.f),
		  cov.unscaled = covmat.unscaled,
		  cov.scaled = covmat))

    if(correlation && p > 0) {
	dd <- sqrt(diag(covmat.unscaled))
	ans$correlation <-
	    covmat.unscaled/outer(dd,dd)
	ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.glm"
    return(ans)
}
