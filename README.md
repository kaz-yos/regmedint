
# regmedint <img src="man/figures/hex.png" align="right" height="140"/>

[![R-CMD-check](https://github.com/kaz-yos/regmedint/workflows/R-CMD-check/badge.svg)](https://github.com/kaz-yos/regmedint/actions)
[![](http://www.r-pkg.org/badges/version/regmedint)](http://www.r-pkg.org/pkg/regmedint)
[![CRAN RStudio mirror
downloads](http://cranlogs.r-pkg.org/badges/regmedint)](http://www.r-pkg.org/pkg/regmedint)

This is an extension of the regression-based causal mediation analysis
first proposed by Valeri and VanderWeele (2013) and Valeri and
VanderWeele (2015). The current version supports including effect
measure modification by covariates (treatment-covariate and
mediator-covariate product terms in mediator and outcome regression
models). It also accommodates the original ‘SAS’ macro (can be found at
Dr. VanderWeele’s [Tools and
Tutorials](https://www.hsph.harvard.edu/tyler-vanderweele/tools-and-tutorials/))
and PROC CAUSALMED procedure in ‘SAS’ when there is no effect measure
modification. Linear and logistic models are supported for the mediator
model. Linear, logistic, loglinear, Poisson, negative binomial, Cox, and
accelerated failure time (exponential and Weibull) models are supported
for the outcome model.

To cite this software, please use: *regmedint* (v1.0.0; Yoshida, Li, &
Mathur, 2021)
<!-- cite the r package:  https://easystats.github.io/report/articles/cite_packages.html -->

# Implemented models

The following grid of models are implemented. `yreg` refers to the
outcome model and `mreg` refers to the mediator model.

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">

| yreg \\\\ mreg       |             linear             |            logistic            |
|:---------------------|:------------------------------:|:------------------------------:|
| linear               |       :heavy_check_mark:       |       :heavy_check_mark:       |
| logistic<sup>1</sup> |       :heavy_check_mark:       |       :heavy_check_mark:       |
| loglinear            | :heavy_check_mark:<sup>2</sup> | :heavy_check_mark:<sup>2</sup> |
| poisson              |       :heavy_check_mark:       |       :heavy_check_mark:       |
| negbin               |       :heavy_check_mark:       |       :heavy_check_mark:       |
| survCox<sup>1</sup>  |       :heavy_check_mark:       |       :heavy_check_mark:       |
| survAFT exp          |       :heavy_check_mark:       |       :heavy_check_mark:       |
| survAFT weibull      |       :heavy_check_mark:       |       :heavy_check_mark:       |

<sup>1</sup> Approximation depends on the rare event assumptions.

<sup>2</sup> Implemented as a modified Poisson model (log link with
robust variance) as in Z2004.

See the corresponding vignettes (Articles on the package website) for
how to perform bootstrapping and multiple imputation.

# Installation

For the developmental version on Github, use the following commands to
install the package.

``` r
install.packages("devtools") # If you do not have devtools already.
```

    ## Error in contrib.url(repos, "source"): trying to use CRAN without setting a mirror

``` r
devtools::install_github("einsley1993/regmedint") 
```

``` r
install.packages("regmedint")
```

# Data Example

Here we show how to use `regmedint()` to analyze an empirical dataset.
We use an heart disease dataset from the UCI Machine Learning
Repository: <http://archive.ics.uci.edu/ml/datasets/Heart+Disease>. We
examine the effect of cholesterol level on the risk of heart disease.
The exposure is cholesterol level, the outcome is heart disease, the
mediator is blood pressure. Covariates include age, sex, chest pain,
blood sugar, and maximum heart rate.

## `regmedint()` to fit models

The `regmedint` function is the user interface for constructing a result
object of class `regmedint`. The interface is similar to the original
SAS macro. For survival outcomes, the indicator variable is an event
indicator (1 for event, 0 for censoring). `c_cond` vector is required be
specified. This vector is the vector of covariate values at which the
conditional effects are evaluated at.

``` r
library(regmedint)
# library(tidyverse)

## Example data: Heart disease dataset
# This data comes from the UCI Machine Learning Repository, containing a collection of demographic and clinical characteristics from 303 patients: http://archive.ics.uci.edu/ml/datasets/Heart+Disease

# install.packages("cheese")
require(cheese)
# dataset 'heart_disease':
# A: Cholesterol
# M: BP
# Y: HeartDisease

# Pre-process data to convert all variables to numeric type:
age <- heart_disease$Age
sex_M <- ifelse(heart_disease$Sex == "Male", 1, 0)
pain_typical <- ifelse(heart_disease$ChestPain == "Typical angina", 1, 0)
pain_atypical <- ifelse(heart_disease$ChestPain == "Atypical angina", 1, 0)
pain_non <- ifelse(heart_disease$ChestPain == "Non-anginal pain", 1, 0)
bp <- heart_disease$BP
cholesterol <- heart_disease$Cholesterol
bloodsugar_T <- ifelse(heart_disease$BloodSugar == "TRUE", 1, 0)
maximumHR <- heart_disease$MaximumHR
HD <- ifelse(heart_disease$HeartDisease == "Yes", 1, 0)

heart.disease <- cbind.data.frame(age, sex_M, pain_typical, pain_atypical, pain_non, 
                         bp, cholesterol, bloodsugar_T, maximumHR, HD)
```

1.  When there is no effect measure modification by covariates,
    `EMM_AC_Mmodel = NULL`, `EMM_AC_Ymodel = NULL`, `EMM_MC = NULL`.

``` r
regmedint_obj1 <- regmedint(data = heart.disease,
                            ## Variables
                            yvar = "HD",
                            avar = "cholesterol",
                            mvar = "bp",
                            cvar = c("age", "sex_M", "pain_non",
                                     "bloodsugar_T", "maximumHR"),
                            EMM_AC_Mmodel = NULL, 
                            EMM_AC_Ymodel = NULL, 
                            EMM_MC = NULL, 
                            eventvar = NULL,
                            ## Values at which effects are evaluated
                            a0 = 211,
                            a1 = 275,
                            m_cde = 131.7,
                            c_cond = c(mean(heart.disease$age), 0, 0, 
                                       mean(heart.disease$bloodsugar_T),
                                       mean(heart.disease$maximumHR)),
                            ## Model types
                            mreg = "linear",
                            yreg = "logistic",
                            ## Additional specification
                            interaction = TRUE,
                            casecontrol = FALSE)
summary(regmedint_obj1)
```

    ## ### Mediator model
    ## 
    ## Call:
    ## lm(formula = bp ~ cholesterol + age + sex_M + pain_non + bloodsugar_T + 
    ##     maximumHR, data = data)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -37.323 -11.853  -1.564  10.521  59.618 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  89.31699   11.86603   7.527 6.29e-13 ***
    ## cholesterol   0.02103    0.01942   1.083  0.27972    
    ## age           0.53800    0.12038   4.469 1.12e-05 ***
    ## sex_M        -1.35956    2.13286  -0.637  0.52433    
    ## pain_non     -2.44490    2.18569  -1.119  0.26422    
    ## bloodsugar_T  7.39882    2.74156   2.699  0.00736 ** 
    ## maximumHR     0.05626    0.04662   1.207  0.22848    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 16.72 on 296 degrees of freedom
    ## Multiple R-squared:  0.1152, Adjusted R-squared:  0.09726 
    ## F-statistic: 6.423 on 6 and 296 DF,  p-value: 0.000002236
    ## 
    ## ### Outcome model
    ## 
    ## Call:
    ## glm(formula = HD ~ cholesterol + bp + cholesterol:bp + age + 
    ##     sex_M + pain_non + bloodsugar_T + maximumHR, family = binomial(link = "logit"), 
    ##     data = data)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -2.4432  -0.7932  -0.2803   0.7902   2.1564  
    ## 
    ## Coefficients:
    ##                   Estimate  Std. Error z value     Pr(>|z|)    
    ## (Intercept)     1.57342995  5.71003005   0.276        0.783    
    ## cholesterol     0.00187241  0.02235911   0.084        0.933    
    ## bp              0.01122220  0.04414423   0.254        0.799    
    ## age             0.01744350  0.01849884   0.943        0.346    
    ## sex_M           1.78266199  0.34921575   5.105 0.0000003312 ***
    ## pain_non       -1.40164684  0.33883294  -4.137 0.0000352352 ***
    ## bloodsugar_T    0.04235784  0.39394071   0.108        0.914    
    ## maximumHR      -0.04507294  0.00806080  -5.592 0.0000000225 ***
    ## cholesterol:bp  0.00003755  0.00016984   0.221        0.825    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 417.98  on 302  degrees of freedom
    ## Residual deviance: 299.44  on 294  degrees of freedom
    ## AIC: 317.44
    ## 
    ## Number of Fisher Scoring iterations: 5
    ## 
    ## ### Mediation analysis 
    ##             est         se         Z          p       lower      upper
    ## cde  0.43635433 0.18584088 2.3479998 0.01887453  0.07211291 0.80059575
    ## pnde 0.45208929 0.20299243 2.2271239 0.02593899  0.05423144 0.84994714
    ## tnie 0.02900112 0.02967081 0.9774292 0.32835673 -0.02915261 0.08715484
    ## tnde 0.45532375 0.20935009 2.1749393 0.02963466  0.04500511 0.86564239
    ## pnie 0.02576666 0.02841159 0.9069065 0.36445626 -0.02991904 0.08145235
    ## te   0.48109041 0.20687280 2.3255372 0.02004325  0.07562716 0.88655365
    ## pm   0.07485022 0.07458535 1.0035512 0.31559497 -0.07133438 0.22103481
    ## 
    ## Evaluated at:
    ## avar: cholesterol
    ##  a1 (intervened value of avar) = 275
    ##  a0 (reference value of avar)  = 211
    ## mvar: bp
    ##  m_cde (intervend value of mvar for cde) = 131.7
    ## cvar: age sex_M pain_non bloodsugar_T maximumHR
    ##  c_cond (covariate vector value) = 54.43894 0 0 0.1485149 149.6073
    ## 
    ## Note that effect estimates can vary over m_cde and c_cond values when interaction = TRUE.

2.  When there is effect measure modification by covariates,
    `EMM_AC_Mmodel`, `EMM_AC_Ymodel` and `EMM_MC` can take a sub-vector
    of covariates in `cvar`.

``` r
regmedint_obj2 <- regmedint(data = heart.disease,
                            ## Variables
                            yvar = "HD",
                            avar = "cholesterol",
                            mvar = "bp",
                            cvar = c("age", "sex_M", "pain_non", "bloodsugar_T",
                                     "maximumHR"),
                            EMM_AC_Mmodel = c("age", "sex_M", "pain_non", "bloodsugar_T",
                                              "maximumHR"), 
                            EMM_AC_Ymodel = c("age", "sex_M", "pain_non", "maximumHR"), 
                            EMM_MC = c("age", "pain_non", "bloodsugar_T", "maximumHR"), 
                            eventvar = NULL,
                            ## Values at which effects are evaluated
                            a0 = 211,
                            a1 = 275,
                            m_cde = 131.7,
                            c_cond = c(mean(heart.disease$age), 0, 0, 
                                       mean(heart.disease$bloodsugar_T),
                                       mean(heart.disease$maximumHR)),
                            ## Model types
                            mreg = "linear",
                            yreg = "logistic",
                            ## Additional specification
                            interaction = TRUE,
                            casecontrol = FALSE)
summary(regmedint_obj2)
```

    ## ### Mediator model
    ## 
    ## Call:
    ## lm(formula = bp ~ cholesterol + age + sex_M + pain_non + bloodsugar_T + 
    ##     maximumHR + cholesterol:age + cholesterol:sex_M + cholesterol:pain_non + 
    ##     cholesterol:bloodsugar_T + cholesterol:maximumHR, data = data)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -37.558 -11.115  -1.517   9.361  59.999 
    ## 
    ## Coefficients:
    ##                            Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)              -13.541915  57.480892  -0.236   0.8139  
    ## cholesterol                0.437803   0.231494   1.891   0.0596 .
    ## age                        1.286279   0.571829   2.249   0.0252 *
    ## sex_M                      9.551806  10.756004   0.888   0.3753  
    ## pain_non                  11.731604  10.516827   1.116   0.2656  
    ## bloodsugar_T               9.131180  13.829918   0.660   0.5096  
    ## maximumHR                  0.387337   0.251201   1.542   0.1242  
    ## cholesterol:age           -0.003008   0.002324  -1.294   0.1965  
    ## cholesterol:sex_M         -0.044940   0.042359  -1.061   0.2896  
    ## cholesterol:pain_non      -0.056694   0.042179  -1.344   0.1799  
    ## cholesterol:bloodsugar_T  -0.009854   0.054340  -0.181   0.8562  
    ## cholesterol:maximumHR     -0.001345   0.001019  -1.319   0.1881  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 16.68 on 291 degrees of freedom
    ## Multiple R-squared:  0.1347, Adjusted R-squared:  0.102 
    ## F-statistic: 4.119 on 11 and 291 DF,  p-value: 0.000012
    ## 
    ## ### Outcome model
    ## 
    ## Call:
    ## glm(formula = HD ~ cholesterol + bp + cholesterol:bp + age + 
    ##     sex_M + pain_non + bloodsugar_T + maximumHR + cholesterol:age + 
    ##     cholesterol:sex_M + cholesterol:pain_non + cholesterol:maximumHR + 
    ##     bp:age + bp:pain_non + bp:bloodsugar_T + bp:maximumHR, family = binomial(link = "logit"), 
    ##     data = data)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -2.5347  -0.7863  -0.2692   0.8258   2.1475  
    ## 
    ## Coefficients:
    ##                          Estimate  Std. Error z value Pr(>|z|)
    ## (Intercept)           -3.27456006 14.15169818  -0.231    0.817
    ## cholesterol            0.00590456  0.04125558   0.143    0.886
    ## bp                     0.03823159  0.10784117   0.355    0.723
    ## age                    0.03842520  0.16834820   0.228    0.819
    ## sex_M                  1.42976060  1.69132559   0.845    0.398
    ## pain_non               0.85755051  3.11124366   0.276    0.783
    ## bloodsugar_T          -1.71590298  2.98998472  -0.574    0.566
    ## maximumHR             -0.01857166  0.06417856  -0.289    0.772
    ## cholesterol:bp         0.00004744  0.00020223   0.235    0.815
    ## cholesterol:age        0.00015077  0.00038859   0.388    0.698
    ## cholesterol:sex_M      0.00145389  0.00655231   0.222    0.824
    ## cholesterol:pain_non  -0.00792115  0.00679099  -1.166    0.243
    ## cholesterol:maximumHR -0.00008135  0.00017001  -0.478    0.632
    ## bp:age                -0.00042252  0.00115497  -0.366    0.714
    ## bp:pain_non           -0.00263279  0.02108868  -0.125    0.901
    ## bp:bloodsugar_T        0.01276240  0.02174713   0.587    0.557
    ## bp:maximumHR          -0.00005575  0.00045093  -0.124    0.902
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 417.98  on 302  degrees of freedom
    ## Residual deviance: 296.99  on 286  degrees of freedom
    ## AIC: 330.99
    ## 
    ## Number of Fisher Scoring iterations: 5
    ## 
    ## ### Mediation analysis 
    ##             est         se        Z         p       lower     upper
    ## cde  0.52417076 0.33182933 1.579640 0.1141894 -0.12620277 1.1745443
    ## pnde 0.53875628 0.33738505 1.596859 0.1102972 -0.12250627 1.2000188
    ## tnie 0.09977534 0.07255725 1.375126 0.1690925 -0.04243426 0.2419849
    ## tnde 0.55263268 0.34806254 1.587740 0.1123452 -0.12955737 1.2348227
    ## pnie 0.08589894 0.07311084 1.174914 0.2400292 -0.05739568 0.2291935
    ## te   0.63853162 0.34394634 1.856486 0.0633843 -0.03559082 1.3126541
    ## pm   0.20121358 0.14544718 1.383413 0.1665381 -0.08385765 0.4862848
    ## 
    ## Evaluated at:
    ## avar: cholesterol
    ##  a1 (intervened value of avar) = 275
    ##  a0 (reference value of avar)  = 211
    ## mvar: bp
    ##  m_cde (intervend value of mvar for cde) = 131.7
    ## cvar: age sex_M pain_non bloodsugar_T maximumHR
    ##  c_cond (covariate vector value) = 54.43894 0 0 0.1485149 149.6073
    ## 
    ## Note that effect estimates can vary over m_cde and c_cond values when interaction = TRUE.

## `print()` to examine simplified results

Implicit printing prints mreg, yreg, and mediation analysis point
estimates. All effect estimates are on the scale of the link function.

``` r
regmedint_obj2
```

    ## ### Mediator model
    ## 
    ## Call:
    ## lm(formula = bp ~ cholesterol + age + sex_M + pain_non + bloodsugar_T + 
    ##     maximumHR + cholesterol:age + cholesterol:sex_M + cholesterol:pain_non + 
    ##     cholesterol:bloodsugar_T + cholesterol:maximumHR, data = data)
    ## 
    ## Coefficients:
    ##              (Intercept)               cholesterol                       age                     sex_M  
    ##               -13.541915                  0.437803                  1.286279                  9.551806  
    ##                 pain_non              bloodsugar_T                 maximumHR           cholesterol:age  
    ##                11.731604                  9.131180                  0.387337                 -0.003008  
    ##        cholesterol:sex_M      cholesterol:pain_non  cholesterol:bloodsugar_T     cholesterol:maximumHR  
    ##                -0.044940                 -0.056694                 -0.009854                 -0.001345  
    ## 
    ## ### Outcome model
    ## 
    ## Call:  glm(formula = HD ~ cholesterol + bp + cholesterol:bp + age + 
    ##     sex_M + pain_non + bloodsugar_T + maximumHR + cholesterol:age + 
    ##     cholesterol:sex_M + cholesterol:pain_non + cholesterol:maximumHR + 
    ##     bp:age + bp:pain_non + bp:bloodsugar_T + bp:maximumHR, family = binomial(link = "logit"), 
    ##     data = data)
    ## 
    ## Coefficients:
    ##           (Intercept)            cholesterol                     bp                    age                  sex_M  
    ##           -3.27456006             0.00590456             0.03823159             0.03842520             1.42976060  
    ##              pain_non           bloodsugar_T              maximumHR         cholesterol:bp        cholesterol:age  
    ##            0.85755051            -1.71590298            -0.01857166             0.00004744             0.00015077  
    ##     cholesterol:sex_M   cholesterol:pain_non  cholesterol:maximumHR                 bp:age            bp:pain_non  
    ##            0.00145389            -0.00792115            -0.00008135            -0.00042252            -0.00263279  
    ##       bp:bloodsugar_T           bp:maximumHR  
    ##            0.01276240            -0.00005575  
    ## 
    ## Degrees of Freedom: 302 Total (i.e. Null);  286 Residual
    ## Null Deviance:       418 
    ## Residual Deviance: 297   AIC: 331
    ## ### Mediation analysis 
    ##        cde       pnde       tnie       tnde       pnie         te         pm 
    ## 0.52417076 0.53875628 0.09977534 0.55263268 0.08589894 0.63853162 0.20121358

## `summary()` to examine extended results

The `summary` method gives the summary for `mreg`, `yreg`, and mediation
analysis results. The `exponentiate` option will add the exponentiated
estimate and confidence limits if the outcome model is not a linear
model. The pure natural direct effect (`pnde`) is what is typically
called the natural direct effect (NDE). The total natural indirect
effect (`tnie`) is the corresponding natural indirect effect (NIE).

``` r
summary(regmedint_obj2, exponentiate = TRUE)
```

    ## ### Mediator model
    ## 
    ## Call:
    ## lm(formula = bp ~ cholesterol + age + sex_M + pain_non + bloodsugar_T + 
    ##     maximumHR + cholesterol:age + cholesterol:sex_M + cholesterol:pain_non + 
    ##     cholesterol:bloodsugar_T + cholesterol:maximumHR, data = data)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -37.558 -11.115  -1.517   9.361  59.999 
    ## 
    ## Coefficients:
    ##                            Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)              -13.541915  57.480892  -0.236   0.8139  
    ## cholesterol                0.437803   0.231494   1.891   0.0596 .
    ## age                        1.286279   0.571829   2.249   0.0252 *
    ## sex_M                      9.551806  10.756004   0.888   0.3753  
    ## pain_non                  11.731604  10.516827   1.116   0.2656  
    ## bloodsugar_T               9.131180  13.829918   0.660   0.5096  
    ## maximumHR                  0.387337   0.251201   1.542   0.1242  
    ## cholesterol:age           -0.003008   0.002324  -1.294   0.1965  
    ## cholesterol:sex_M         -0.044940   0.042359  -1.061   0.2896  
    ## cholesterol:pain_non      -0.056694   0.042179  -1.344   0.1799  
    ## cholesterol:bloodsugar_T  -0.009854   0.054340  -0.181   0.8562  
    ## cholesterol:maximumHR     -0.001345   0.001019  -1.319   0.1881  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 16.68 on 291 degrees of freedom
    ## Multiple R-squared:  0.1347, Adjusted R-squared:  0.102 
    ## F-statistic: 4.119 on 11 and 291 DF,  p-value: 0.000012
    ## 
    ## ### Outcome model
    ## 
    ## Call:
    ## glm(formula = HD ~ cholesterol + bp + cholesterol:bp + age + 
    ##     sex_M + pain_non + bloodsugar_T + maximumHR + cholesterol:age + 
    ##     cholesterol:sex_M + cholesterol:pain_non + cholesterol:maximumHR + 
    ##     bp:age + bp:pain_non + bp:bloodsugar_T + bp:maximumHR, family = binomial(link = "logit"), 
    ##     data = data)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -2.5347  -0.7863  -0.2692   0.8258   2.1475  
    ## 
    ## Coefficients:
    ##                          Estimate  Std. Error z value Pr(>|z|)
    ## (Intercept)           -3.27456006 14.15169818  -0.231    0.817
    ## cholesterol            0.00590456  0.04125558   0.143    0.886
    ## bp                     0.03823159  0.10784117   0.355    0.723
    ## age                    0.03842520  0.16834820   0.228    0.819
    ## sex_M                  1.42976060  1.69132559   0.845    0.398
    ## pain_non               0.85755051  3.11124366   0.276    0.783
    ## bloodsugar_T          -1.71590298  2.98998472  -0.574    0.566
    ## maximumHR             -0.01857166  0.06417856  -0.289    0.772
    ## cholesterol:bp         0.00004744  0.00020223   0.235    0.815
    ## cholesterol:age        0.00015077  0.00038859   0.388    0.698
    ## cholesterol:sex_M      0.00145389  0.00655231   0.222    0.824
    ## cholesterol:pain_non  -0.00792115  0.00679099  -1.166    0.243
    ## cholesterol:maximumHR -0.00008135  0.00017001  -0.478    0.632
    ## bp:age                -0.00042252  0.00115497  -0.366    0.714
    ## bp:pain_non           -0.00263279  0.02108868  -0.125    0.901
    ## bp:bloodsugar_T        0.01276240  0.02174713   0.587    0.557
    ## bp:maximumHR          -0.00005575  0.00045093  -0.124    0.902
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 417.98  on 302  degrees of freedom
    ## Residual deviance: 296.99  on 286  degrees of freedom
    ## AIC: 330.99
    ## 
    ## Number of Fisher Scoring iterations: 5
    ## 
    ## ### Mediation analysis 
    ##             est         se        Z         p       lower     upper exp(est) exp(lower) exp(upper)
    ## cde  0.52417076 0.33182933 1.579640 0.1141894 -0.12620277 1.1745443 1.689058  0.8814361   3.236668
    ## pnde 0.53875628 0.33738505 1.596859 0.1102972 -0.12250627 1.2000188 1.713874  0.8847004   3.320179
    ## tnie 0.09977534 0.07255725 1.375126 0.1690925 -0.04243426 0.2419849 1.104923  0.9584535   1.273775
    ## tnde 0.55263268 0.34806254 1.587740 0.1123452 -0.12955737 1.2348227 1.737822  0.8784842   3.437769
    ## pnie 0.08589894 0.07311084 1.174914 0.2400292 -0.05739568 0.2291935 1.089696  0.9442204   1.257585
    ## te   0.63853162 0.34394634 1.856486 0.0633843 -0.03559082 1.3126541 1.893698  0.9650351   3.716023
    ## pm   0.20121358 0.14544718 1.383413 0.1665381 -0.08385765 0.4862848       NA         NA         NA
    ## 
    ## Evaluated at:
    ## avar: cholesterol
    ##  a1 (intervened value of avar) = 275
    ##  a0 (reference value of avar)  = 211
    ## mvar: bp
    ##  m_cde (intervend value of mvar for cde) = 131.7
    ## cvar: age sex_M pain_non bloodsugar_T maximumHR
    ##  c_cond (covariate vector value) = 54.43894 0 0 0.1485149 149.6073
    ## 
    ## Note that effect estimates can vary over m_cde and c_cond values when interaction = TRUE.

Use `coef` to extract the mediation analysis results only.

``` r
coef(summary(regmedint_obj2, exponentiate = TRUE))
```

    ##             est         se        Z         p       lower     upper exp(est) exp(lower) exp(upper)
    ## cde  0.52417076 0.33182933 1.579640 0.1141894 -0.12620277 1.1745443 1.689058  0.8814361   3.236668
    ## pnde 0.53875628 0.33738505 1.596859 0.1102972 -0.12250627 1.2000188 1.713874  0.8847004   3.320179
    ## tnie 0.09977534 0.07255725 1.375126 0.1690925 -0.04243426 0.2419849 1.104923  0.9584535   1.273775
    ## tnde 0.55263268 0.34806254 1.587740 0.1123452 -0.12955737 1.2348227 1.737822  0.8784842   3.437769
    ## pnie 0.08589894 0.07311084 1.174914 0.2400292 -0.05739568 0.2291935 1.089696  0.9442204   1.257585
    ## te   0.63853162 0.34394634 1.856486 0.0633843 -0.03559082 1.3126541 1.893698  0.9650351   3.716023
    ## pm   0.20121358 0.14544718 1.383413 0.1665381 -0.08385765 0.4862848       NA         NA         NA

Note that the estimates can be re-evaluated at different `m_cde` and
`c_cond` without re-fitting the underlyng models.

``` r
coef(summary(regmedint_obj2, exponentiate = TRUE, 
             m_cde = min(heart.disease$bp), 
             c_cond = c(min(heart.disease$age), 0, 1, 
                                       max(heart.disease$bloodsugar_T),
                                       min(heart.disease$maximumHR))))
```

    ##             est        se          Z         p     lower    upper exp(est) exp(lower) exp(upper)
    ## cde  0.06653441 1.3668539 0.04867705 0.9611767 -2.612450 2.745519 1.068798 0.07335460  15.572693
    ## pnde 0.15779864 1.3416660 0.11761395 0.9063736 -2.471818 2.787416 1.170930 0.08443119  16.239000
    ## tnie 0.54539895 0.8725055 0.62509517 0.5319086 -1.164680 2.255478 1.725297 0.31202238   9.539855
    ## tnde 0.19443891 1.3683849 0.14209373 0.8870060 -2.487546 2.876424 1.214629 0.08311366  17.750683
    ## pnie 0.50875868 0.8269748 0.61520461 0.5384196 -1.112082 2.129599 1.663225 0.32887352   8.411497
    ## te   0.70319759 1.4715535 0.47786071 0.6327493 -2.180994 3.587389 2.020202 0.11292920  36.139607
    ## pm   0.83245440 1.1787537 0.70621572 0.4800540 -1.477860 3.142769       NA         NA         NA

# Formulas

See
<https://github.com/kaz-yos/regmedint-supplement/blob/master/supplement.pdf>
for the following formulas.

## Effect formulas in the supplementary document

| yreg \\\\ mreg  | linear               | logistic             |
|-----------------|----------------------|----------------------|
| linear          | Formulas (1) - (5)   | Formulas (11) - (15) |
|                 |                      |                      |
| logistic        | Formulas (21) - (25) | Formulas (31) - (35) |
| loglinear       | Formulas (21) - (25) | Formulas (31) - (35) |
| poisson         | Formulas (21) - (25) | Formulas (31) - (35) |
| negbin          | Formulas (21) - (25) | Formulas (31) - (35) |
|                 |                      |                      |
| survCox         | Formulas (21) - (25) | Formulas (31) - (35) |
| survAFT exp     | Formulas (21) - (25) | Formulas (31) - (35) |
| survAFT weibull | Formulas (21) - (25) | Formulas (31) - (35) |

## Standard error formulas in the supplementary document

| yreg \\\\ mreg  | linear               | logistic             |
|-----------------|----------------------|----------------------|
| linear          | Formulas (6) - (10)  | Formulas (16) - (20) |
|                 |                      |                      |
| logistic        | Formulas (26) - (30) | Formulas (36) - (40) |
| loglinear       | Formulas (26) - (30) | Formulas (36) - (40) |
| poisson         | Formulas (26) - (30) | Formulas (36) - (40) |
| negbin          | Formulas (26) - (30) | Formulas (36) - (40) |
|                 |                      |                      |
| survCox         | Formulas (26) - (30) | Formulas (36) - (40) |
| survAFT exp     | Formulas (26) - (30) | Formulas (36) - (40) |
| survAFT weibull | Formulas (26) - (30) | Formulas (36) - (40) |

Note: The point estimate and standard error formulas (multivariate delta
method) were derived based on the following references.

-   V2015: VanderWeele (2015) Explanation in Causal Inference.
-   VV2013A: Valeri & VanderWeele (2013) Appendix
-   VV2015A: Valeri & VanderWeele (2015) Appendix

## Effect formulas are based on the following propositions

| yreg \\\\ mreg  | linear                               | logistic                             |
|-----------------|--------------------------------------|--------------------------------------|
| linear          | V2015 p466 Proposition 2.3           | V2015 p471 Proposition 2.5           |
|                 |                                      |                                      |
| logistic        | V2015 p468 Proposition 2.4           | V2015 p473 Proposition 2.6           |
| loglinear       | VV2013A p8 Use Proposition 2.4       | VV2013A p8 Use Proposition 2.6       |
| poisson         | VV2013A p8 Use Proposition 2.4       | VV2013A p8 Use Proposition 2.6       |
| negbin          | VV2013A p8 Use Proposition 2.4       | VV2013A p8 Use Proposition 2.6       |
|                 |                                      |                                      |
| survCox         | V2015 p496 Proposition 4.4 (Use 2.4) | V2015 p499 Proposition 4.6 (Use 2.6) |
| survAFT exp     | V2015 p494 Proposition 4.1 (Use 2.4) | V2015 p495 Proposition 4.3 (Use 2.6) |
| survAFT weibull | V2015 p494 Proposition 4.1 (Use 2.4) | V2015 p495 Proposition 4.3 (Use 2.6) |

## Standard error formulas are based on the following propositions

| yreg \\\\ mreg  | linear                         | logistic                       |
|-----------------|--------------------------------|--------------------------------|
| linear          | V2015 p466 Proposition 2.3     | V2015 p471 Proposition 2.5     |
|                 |                                |                                |
| logistic        | V2015 p468 Proposition 2.4     | V2015 p473 Proposition 2.6     |
| loglinear       | VV2013A p8 Use Proposition 2.4 | VV2013A p8 Use Proposition 2.6 |
| poisson         | VV2013A p8 Use Proposition 2.4 | VV2013A p8 Use Proposition 2.6 |
| negbin          | VV2013A p8 Use Proposition 2.4 | VV2013A p8 Use Proposition 2.6 |
|                 |                                |                                |
| survCox         | V2015 p496 Use Proposition 2.4 | V2015 p499 Use Proposition 2.6 |
| survAFT exp     | V2015 p494 Use Proposition 2.4 | V2015 p495 Use Proposition 2.6 |
| survAFT weibull | V2015 p494 Use Proposition 2.4 | V2015 p495 Use Proposition 2.6 |

<!-- # Design -->
<!-- The software design is outlined here for those who may be interested. -->
<!-- - Call structure -->
<!--   - regmedint UI function -->
<!--     - new_regmedint internal constructor -->
<!--       - fit_mreg -->
<!--       - fit_yreg -->
<!--       - calc_myreg calls a specialized worker function, which return two functions, one for point estimates and the other for standard error estimate. -->
<!--         - calc_myreg_mreg_linear_yreg_linear -->
<!--         - calc_myreg_mreg_linear_yreg_logistic -->
<!--         - calc_myreg_mreg_logistic_yreg_linear -->
<!--         - calc_myreg_mreg_logistic_yreg_logistic -->
<!-- - regmedint object structure -->
<!--   - mreg_fit mediator regression model object as is -->
<!--   - yreg_fit outcome regression model object as is -->
<!--   - myreg_funs list -->
<!--     - est_fun: (a0,a1,m_cde,c_cond) → (cde,pnde,tnie,tnde,pnie,te,pm) -->
<!--     - se_fun: (a0,a1,m_cde,c_cond) → se for (cde,pnde,tnie,tnde,pnie,te,pm) -->
<!--     - args preserves arguments given to the UI -->
<!-- - User methods for the regmedint object -->
<!--   - print.regmedint: prints coefficients for mreg, yreg, and mediation analysis -->
<!--   - summary.regmedint: regmedint → summary_regmedint -->
<!--     - print.summary_regmedint: prints summary objects for mreg, yreg, and mediation analysis -->
<!--     - coef.summary_regmedint: -->
<!--   - coef.regmedint: regmedint → vector (cde,pnde,tnie,tnde,pnie,te,pm) -->
<!--   - vcov.regmedint: regmedint → matrix (cde,pnde,tnie,tnde,pnie,te,pm). Off-diagonals are NA. -->
<!--   - confint.regmedint: regmedint → matrix of (lower,upper) -->

# Similar or related projects for counterfactual-based causal mediation analysis

## R

-   mediation (simulation-based):
    <https://CRAN.R-project.org/package=mediation>
-   medflex (natural effect model):
    <https://CRAN.R-project.org/package=medflex>
-   intmed (interventional analogue):
    <https://CRAN.R-project.org/package=intmed>
-   CMAverse (regression-based approach, weighting-based approach,
    inverse odd-ratio weighting, natural effect model, marginal
    structural model, g-formula approach):
    <https://bs1125.github.io/CMAverse/>
-   mediator (regression-based): <https://github.com/GerkeLab/mediator>
-   causalMediation (regression-based):
    <https:/github.com/harvard-P01/causalMediation>

## Other statistical environment

-   SAS macro (original regression-based)
    <https://www.hsph.harvard.edu/tyler-vanderweele/tools-and-tutorials/>
-   SAS PROC CAUSALMED (regression-based)
    <https://support.sas.com/rnd/app/stat/procedures/causalmed.html>

# References

-   V2015: VanderWeele (2015) Explanation in Causal Inference.
-   VV2013: Valeri & VanderWeele (2013) Psych Method. 18:137.
-   VV2015: Valeri & VanderWeele (2015) Epidemiology. 26:e23.
-   Z2004: Zou (2004) Am J Epidemiol 159:702.
