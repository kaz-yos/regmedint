/**/
/** */
/*  */
/* Created on 2020-03-08 */
/**/

/** Set libname */
libname w './';

/* Load SAS macro */
%include './mediation.sas';

/* Load data */
%include './data-valeri-vanderweele-2015.sas';

/* Invoke macro */
%mediation(
    data = data,
    yvar = y,
    avar = x,
    mvar = m,
    cvar = c,
    a0 = 0,
    a1 = 1,
    m = 1,
    yreg = survAFT_exp,
    mreg = logistic,
    interaction = true,
    casecontrol = false,
    output = full,
    c = 0.5,
    boot = ,
    cens = cens);
run;
