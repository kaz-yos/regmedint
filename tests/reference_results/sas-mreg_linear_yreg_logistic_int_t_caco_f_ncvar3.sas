
/** Set libname */
libname w './';

/* Load SAS macro */
%include './mediation.sas';

/* Load data */
proc import datafile = './data-pbc_cc.csv'
    out = pbc_cc
    dbms = csv
    replace;
run;

%mediation(
    data = pbc_cc,
    yvar = spiders,
    avar = trt,
    mvar = bili,
    cvar = age male stage,
    a0 = 1,
    a1 = 2,
    m = 1,
    yreg = logistic,
    mreg = linear,
    interaction = true,
    casecontrol = false,
    output = full,
    c = 50 1 2,
    boot = ,
    cens = cens);
run;
