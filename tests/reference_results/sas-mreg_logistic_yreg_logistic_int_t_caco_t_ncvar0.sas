
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
    mvar = bili_bin,
    cvar = ,
    a0 = 1.1,
    a1 = 2.3,
    m = 1.4,
    yreg = logistic,
    mreg = logistic,
    interaction = true,
    casecontrol = true,
    output = full,
    c = ,
    boot = ,
    cens = cens);
run;
