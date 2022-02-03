/* Downloaded from https://www.hsph.harvard.edu/tyler-vanderweele/tools-and-tutorials/ */
/* on 2020-03-08 */
/* https://cdn1.sph.harvard.edu/wp-content/uploads/sites/603/2019/03/MediationPsychMethods.zip */

/* Valeri L, VanderWeele TJ.
Epidemiology. 2015 Mar;26(2):e23-4. doi: 10.1097/EDE.0000000000000253.
    SAS macro for causal mediation analysis with survival data. */

/* Valeri L, Vanderweele TJ.
Psychol Methods. 2013 Jun;18(2):137-50. doi: 10.1037/a0031034. Epub 2013 Feb 4.
Mediation analysis allowing for exposure-mediator interactions and causal interpretation: theoretical assumptions and implementation with SAS and SPSS macros. */


*mediation analysis*;
%macro mediation(data=,yvar=,avar=,mvar=,cvar=,a0=,a1=,m=,nc=, yreg=,mreg=,
interaction=,casecontrol=false,output=reduced,c=,boot=,cens=);

*/data house keeping*/;

/* 2020-03-30 Modified by @kaz-yos */
/* https://communities.sas.com/t5/SAS-Programming/Applying-countw-to-an-empty-string/td-p/15285 */
%if (&cvar^=) %then %do;
/* This only works with an non-empty &cvar. */
%let nc=%sysfunc(countw(&cvar));
%end;
%else %do;
/* Empty case */
%let nc=0;
%end;

 %put There are &nc confounders in the string "&cvar";



data data1;
set &data (keep=&yvar &mvar &avar &cvar &cens);
run;
	%if &interaction=true %then %do;
data data1;
set data1;
int=&avar*&mvar;
run;
	%end;

	%if (&cvar^= & &casecontrol=false) | (&cvar^= & &casecontrol=) %then %do;
		%LET cvars= &cvar;
		%LET i =1;
		%DO %UNTIL(NOT %LENGTH(%SCAN(&cvars,&i))) ;
proc means noprint data=data1;
var %SCAN(&cvars,&i);
output out=data2&i mean=/autoname;
run;
data data2&i;
set data2&i;
drop _TYPE_ _FREQ_;
run;
proc iml;
use data2&i;
read all into vb;
mean=vb[1,1];
cname1 = {"mean"};
create data2new&i from mean [colname=cname1];
append from mean;
quit;
proc append base=data3 data=data2new&i;
run;
proc sql;
		%LET i=%EVAL(&i+1);
		%END;
proc iml;
use data3;
read all into vb;
data3=t(vb);
create data2 from data3;
append from data3;
quit;
		%if &c^= %then %do;
			%LET cval= &c;
			%LET i =1 ;
			%DO %UNTIL(NOT %LENGTH(%SCAN(&cval,&i))) ;
proc sql;
create table data2 as
select *, mean(%SCAN(&cval,&i)) as cval&i
from data2
run;
			%LET i=%EVAL(&i+1);
			%END;
		%end;
	%end;




		%if (&cvar^= & &casecontrol=true) %then %do;
		%LET cvars= &cvar;
		%LET i =1;
		%DO %UNTIL(NOT %LENGTH(%SCAN(&cvars,&i))) ;
proc means noprint data=data1;
where &yvar=0;
var %SCAN(&cvars,&i);
output out=data2&i mean=/autoname;
run;
data data2&i;
set data2&i;
drop _TYPE_ _FREQ_;
run;
proc iml;
use data2&i;
read all into vb;
mean=vb[1,1];
cname1 = {"mean"};
create data2new&i from mean [colname=cname1];
append from mean;
quit;
proc append base=data3 data=data2new&i;
run;
proc sql;
		%LET i=%EVAL(&i+1);
		%END;
proc iml;
use data3;
read all into vb;
data3=t(vb);
create data2 from data3;
append from data3;
quit;
		%if &c^= %then %do;
			%LET cval= &c;
			%LET i =1 ;
			%DO %UNTIL(NOT %LENGTH(%SCAN(&cval,&i))) ;
proc sql;
create table data2 as
select *, mean(%SCAN(&cval,&i)) as cval&i
from data2
run;
			%LET i=%EVAL(&i+1);
			%END;
		%end;
	%end;






***************************   BOOTSTRAP PROCEDURE   ******************************************************************;

	%if (&boot^= & &boot^=false) %then %do;
	*DMSLOGSIZE=MAX;
	%if &boot=true %then %do;
%LET n = 1000;
	%end;
%if &boot^=true %then %do;
%LET n = &boot;
	%end;
******************* bootstrap samples******************************;
data data1;

do sample = 1 to &n; /* To create b bootstrap replications */
do i = 1 to nobs;
indexbootstrap = round(ranuni(0) * nobs);
set data1
nobs = nobs
point = indexbootstrap;
output;
end;
end;
stop;
run;
		%if &interaction=true %then %do;
data data1;
set data1;
int=&avar*&mvar;
run;
		%end;


  			%do t=1 %to &n;

   data data1&t;
   set data1(where=(sample=&t));
   run;

		%end;

	%end;



***************** regression-for bootstrap *************************;
	%if (&boot^= & &boot^=false) %then %do;

		%do t=1 %to &n;

************************************************************************************************************************;
			%if &yreg=linear  %then %do;
************************************************************************************************************************;
				%if &interaction=false & &cvar^= %then %do;
proc reg data=data1&t covout noprint
outest=out1&t(drop=_model_ _type_ _name_ _depvar_ _rmse_ &yvar) ;
model  &yvar=&avar &mvar &cvar ;
proc print;
run;

				%end;
				%if &interaction=false & &cvar= %then %do;
proc reg data=data1&t covout noprint
outest=out1&t(drop=_model_ _type_ _name_ _depvar_ _rmse_ &yvar) ;
model  &yvar=&avar &mvar;
run;
				%end;
				%if &interaction=true & &cvar^= %then %do;
proc reg data=data1&t covout noprint
outest=out1&t(drop=_model_ _type_ _name_ _depvar_ _rmse_ &yvar) ;
model  &yvar=&avar &mvar int &cvar ;
run;

				%end;

				%if &interaction=true & &cvar= %then %do;
proc reg data=data1&t covout noprint
outest=out1&t(drop=_model_ _type_ _name_ _depvar_ _rmse_ &yvar)  ;
model  &yvar=&avar &mvar int ;

run;
				%end; * if interaction;
			%end; *ylinear;

***********************************************************************************************************************************************************************************************************************************;
			%if &yreg=logistic  | &yreg=loglinear  |&yreg=poisson | &yreg=negbin |&yreg=survCox |&yreg=survAdd |&yreg=survAFT_weibull |&yreg=survAFT_exp |&yreg=survAFT_gamma |&yreg=survAFT_loglogistic  |&yreg=survAFT_normal %then %do;
***********************************************************************************************************************************************************************************************************************************;

*need to include models for survival outcome!!!;

				%if &interaction=true & &cvar^= %then %do;
					%if &yreg=logistic %then %do;
proc logistic  data=data1&t descending covout noprint
outest=out1&t(drop=_link_ _status_ _LNLIKE_ _type_ _name_ )  ;
model  &yvar=&avar &mvar int &cvar ;
run;
					%end;
					%if &yreg=loglinear %then %do;
proc genmod data=data1&t descending ;
model &yvar=&avar &mvar  int &cvar/dist=binomial link=log covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[,];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-1,2];
par=t(par);
out=par//cov;
create out1&t from out;
append from out;
quit;
					%end;
					%if &yreg=poisson %then %do;
proc genmod data=data1&t  ;
model &yvar=&avar &mvar  int &cvar/dist=poisson covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[,];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-1,2];
par=t(par);
out=par//cov;
create out1&t from out;
append from out;
quit;
					%end;
					%if &yreg=negbin %then %do;
proc genmod data=data1&t  ;
model &yvar=&avar &mvar  int &cvar/dist=negbin covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[1:nrow(vb)-1,1:ncol(vb)-1];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-1,2];
par=t(par);
out=par//cov;
create out1&t from out;
append from out;
quit;
					%end;
*/add survival*/;
*/note COX does not have intercept!!!! I add a column of zeros and a row of zeros*/;
	%if &yreg=survCox %then %do;
	proc phreg data = data1&t ;
model &yvar*&cens(1) = &avar &mvar int &cvar/covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
			run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[,];
use gmparms;
read all into vb;
par=vb[1:nrow(vb),2];
par=t(par);
x=vb[1:nrow(vb),2]*0//0//0;
rzero=t(vb[1:nrow(vb),2]*0);
out=par//rzero//cov;
out=x||out;
create out1&t from out;
append from out;
quit;

%end;


%if &yreg=survAFT_exp %then %do;
proc lifereg data = data1&t ;
   model &yvar*&cens(1) = &avar &mvar int &cvar /covb DISTRIBUTION=exponential ;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
			run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[1:nrow(vb)-1,1:nrow(vb)-1];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-2,2];
par=t(par);
out=par//cov;
create out1&t from out;
append from out;
quit;
%end;

%if &yreg=survAFT_weibull %then %do;
proc lifereg data = data1&t ;
   model &yvar*&cens(1) = &avar &mvar int &cvar /covb DISTRIBUTION=weibull;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
			run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[1:nrow(vb)-1,1:nrow(vb)-1];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-2,2];
par=t(par);
out=par//cov;
create out1&t from out;
append from out;
quit;
%end;

				%end;

				%if &interaction=true & &cvar= %then %do;
					%if &yreg=logistic %then %do;
proc logistic  data=data1&t descending covout noprint
outest=out1&t(drop=_link_ _status_ _LNLIKE_ _type_ _name_ )  ;
model  &yvar=&avar &mvar int ;
run;
					%end;
					%if &yreg=loglinear %then %do;
proc genmod data=data1&t descending ;
model &yvar=&avar &mvar int /dist=binomial link=log covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[,];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-1,2];
par=t(par);
out=par//cov;
create out1&t from out;
append from out;
quit;
					%end;
					%if &yreg=poisson %then %do;
proc genmod data=data1&t ;
model &yvar=&avar &mvar  int /dist=poisson covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[,];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-1,2];
par=t(par);
out=par//cov;
create out1&t from out;
append from out;
quit;
					%end;
					%if &yreg=negbin %then %do;
proc genmod data=data1&t  ;
model &yvar=&avar &mvar  int /dist=negbin covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[1:nrow(vb)-1,1:ncol(vb)-1];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-1,2];
par=t(par);
out=par//cov;
create out1&t from out;
append from out;
quit;
					%end;
*/add survival*/;
*/note COX does not have intercept!!!! I add a column of zeros*/;
	%if &yreg=survCox %then %do;
	proc phreg data = data1&t ;
model &yvar*&cens(1) = &avar &mvar int /covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
			run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[,];
use gmparms;
read all into vb;
par=vb[1:nrow(vb),2];
par=t(par);
x=vb[1:nrow(vb),2]*0//0//0;
rzero=t(vb[1:nrow(vb),2]*0);
out=par//rzero//cov;
out=x||out;
create out1&t from out;
append from out;
quit;

%end;

%if &yreg=survAFT_exp %then %do;
proc lifereg data = data1&t ;
   model &yvar*&cens(1) = &avar &mvar int /covb DISTRIBUTION=exponential ;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
			run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[1:nrow(vb)-1,1:nrow(vb)-1];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-2,2];
par=t(par);
out=par//cov;
create out1&t from out;
append from out;
quit;
%end;

%if &yreg=survAFT_weibull %then %do;
proc lifereg data = data1&t ;
   model &yvar*&cens(1) = &avar &mvar int /covb DISTRIBUTION=weibull;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
			run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[1:nrow(vb)-1,1:nrow(vb)-1];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-2,2];
par=t(par);
out=par//cov;
create out1&t from out;
append from out;
quit;
%end;
				%end;

				%if &interaction=false & &cvar^= %then %do;
					%if &yreg=logistic %then %do;
proc logistic  data=data1&t descending covout noprint
outest=out1&t(drop=_link_ _status_ _LNLIKE_ _type_ _name_ )  ;
model  &yvar=&avar &mvar  &cvar ;
run;
					%end;
					%if &yreg=loglinear %then %do;
proc genmod data=data1&t descending ;
model &yvar=&avar &mvar  &cvar/dist=binomial link=log covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[,];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-1,2];
par=t(par);
out=par//cov;
create out1&t from out;
append from out;
quit;
					%end;
					%if &yreg=poisson %then %do;
proc genmod data=data1&t  ;
model &yvar=&avar &mvar  &cvar/dist=poisson covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[,];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-1,2];
par=t(par);
out=par//cov;
create out1&t from out;
append from out;
quit;
					%end;
					%if &yreg=negbin %then %do;
proc genmod data=data1&t ;
model &yvar=&avar &mvar  &cvar/dist=negbin covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[1:nrow(vb)-1,1:ncol(vb)-1];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-1,2];
par=t(par);
out=par//cov;
create out1&t from out;
append from out;
quit;
					%end;
*/add survival*/;
*/note COX does not have intercept!!!! I add a column of zeros*/;
	%if &yreg=survCox %then %do;
	proc phreg data = data1&t ;
model &yvar*&cens(1) = &avar &mvar &cvar/covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
			run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[,];
use gmparms;
read all into vb;
par=vb[1:nrow(vb),2];
par=t(par);
x=vb[1:nrow(vb),2]*0//0//0;
rzero=t(vb[1:nrow(vb),2]*0);
out=par//rzero//cov;
out=x||out;
create out1&t from out;
append from out;
quit;

%end;

%if &yreg=survAFT_exp %then %do;
proc lifereg data = data1&t ;
   model &yvar*&cens(1) = &avar &mvar &cvar/covb DISTRIBUTION=exponential ;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
			run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[1:nrow(vb)-1,1:nrow(vb)-1];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-2,2];
par=t(par);
out=par//cov;
create out1&t from out;
append from out;
quit;
%end;

%if &yreg=survAFT_weibull %then %do;
proc lifereg data = data1&t ;
   model &yvar*&cens(1) = &avar &mvar &cvar/covb DISTRIBUTION=weibull ;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
			run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[1:nrow(vb)-1,1:nrow(vb)-1];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-2,2];
par=t(par);
out=par//cov;
create out1&t from out;
append from out;
quit;
%end;
				%end;

				%if &interaction=false & &cvar= %then %do;
					%if &yreg=logistic %then %do;
proc logistic  data=data1&t descending covout noprint
outest=out1&t(drop=_link_ _status_ _LNLIKE_ _type_ _name_ )  ;
model  &yvar=&avar &mvar  ;
run;
					%end;
					%if &yreg=loglinear %then %do;
proc genmod data=data1 descending ;
model &yvar=&avar &mvar /dist=binomial link=log covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[,];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-1,2];
par=t(par);
out=par//cov;
create out1&t from out;
append from out;
quit;
					%end;
					%if &yreg=poisson %then %do;
proc genmod data=data1&t  ;
model &yvar=&avar &mvar /dist=poisson covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[,];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-1,2];
par=t(par);
out=par//cov;
create out1&t from out;
append from out;
quit;
					%end;
					%if &yreg=negbin %then %do;
proc genmod data=data1&t ;
model &yvar=&avar &mvar /dist=negbin covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[1:nrow(vb)-1,1:ncol(vb)-1];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-1,2];
par=t(par);
out=par//cov;
create out1&t from out;
append from out;
quit;
					%end;
*/add survival*/;
*/note COX does not have intercept!!!! I add a column of zeros*/;
	%if &yreg=survCox %then %do;
	proc phreg data = data1&t ;
model &yvar*&cens(1) = &avar &mvar /covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
			run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[,];
use gmparms;
read all into vb;
par=vb[1:nrow(vb),2];
par=t(par);
x=vb[1:nrow(vb),2]*0//0//0;
rzero=t(vb[1:nrow(vb),2]*0);
out=par//rzero//cov;
out=x||out;
create out1&t from out;
append from out;
quit;

%end;

%if &yreg=survAFT_exp %then %do;
proc lifereg data = data1&t ;
   model &yvar*&cens(1) = &avar &mvar /covb DISTRIBUTION=exponential ;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
			run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[1:nrow(vb)-1,1:nrow(vb)-1];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-2,2];
par=t(par);
out=par//cov;
create out1&t from out;
append from out;
quit;
%end;

%if &yreg=survAFT_weibull %then %do;
proc lifereg data = data1&t ;
   model &yvar*&cens(1) = &avar &mvar /covb DISTRIBUTION=weibull ;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
			run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[1:nrow(vb)-1,1:nrow(vb)-1];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-2,2];
par=t(par);
out=par//cov;
create out1&t from out;
append from out;
quit;
%end;
				%end;
			%end;


************************************************************************************************************************;
			%if  &mreg=linear & &yreg=linear %then %do;
************************************************************************************************************************;

				%if &cvar^= %then %do;
					%if &casecontrol^=true %then %do;
proc reg data=data1&t covout noprint
outest=out2&t(drop=_model_ _type_ _name_ _depvar_ _rmse_ &mvar)  ;
model  &mvar=&avar &cvar ;
proc print;
run;
					%end;
					%if &casecontrol=true %then %do;
proc reg data=data1&t covout noprint
outest=out2&t(drop=_model_ _type_ _name_ _depvar_ _rmse_ &mvar)  ;
where &yvar=0;
model  &mvar=&avar &cvar;
proc print;
run;
					%end;
				%end;
				%if  &cvar= %then %do;
					%if &casecontrol^=true %then %do;
proc reg data=data1&t covout noprint
outest=out2&t(drop=_model_ _type_ _name_ _depvar_ _rmse_ &mvar)   ;
model  &mvar=&avar;

run;
					%end;
					%if &casecontrol=true %then %do;
proc reg data=data1&t covout noprint
outest=out2&t(drop=_model_ _type_ _name_ _depvar_ _rmse_ &mvar)  ;
where &yvar=0;
model  &mvar=&avar;
proc print;
run;
					%end;
				%end;
			%end;
************************************************************************************************************************;
			%if  &mreg=linear & &yreg^=linear %then %do;
************************************************************************************************************************;

				%if &cvar^= %then %do;
					%if &casecontrol^=true %then %do;
proc reg data=data1&t covout noprint
outest=out2&t(drop=_model_ _type_ _name_ _depvar_  &mvar)  ;
model  &mvar=&avar &cvar ;
proc print;
run;
					%end;
					%if &casecontrol=true %then %do;
proc reg data=data1&t covout noprint
outest=out2&t(drop=_model_ _type_ _name_ _depvar_  &mvar)  ;
where &yvar=0;
model  &mvar=&avar &cvar;
proc print;
run;
					%end;
				%end;
				%if  &cvar= %then %do;
					%if &casecontrol^=true %then %do;
proc reg data=data1&t covout noprint
outest=out2&t(drop=_model_ _type_ _name_ _depvar_  &mvar)  ;
model  &mvar=&avar;
proc print;
run;
					%end;
					%if &casecontrol=true %then %do;
proc reg data=data1&t covout noprint
outest=out2&t(drop=_model_ _type_ _name_ _depvar_ &mvar)  ;
where &yvar=0;
model  &mvar=&avar;
proc print;
run;
					%end;
				%end;
			%end;



************************************************************************************************************************;
			%if &mreg=logistic %then %do;
************************************************************************************************************************;

				%if &cvar^= %then %do;
					%if &casecontrol^=true %then %do;
proc logistic data=data1&t descending covout noprint
outest=out2&t(drop=_link_ _status_ _LNLIKE_ _type_ _name_ )  ;
model  &mvar=&avar &cvar ;

run;
					%end;
					%if &casecontrol=true %then %do;
proc logistic data=data1&t descending covout noprint
outest=out2&t(drop=_link_ _status_ _LNLIKE_ _type_ _name_ )  ;
where &yvar=0;
model  &mvar=&avar &cvar ;
run;
					%end;
				%end;
				%if &cvar= %then %do;
					%if &casecontrol^=true %then %do;
proc logistic data=data1&t descending covout noprint
outest=out2&t(drop=_link_ _status_ _LNLIKE_ _type_ _name_ )  ;
model  &mvar=&avar  ;
run;
					%end;
					%if &casecontrol=true %then %do;
proc logistic data=data1&t descending covout noprint
outest=out2&t(drop=_link_ _status_ _LNLIKE_ _type_ _name_ )  ;
where &yvar=0;
model  &mvar=&avar  ;
run;
					%end;
				%end;
			%end;
		%end;

***************** regression-for bootstrap 	END *************************;




***************** causal effects for bootstrap  *************************;

 */create objects in which we save the bootstrap samples of causal effects*/;
proc iml;


		%if &mreg=linear & &interaction=false  %then %do;
bootsample=J(&n,3,0);
		%end;

		%if &mreg=linear & &interaction=true %then %do;
			%if &cvar^= & &output=full %then %do;
bootsample=J(&n,12,0);
			%end;
			%if  &cvar= | (&cvar^= & &output^=full) %then %do;
bootsample=J(&n,6,0);
			%end;
		%end;

		%if &mreg=logistic %then %do;

			%if &cvar^= & &output=full %then %do;
bootsample=J(&n,12,0);
			%end;
			%if &interaction=false & &cvar= & &yreg=linear %then %do;
bootsample=J(&n,3,0);
			%end;
			%if (&interaction=true & &cvar=) | (&cvar^= & &output^=full) | (&interaction=false & &cvar= & &yreg^=linear ) %then %do;
bootsample=J(&n,6,0);
			%end;

		%end;





*/compute the causal effects*/;
		%if (&mreg=linear & &interaction=false ) | (&yreg=linear & &mreg=logistic & &interaction=false & &cvar=)  %then %do;


			%do t=1 %to &n;

USE out2&t;
READ ALL INTO VB;
				%if (&yreg=linear) %then %do;
beta0= VB[1,1];
beta1=VB[1,2];
				%end;
				%if (&yreg^=linear) %then %do;
beta1=VB[1,3];
				%end;
USE out1&t;
READ ALL INTO VB;
theta1=VB[1,2];
theta2=VB[1,3];
*/cde and nde*/;
				%if (&yreg=linear & &mreg=logistic) %then %do;
bootsample[&t,1]=(theta1)*(&a1-&a0);
*/nie*/;
bootsample[&t,2]=(theta2)*(exp(beta0+beta1*&a1)/(1+exp(beta0+beta1*&a1))-exp(beta0+beta1*&a0)/(1+exp(beta0+beta1*&a0)));
*/te*/;
bootsample[&t,3]=bootsample[&t,1]+bootsample[&t,2];
				%end;
				%if (&yreg=linear & &mreg=linear & &interaction=false) %then %do;
bootsample[&t,1]=((theta1)*(&a1-&a0));
*/nie*/;
bootsample[&t,2]=((theta2*beta1)*(&a1-&a0));
*/te*/;
bootsample[&t,3]=((theta1+theta2*beta1)*(&a1-&a0));
				%end;
				%if (&yreg^=linear & &mreg=linear ) %then %do;
bootsample[&t,1]=exp((theta1)*(&a1-&a0));
bootsample[&t,2]=exp((theta2*beta1)*(&a1-&a0));
bootsample[&t,3]=bootsample[&t,1]*bootsample[&t,2];
				%end;
			%end;
x=bootsample;
cname1 = { "boot1" "boot2" "boot3"};
create bootdata from x [colname=cname1];
append from x;
		%end;* noint;


		%if (&interaction=true) | (&mreg=logistic & &interaction=false & &cvar^=) |(&yreg^=linear & &mreg=logistic & &interaction=false) %then %do;

			%if &cvar^= %then %do;
USE data2;
read all into vb;
				%if &c= %then %do;
cmean=VB[1,1:ncol(vb)];
				%end;

				%if &c^= %then %do;
cmean=VB[1,1:ncol(vb)-&nc];
c=VB[1,ncol(vb)-&nc+1:ncol(vb)] ;
				%end;
			%end;

			%if (&cvar^=) %then %do;
				%do t=1 %to &n;

USE out1&t;
READ ALL INTO VB;
theta1=VB[1,2];
theta2=VB[1,3];
					%if &interaction=true %then %do;
theta3=VB[1,4] ;
					%end;
USE out2&t;
READ ALL INTO VB;
					%if (&yreg=linear & &mreg=linear) | (&mreg=logistic) %then %do;
beta0=VB[1,1];
beta1=VB[1,2];
beta2= VB[1,3:ncol(vb)];
					%end;
					%if (&yreg^=linear & &mreg=linear) %then %do;
s2=VB[1,1];
s2=s2**2;
beta0=VB[1,2];
beta1=VB[1,3];
beta2= VB[1,4:ncol(vb)];
tsq=(theta3**2);
rm=s2;
asq=(&a1**2);
a1sq=(&a0**2);
					%end;


					%if (&yreg=linear & &mreg=linear & &interaction=true) %then %do;
print cmean;
*/MARGINAL CDE*/;
bootsample[&t,1]=(theta1+theta3*&m)*(&a1-&a0);
*/MARGINAL NDE*/;
bootsample[&t,2]=(theta1+theta3*beta0+theta3*beta1*&a0+(theta3*beta2*t(cmean)))*(&a1-&a0);
*/MARGINAL NIE*/;
bootsample[&t,3]=(theta2*beta1+theta3*beta1*&a0)*(&a1-&a0);
*/ MARGINAL TNDE*/;
bootsample[&t,4]=(theta1+theta3*beta0+theta3*beta1*&a1+(theta3*beta2*t(cmean)))*(&a1-&a0);
*/ MARGINAL TNIE*/;
bootsample[&t,5]=(theta2*beta1+theta3*beta1*&a1)*(&a1-&a0);
*/te marginal*/;
bootsample[&t,6]=(theta1+theta3*beta0+theta3*beta1*&a0+(theta3*beta2*t(cmean))+theta2*beta1+theta3*beta1*&a1)*(&a1-&a0);
						%if &c^= %then %do;
*/CONDITIONAL CDE*/;
bootsample[&t,7]=(theta1)*(&a1-&a0)+(theta3*(&m))*(&a1-&a0);
*/CONDITIONAL NDE*/;
bootsample[&t,8]=(theta1+theta3*beta0+theta3*beta1*&a0+(theta3*beta2*t(c)))*(&a1-&a0);
*/CONDITIONAL NIE*/;
bootsample[&t,9]=(theta2*beta1+theta3*beta1*&a0)*(&a1-&a0);
*/CONDITIONAL TNDE*/;
bootsample[&t,10]=(theta1+theta3*beta0+theta3*beta1*&a1+(theta3*beta2*t(c)))*(&a1-&a0);
*/ CONDITIONAL TNIE*/;
bootsample[&t,11]=(theta2*beta1+theta3*beta1*&a1)*(&a1-&a0);
*/te conditional*/;
bootsample[&t,12]=(theta1+theta3*beta0+theta3*beta1*&a0+(theta3*beta2*t(c))+theta2*beta1+theta3*beta1*&a1)*(&a1-&a0);
						%end;
					%end;
					%if (&yreg=linear & &mreg=logistic & &interaction=true)  %then %do;
*/MARGINAL CDE*/;
bootsample[&t,1]=(theta1+theta3*&m)*(&a1-&a0);
*/MARGINAL NDE*/;
bootsample[&t,2]=(theta1+theta3*exp(beta0+beta1*&a0+sum(beta2*t(cmean)))/(1+exp(beta0+beta1*&a0+sum(beta2*t(cmean)))))*(&a1-&a0);
*/MARGINAL NIE*/;
bootsample[&t,3]=(theta2+theta3*&a0)*(
exp(beta0+beta1*&a1+sum(beta2*t(cmean)))/
(1+exp(beta0+beta1*&a1+sum(beta2*t(cmean))))-
exp(beta0+beta1*&a0+sum(beta2*t(cmean)))/
(1+exp(beta0+beta1*&a0+sum(beta2*t(cmean))))
);
*/ MARGINAL TNDE*/;
bootsample[&t,4]=(theta1+theta3*exp(beta0+beta1*&a1+sum(beta2*t(cmean)))/(1+exp(beta0+beta1*&a1+sum(beta2*t(cmean)))))*(&a1-&a0);
*/ MARGINAL TNIE*/;
bootsample[&t,5]=(theta2+theta3*&a1)*(exp(beta0+beta1*&a1+sum(beta2*t(cmean)))/(1+exp(beta0+beta1*&a1+sum(beta2*t(cmean))))-exp(beta0+beta1*&a0+sum(beta2*t(cmean)))/(1+exp(beta0+beta1*&a0+sum(beta2*t(cmean)))));
*/te marginal*/;
bootsample[&t,6]=bootsample[&t,2]+bootsample[&t,5];
						%if &c^= %then %do;
*/CONDITIONAL CDE*/;
bootsample[&t,7]=(theta1)*(&a1-&a0)+(theta3*(&m))*(&a1-&a0);
*/CONDITIONAL NDE*/;
bootsample[&t,8]=(theta1+theta3*exp(beta0+beta1*&a0+sum(beta2*t(c)))/(1+exp(beta0+beta1*&a0+sum(beta2*t(c)))))*(&a1-&a0);
*/CONDITIONAL NIE*/;
bootsample[&t,9]=(theta2+theta3*&a0)*(
exp(beta0+beta1*&a1+sum(beta2*t(c)))/
(1+exp(beta0+beta1*&a1+sum(beta2*t(c))))-
exp(beta0+beta1*&a0+sum(beta2*t(c)))/
(1+exp(beta0+beta1*&a0+sum(beta2*t(c))))
);
*/CONDITIONAL TNDE*/;
bootsample[&t,10]=(theta1+theta3*exp(beta0+beta1*&a1+sum(beta2*t(c)))/(1+exp(beta0+beta1*&a1+sum(beta2*t(c)))))*(&a1-&a0);
*/ CONDITIONAL TNIE*/;
bootsample[&t,11]=(theta2+theta3*&a1)*(
exp(beta0+beta1*&a1+sum(beta2*t(c)))/
(1+exp(beta0+beta1*&a1+sum(beta2*t(c))))-
exp(beta0+beta1*&a0+sum(beta2*t(c)))/
(1+exp(beta0+beta1*&a0+sum(beta2*t(c))))
);
*/te conditional*/;
bootsample[&t,12]=bootsample[&t,8]+bootsample[&t,11];
						%end;
					%end;
					%if (&yreg=linear & &mreg=logistic & &interaction=false) %then %do;
*/MARGINAL CDE*/;
bootsample[&t,1]=(theta1)*(&a1-&a0);
*/MARGINAL NDE*/;
bootsample[&t,2]=(theta1)*(&a1-&a0);
*/MARGINAL NIE*/;
bootsample[&t,3]=(theta2)*(exp(beta0+beta1*&a1+sum(beta2*t(cmean)))/(1+exp(beta0+beta1*&a1+sum(beta2*t(cmean))))-exp(beta0+beta1*&a0+sum(beta2*t(cmean)))/(1+exp(beta0+beta1*&a0+sum(beta2*t(cmean)))));
*/ MARGINAL TNDE*/;
bootsample[&t,4]=(theta1)*(&a1-&a0);
*/ MARGINAL TNIE*/;
bootsample[&t,5]=(theta2)*(exp(beta0+beta1*&a1+sum(beta2*t(cmean)))/(1+exp(beta0+beta1*&a1+sum(beta2*t(cmean))))-exp(beta0+beta1*&a0+sum(beta2*t(cmean)))/(1+exp(beta0+beta1*&a0+sum(beta2*t(cmean)))));
*/te marginal*/;
bootsample[&t,6]=bootsample[&t,2]+bootsample[&t,5];
						%if &c^= %then %do;
*/CONDITIONAL CDE*/;
bootsample[&t,7]=(theta1)*(&a1-&a0);
*/CONDITIONAL NDE*/;
bootsample[&t,8]=(theta1)*(&a1-&a0);
*/CONDITIONAL NIE*/;
bootsample[&t,9]=(theta2)*(exp(beta0+beta1*&a1+sum(beta2*t(c)))/(1+exp(beta0+beta1*&a1+sum(beta2*t(c))))-exp(beta0+beta1*&a0+sum(beta2*t(c)))/(1+exp(beta0+beta1*&a0+sum(beta2*t(c)))));
*/CONDITIONAL TNDE*/;
bootsample[&t,10]=(theta1)*(&a1-&a0);
*/ CONDITIONAL TNIE*/;
bootsample[&t,11]=(theta2)*(exp(beta0+beta1*&a1+sum(beta2*t(c)))/(1+exp(beta0+beta1*&a1+sum(beta2*t(c))))-exp(beta0+beta1*&a0+sum(beta2*t(c)))/(1+exp(beta0+beta1*&a0+sum(beta2*t(c)))));
*/te conditional*/;
bootsample[&t,12]=bootsample[&t,8]+bootsample[&t,11];
						%end;
					%end;
					%if (&yreg^=linear & &mreg=linear & &interaction=true) %then %do;
  	  */MARGINAL CDE*/;
	  x6=(theta1+theta3*&m)*(&a1-&a0);
bootsample[&t,1]=exp(x6);
	   */MARGINAL NDE*/;
	  x7=(theta1+theta3*beta0+theta3*beta1*&a0+sum(theta3*beta2*t(cmean))+theta3*theta2*rm)*(&a1-&a0)+1/2*tsq*rm*(asq-a1sq);
bootsample[&t,2]=exp(x7);
	  */MARGINAL NIE*/;
	  x8=(theta2*beta1+theta3*beta1*&a0)*(&a1-&a0);
bootsample[&t,3]=exp(x8);
	  */ MARGINAL TNDE*/;
	  x9=(theta1+theta3*beta0+theta3*beta1*&a1+sum(theta3*beta2*t(cmean))+theta3*theta2*rm)*(&a1-&a0)+1/2*tsq*rm*(asq-a1sq);
bootsample[&t,4]=exp(x9);
	  */ MARGINAL TNIE*/;
	  x10=(theta2*beta1+theta3*beta1*&a1)*(&a1-&a0);
bootsample[&t,5]=exp(x10);
	  */te*/;
bootsample[&t,6]=bootsample[&t,2]*bootsample[&t,5];
						%if &c^= %then %do;
*/CONDITIONAL CDE*/;
bootsample[&t,7]=exp((theta1+theta3*&m)*(&a1-&a0));
      */CONDITIONAL NDE*/;
bootsample[&t,8]=exp((theta1+theta3*beta0+theta3*beta1*&a0+sum(theta3*beta2*t(c))+theta3*theta2*rm)*(&a1-&a0)+(1/2)*tsq*rm*(asq-a1sq)
);
      */CONDITIONAL NIE*/;
	  x3=(theta2*beta1+theta3*beta1*&a0)*(&a1-&a0);
bootsample[&t,9]=exp(x3);
      */CONDITIONAL TNDE*/;
	  x4=(theta1+theta3*beta0+theta3*beta1*&a1+sum(theta3*beta2*t(c))+theta3*theta2*rm)*(&a1-&a0)+(1/2)*tsq*rm*(asq-a1sq);
bootsample[&t,10]=exp(x4);
      */ CONDITIONAL TNIE*/;
	  x5=(theta2*beta1+theta3*beta1*&a1)*(&a1-&a0);
bootsample[&t,11]=exp(x5);
	  */te*/;
bootsample[&t,12]=bootsample[&t,8]*bootsample[&t,11];

						%end;
					%end;
					%if (&yreg^=linear & &mreg=logistic & &interaction=false) %then %do;
*/MARGINAL CDE*/;
	  x6=(theta1)*(&a1-&a0);
bootsample[&t,1]=exp(x6);
	   */MARGINAL NDE*/;
bootsample[&t,2]=exp((theta1)*(&a1-&a0))*
(1+exp(theta2+beta0+beta1*&a0+sum(beta2*t(cmean))))/
(1+exp(theta2+beta0+beta1*&a0+sum(beta2*t(cmean))));
	  */MARGINAL NIE*/;
bootsample[&t,3]=(
(1+exp(beta0+beta1*&a0+sum(beta2*t(cmean))))*
(1+exp(theta2+beta0+beta1*&a1+sum(beta2*t(cmean))))
)/
(
(1+exp(beta0+beta1*&a1+sum(beta2*t(cmean))))*
(1+exp(theta2+beta0+beta1*&a0+sum(beta2*t(cmean))))
);
	  */ MARGINAL TNDE*/;
bootsample[&t,4]=exp((theta1)*(&a1-&a0))*
(1+exp(theta2+beta0+beta1*&a1+sum(beta2*t(cmean))))/
(1+exp(theta2+beta0+beta1*&a1+sum(beta2*t(cmean))));
	  */ MARGINAL TNIE*/;
bootsample[&t,5]=(
(1+exp(beta0+beta1*&a0+sum(beta2*t(cmean))))*
(1+exp(theta2+beta0+beta1*&a1+sum(beta2*t(cmean))))
)/
(
(1+exp(beta0+beta1*&a1+sum(beta2*t(cmean))))*
(1+exp(theta2+beta0+beta1*&a0+sum(beta2*t(cmean))))
);
bootsample[&t,6]=bootsample[&t,2]*bootsample[&t,5];
						%if &c^= %then %do;
*/CONDITIONAL CDE*/;
	  x1=exp(theta1*(&a1-&a0));
bootsample[&t,7]=x1;
      */CONDITIONAL NDE*/;
bootsample[&t,8]=exp((theta1)*(&a1-&a0))*
(1+exp(theta2+beta0+beta1*&a0+sum(beta2*t(c))))/
(1+exp(theta2+beta0+beta1*&a0+sum(beta2*t(c))));
      */CONDITIONAL NIE*/;
bootsample[&t,9]=(
(1+exp(beta0+beta1*&a0+beta2*t(c)))*
(1+exp(theta2+beta0+beta1*&a1+sum(beta2*t(c))))
)/
(
(1+exp(beta0+beta1*&a1+sum(beta2*t(c))))*
(1+exp(theta2+beta0+beta1*&a0+sum(beta2*t(c))))
);
      */CONDITIONAL TNDE*/;
bootsample[&t,10]=exp((theta1)*(&a1-&a0))*
(1+exp(theta2+beta0+beta1*&a1+sum(beta2*t(c))))/
(1+exp(theta2+beta0+beta1*&a1+sum(beta2*t(c))));
      */ CONDITIONAL TNIE*/;
bootsample[&t,11]=(
(1+exp(beta0+beta1*&a0+sum(beta2*t(c))))*
(1+exp(theta2+beta0+beta1*&a1+sum(beta2*t(c))))
)/
(
(1+exp(beta0+beta1*&a1+sum(beta2*t(c))))*
(1+exp(theta2+beta0+beta1*&a0+sum(beta2*t(c))))
);
bootsample[&t,12]=bootsample[&t,8]*bootsample[&t,11];
						%end;
					%end;
					%if (&yreg^=linear & &mreg=logistic & &interaction=true) %then %do;
*/MARGINAL CDE*/;
x6=(theta1+theta3*&m)*(&a1-&a0);
bootsample[&t,1]=exp(x6);
*/MARGINAL NDE*/;
bootsample[&t,2]=exp(theta1*(&a1-&a0))*
(1+exp(theta2+theta3*&a1+beta0+beta1*&a0+sum(beta2*t(cmean))))/
(1+exp(theta2+theta3*&a0+beta0+beta1*&a0+sum(beta2*t(cmean))));
*/MARGINAL NIE*/;
bootsample[&t,3]=(
(1+exp(beta0+beta1*&a0+sum(beta2*t(cmean))))*
(1+exp(theta2+theta3*&a0+beta0+beta1*&a1+sum(beta2*t(cmean))))
)/
(
(1+exp(beta0+beta1*&a1+sum(beta2*t(cmean))))*
(1+exp(theta2+theta3*&a0+beta0+beta1*&a0+sum(beta2*t(cmean))))
);
*/ MARGINAL TNDE*/;
bootsample[&t,4]=exp(theta1*(&a1-&a0))*
(1+exp(theta2+theta3*&a1+beta0+beta1*&a1+sum(beta2*t(cmean))))/
(1+exp(theta2+theta3*&a0+beta0+beta1*&a1+sum(beta2*t(cmean))));
*/ MARGINAL TNIE*/;
bootsample[&t,5]=(
(1+exp(beta0+beta1*&a0+sum(beta2*t(cmean))))*
(1+exp(theta2+theta3*&a1+beta0+beta1*&a1+sum(beta2*t(cmean))))
)/
(
(1+exp(beta0+beta1*&a1+sum(beta2*t(cmean))))*
(1+exp(theta2+theta3*&a1+beta0+beta1*&a0+sum(beta2*t(cmean))))
);
bootsample[&t,6]=bootsample[&t,2]*bootsample[&t,5];
						%if &c^= %then %do;
*/CONDITIONAL CDE*/;
x1=exp((theta1+theta3*&m)*(&a1-&a0));
bootsample[&t,7]=x1;
*/CONDITIONAL NDE*/;
bootsample[&t,8]=exp(theta1*(&a1-&a0))*
(1+exp(theta2+theta3*&a1+beta0+beta1*&a0+sum(beta2*t(c))))/
(1+exp(theta2+theta3*&a0+beta0+beta1*&a0+sum(beta2*t(c))));
*/CONDITIONAL NIE*/;
bootsample[&t,9]=(
(1+exp(beta0+beta1*&a0+sum(beta2*t(c))))*
(1+exp(theta2+theta3*&a0+beta0+beta1*&a1+sum(beta2*t(c))))
)/
(
(1+exp(beta0+beta1*&a1+sum(beta2*t(c))))*
(1+exp(theta2+theta3*&a0+beta0+beta1*&a0+sum(beta2*t(c))))
);
*/CONDITIONAL TNDE*/;
bootsample[&t,10]=exp(theta1*(&a1-&a0))*
(1+exp(theta2+theta3*&a1+beta0+beta1*&a1+sum(beta2*t(c))))/
(1+exp(theta2+theta3*&a0+beta0+beta1*&a1+sum(beta2*t(c))));
*/ CONDITIONAL TNIE*/;
bootsample[&t,11]=(
(1+exp(beta0+beta1*&a0+sum(beta2*t(c))))*
(1+exp(theta2+theta3*&a1+beta0+beta1*&a1+sum(beta2*t(c))))
)/
(
(1+exp(beta0+beta1*&a1+sum(beta2*t(c))))*
(1+exp(theta2+theta3*&a1+beta0+beta1*&a0+sum(beta2*t(c))))
);
bootsample[&t,12]=bootsample[&t,8]*bootsample[&t,11];
						%end;
					%end;
				%end;*t loop;
				%if &c^= %then %do;
x=bootsample;
cname1 = { "boot1" "boot2" "boot3" "boot4" "boot5" "boot6" "boot7" "boot8" "boot9" "boot10" "boot11" "boot12"};
create bootdata from x [ colname=cname1 ];
append from x;
				%end;
				%if &c= %then %do;
x=bootsample;
cname1 = { "boot1" "boot2" "boot3" "boot4" "boot5" "boot6"};
create bootdata from x [ colname=cname1 ];
append from x;
				%end;
			%end;

			%if &cvar= %then %do;

				%do t=1 %to &n;


USE out1&t;
READ ALL INTO VB;
NVB1= NROW(VB);
V1=VB[2:NVB1,];
theta1=VB[1,2];
theta2=VB[1,3];
					%if &interaction=true %then %do;
theta3=VB[1,4] ;
					%end;
USE out2&t;
READ ALL INTO VB;
NVB2= NROW(VB);
V2=VB[2:NVB2,];
					%if (&yreg=linear & &mreg=linear) | (&mreg=logistic) %then %do;
beta0=VB[1,1];
beta1=VB[1,2];
					%end;
					%if &yreg^=linear & &mreg=linear %then %do;
s2=VB[1,1];
s2=s2**2;
beta0=VB[1,2];
beta1=VB[1,3];
tsq=(theta3**2);
rm=s2;
asq=(&a1**2);
a1sq=(&a0**2);
					%end;
					%if &yreg=linear & &mreg=linear & &interaction=true %then %do;
*/CONDITIONAL=MARGINAL CDE*/;
bootsample[&t,1]=(theta1)*(&a1-&a0)+(theta3*(&m))*(&a1-&a0);
*/CONDITIONAL=MARGINAL NDE*/;
bootsample[&t,2]=(theta1+theta3*beta0+theta3*beta1*&a0)*(&a1-&a0);
*/CONDITIONAL=MARGINAL NIE*/;
bootsample[&t,3]=(theta2*beta1+theta3*beta1*&a0)*(&a1-&a0);
*/CONDITIONAL=MARGINAL TNDE*/;
bootsample[&t,4]=(theta1+theta3*beta0+theta3*beta1*&a1)*(&a1-&a0);
*/ CONDITIONAL=MARGINAL TNIE*/;
bootsample[&t,5]=(theta2*beta1+theta3*beta1*&a1)*(&a1-&a0);
*/te*/;
bootsample[&t,6]=(theta1+theta3*beta0+theta3*beta1*&a0+theta2*beta1+theta3*beta1*&a1)*(&a1-&a0);
					%end;
					%if &yreg=linear & &mreg=logistic & &interaction=true %then %do;
bootsample[&t,1]=(theta1+theta3*&m)*(&a1-&a0);
*/CONDITIONAL=MARGINAL NDE*/;
bootsample[&t,2]=(theta1+theta3*exp(beta0+beta1*&a0)/(1+exp(beta0+beta1*&a0)))*(&a1-&a0);
*/ CONDITIONAL=MARGINAL TNIE*/;
bootsample[&t,3]=(theta2+theta3*&a0)*(exp(beta0+beta1*&a1)/(1+exp(beta0+beta1*&a1))-exp(beta0+beta1*&a0)/(1+exp(beta0+beta1*&a0)));
*/CONDITIONAL=MARGINAL TNDE*/;
bootsample[&t,4]=(theta1+theta3*exp(beta0+beta1*&a1)/(1+exp(beta0+beta1*&a1)))*(&a1-&a0);
*/ CONDITIONAL=MARGINAL TNIE*/;
bootsample[&t,5]=(theta2+theta3*&a1)*(exp(beta0+beta1*&a1)/(1+exp(beta0+beta1*&a1))-exp(beta0+beta1*&a0)/(1+exp(beta0+beta1*&a0)));
*/te*/;
bootsample[&t,6]=bootsample[&t,2]+bootsample[&t,5];
					%end;
					%if &yreg^=linear & &mreg=linear & &interaction=true %then %do;
			 */MARGINAL=CONDITIONAL CDE*/;
x1=(theta1+theta3*&m)*(&a1-&a0);
bootsample[&t,1]=exp(x1);
      */MARGINAL=CONDITIONAL NDE*/;
x2=(theta1+theta3*beta0+theta3*beta1*&a0+theta3*theta2*rm)*(&a1-&a0)+(1/2)*tsq*rm*(asq-a1sq);
bootsample[&t,2]=exp(x2);
      */MARGINAL=CONDITIONAL NIE*/;
x3=(theta2*beta1+theta3*beta1*&a0)*(&a1-&a0);
bootsample[&t,3]=exp(x3);
      */MARGINAL=CONDITIONAL TNDE*/;
x4=(theta1+theta3*beta0+theta3*beta1*&a1+theta3*theta2*rm)*(&a1-&a0)+(1/2)*tsq*rm*(asq-a1sq);
      bootsample[&t,4]=exp(x4);
      */ MARGINAL=CONDITIONAL TNIE*/;
x5=(theta2*beta1+theta3*beta1*&a1)*(&a1-&a0);
bootsample[&t,5]=exp(x5);
 */ MARGINAL=CONDITIONAL TE*/;
bootsample[&t,6]=bootsample[&t,2]*bootsample[&t,5];
					%end;
					%if &yreg^=linear & &mreg=logistic & &interaction=false %then %do;
	  */MARGINAL=CONDITIONAL CDE*/;
	  x1=exp(theta1*(&a1-&a0));
bootsample[&t,1]=x1;
      */MARGINAL=CONDITIONAL NDE*/;
bootsample[&t,2]=exp((theta1)*(&a1-&a0))*
(1+exp(theta2+beta0+beta1*&a0))/
(1+exp(theta2+beta0+beta1*&a0));
      */MARGINAL=CONDITIONAL NIE*/;
bootsample[&t,3]=(
(1+exp(beta0+beta1*&a0))*
(1+exp(theta2+beta0+beta1*&a1))
)/
(
(1+exp(beta0+beta1*&a1))*
(1+exp(theta2+beta0+beta1*&a0))
);
      */MARGINAL=CONDITIONAL TNDE*/;
bootsample[&t,4]=exp((theta1)*(&a1-&a0))*
(1+exp(theta2+beta0+beta1*&a1))/
(1+exp(theta2+beta0+beta1*&a1));
      */ MARGINAL=CONDITIONAL TNIE*/;
bootsample[&t,5]=(
(1+exp(beta0+beta1*&a0))*
(1+exp(theta2+beta0+beta1*&a1))
)/
(
(1+exp(beta0+beta1*&a1))*(1+exp(theta2+beta0+beta1*&a0))
);
bootsample[&t,6]=bootsample[&t,2]*bootsample[&t,5];
					%end;
					%if &yreg^=linear & &mreg=logistic & &interaction=true %then %do;
*/MARGINAL CDE*/;
x6=(theta1+theta3*&m)*(&a1-&a0);
bootsample[&t,1]=exp(x6);
*/MARGINAL NDE*/;
bootsample[&t,2]=exp(theta1*(&a1-&a0))*
(1+exp(theta2+theta3*&a1+beta0+beta1*&a0))/
(1+exp(theta2+theta3*&a0+beta0+beta1*&a0));
*/MARGINAL NIE*/;
bootsample[&t,3]=(
(1+exp(beta0+beta1*&a0))*
(1+exp(theta2+theta3*&a0+beta0+beta1*&a1))
)/
(
(1+exp(beta0+beta1*&a1))*
(1+exp(theta2+theta3*&a0+beta0+beta1*&a0))
);
*/ MARGINAL TNDE*/;
bootsample[&t,4]=exp(theta1*(&a1-&a0))*
(1+exp(theta2+theta3*&a1+beta0+beta1*&a1))/
(1+exp(theta2+theta3*&a0+beta0+beta1*&a1));
*/ MARGINAL TNIE*/;
bootsample[&t,5]=(
(1+exp(beta0+beta1*&a0))*
(1+exp(theta2+theta3*&a1+beta0+beta1*&a1))
)/
(
(1+exp(beta0+beta1*&a1))*
(1+exp(theta2+theta3*&a1+beta0+beta1*&a0))
);
bootsample[&t,6]=bootsample[&t,2]*bootsample[&t,5];
					%end;
				%end;*t loop;

x=bootsample;
cname1 = { "boot1" "boot2" "boot3" "boot4" "boot5" "boot6"};
create bootdata from x [ colname=cname1 ];
append from x;
			%end;
		%end;*end linear linear int;

***************** causal effects for bootstrap END  *************************;

***************** causal effects, standard errors and confidence intervals from bootstrap *************************;

*/ no interaction ;
		%if (&mreg=linear & &interaction=false )| (&yreg=linear & &mreg=logistic & &interaction=false & &cvar=)  %then %do;

		*effects*;
use bootdata;
read all into bootdata;
effect=J(1,3);
			%do j=1 %to 3;

effect[,&j]=sum((bootdata[,&j]))/&n;

			%end;
x=(effect);
cname1 = {"effect1" "effect2" "effect3"};
create effect from x [colname=cname1] ;
append from x;
use bootdata;
read all into bootdata;
use effect;
read all into effect;
se=J(3,1);
square=J(&n,3);


*standard errors*;
			%do j=1 %to 3;

				%do t=1 %to &n;
square[&t,&j]=((bootdata[&t,&j])-effect[,&j])**2;
				%end;*t loop;

				se[&j,]=(
sqrt(
sum(
(square[,&j])
)
))/sqrt(&n);

			%end;
y=se;
create se from y;
append from y;
quit;
*Percentile confidence intervals*;

			%let alphalev = .05;
			%let a1 = %sysevalf(&alphalev/2*100);
			%let a2 = %sysevalf((1 - &alphalev/2)*100);

			%do j=1 %to 3;
proc univariate data = bootdata alpha = .05 noprint;
var boot&j;
output out=pmethod&j mean = effect&j pctlpts=&a1 &a2 pctlpre = p pctlname = _cil&j _ciu&j ;
run;
			%end;

proc iml;
			%do j=1 %to 3;
use pmethod&j;
read all into vb;
cil&j=vb[1,2];
ciu&j=vb[1,3];
			%end;
cil=cil1||cil2||cil3;
ciu=ciu1||ciu2||ciu3;
x= t(cil)||t(ciu) ;
create ci from x;
append from x;
quit;
		%end;*end  noint;




*/effects, standard errors, confidence intervals and p-value interaction*/;


*/ other ;

		%if (&interaction=true) | (&mreg=logistic & &interaction=false & &cvar^=) |(&yreg^=linear & &mreg=logistic & &interaction=false) %then %do;

use bootdata;
read all into bootdata;
			%if &c^= & &cvar^= %then %do;
effect=J(1,12);
				%do j=1 %to 12;

effect[,&j]=sum((bootdata[,&j]))/&n;

				%end;
x=effect;
cname1 = {"effect1" "effect2" "effect3" "effect4" "effect5" "effect6" "effect7" "effect8" "effect9" "effect10" "effect11" "effect12"};
create effect from x [colname=cname1] ;
append from x;
use bootdata;
read all into bootdata;
use effect;
read all into effect;
se=J(12,1);

square=J(&n,12);

*standard errors*;
				%do j=1 %to 12;
					%do t=1 %to &n;

square[&t,&j]=((bootdata[&t,&j])-effect[,&j])**2;
					%end;*t loop;


				se[&j,]=(
sqrt(
sum(
(square[,&j])
)
))/sqrt(&n);

				%end;
y=se;
create se from y;
append from y;
			%end;


			%if &cvar= | (&cvar^= & &c=) %then %do;
effect=J(1,6);
				%do j=1 %to 6;

effect[,&j]=sum((bootdata[,&j]))/&n;

				%end;
x=effect;
cname1 = {"effect1" "effect2" "effect3" "effect4" "effect5" "effect6"};
create effect from x [colname=cname1] ;
append from x;
use bootdata;
read all into bootdata;
use effect;
read all into effect;
se=J(6,1);

square=J(&n,6);

*standard errors*;
				%do j=1 %to 6;

					%do t=1 %to &n;

square[&t,&j]=((bootdata[&t,&j])-effect[,&j])**2;
					%end;*t loop;


				se[&j,]=(
sqrt(
sum(
(square[,&j])
)
))/sqrt(&n);

				%end;
y=se;
create se from y;
append from y;
quit;
			%end;


*Percentile confidence intervals*;
			%let alphalev = .05;
			%let a1 = %sysevalf(&alphalev/2*100);
			%let a2 = %sysevalf((1 - &alphalev/2)*100);

			%if &c^= %then %do;
		  	  %do j=1 %to 12;
proc univariate data = bootdata alpha = .05 noprint;
  var boot&j;
  output out=pmethod&j mean = effect&j pctlpts=&a1 &a2 pctlpre = p pctlname = _cil&j _ciu&j ;
run;
				%end;

			%end;
			%if (&cvar^= & &c=) | &cvar= %then %do;
				%do j=1 %to 6;
proc univariate data = bootdata alpha = .05 noprint;
var boot&j;
output out=pmethod&j mean = effect&j pctlpts=&a1 &a2 pctlpre = p pctlname = _cil&j _ciu&j ;
run;
				%end;
			%end;



proc iml;
			%if &c^= & &cvar^=  %then %do;

				%do j=1 %to 12;
USE pmethod&j;
read all into vb;
cil&j = vb[1,2];
ciu&j =vb[1,3] ;

				%end;
cil=cil1||cil2||cil3 ||cil4||cil5||cil6||cil7||cil8||cil9||cil10||cil11||cil12;
ciu=ciu1||ciu2||ciu3||ciu4||ciu5||ciu6||ciu7||ciu8||ciu9||ciu10||ciu11||ciu12;
x= t(cil)||t(ciu) ;
create ci from x;
append from x;

			%end;
			%if &cvar= | (&cvar^= & &c=) %then %do;

				%do j=1 %to 6;
USE pmethod&j;
read all into vb;
cil&j = vb[1,2];
ciu&j = vb[1,3];

				%end; *j;
cil=cil1||cil2||cil3 ||cil4||cil5||cil6;
ciu=ciu1||ciu2||ciu3||ciu4||ciu5||ciu6;
x= t(cil)||t(ciu) ;
create ci from x;
append from x;

			%end; *c;
		%end; *yreg etc;
	%end; * boot;

***************************   BOOTSTRAP PROCEDURE -END-  ***************************************************************;

dm "out;clear";


************* regression to print **************************************************************************************;


************************************************************************************************************************;
	%if &yreg=linear %then %do;
************************************************************************************************************************;

		%if &interaction=false & &cvar^= %then %do;
proc reg data=data1 covout
outest=out1(drop=_model_ _type_ _name_ _depvar_ _rmse_ &yvar) ;
model  &yvar=&avar &mvar &cvar ;
proc print;
run;
		%end;
		%if &interaction=false & &cvar= %then %do;
proc reg data=data1 covout
outest=out1(drop=_model_ _type_ _name_ _depvar_ _rmse_ &yvar) ;
model  &yvar=&avar &mvar;
proc print;
run;
		%end;
		%if &interaction=true & &cvar^= %then %do;
proc reg data=data1 covout
outest=out1(drop=_model_ _type_ _name_ _depvar_ _rmse_ &yvar) ;
model  &yvar=&avar &mvar int &cvar ;
proc print;
run;
		%end;
		%if &interaction=true & &cvar= %then %do;
proc reg data=data1 covout
outest=out1(drop=_model_ _type_ _name_ _depvar_ _rmse_ &yvar) ;
model  &yvar=&avar &mvar int ;
proc print;
run;
		%end;
	%end;

***********************************************************************************************************************************************************************************************************************************;
	%if &yreg=logistic  | &yreg=loglinear  |&yreg=poisson | &yreg=negbin |&yreg=survCox |&yreg=survAdd |&yreg=survAFT_weibull |&yreg=survAFT_exp |&yreg=survAFT_gamma |&yreg=survAFT_loglogistic  |&yreg=survAFT_normal %then %do;
***********************************************************************************************************************************************************************************************************************************;

*need to include output to print for survival outcome!!!;


		%if &interaction=true & &cvar^= %then %do;
			%if &yreg=logistic %then %do;
proc logistic  data=data1 descending covout
outest=out1(drop=_link_ _status_ _LNLIKE_ _type_ _name_ )  ;
model  &yvar=&avar &mvar int &cvar ;
run;
			%end;
			%if &yreg=loglinear %then %do;
proc genmod data=data1 descending;
model &yvar=&avar &mvar  int &cvar/dist=binomial link=log covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[,];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-1,2];
par=t(par);
out=par//cov;
create out1 from out;
append from out;
quit;
			%end;
			%if &yreg=poisson %then %do;
proc genmod data=data1 ;
model &yvar=&avar &mvar  int &cvar/dist=poisson covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[,];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-1,2];
par=t(par);
out=par//cov;
create out1 from out;
append from out;
quit;
			%end;
			%if &yreg=negbin %then %do;
proc genmod data=data1 ;
model &yvar=&avar &mvar  int &cvar/dist=negbin covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[1:nrow(vb)-1,1:ncol(vb)-1];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-1,2];
par=t(par);
out=par//cov;
create out1 from out;
append from out;
quit;
			%end;

*/add survival*/;
*/cen and non cens case*/;
*/note COX does not have intercept!!!! I add a column of zeros*/;

				%if &yreg=survCox %then %do;
	proc phreg data = data1;
model &yvar*&cens(1) = &avar &mvar int &cvar/covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
			run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[,];
use gmparms;
read all into vb;
par=vb[1:nrow(vb),2];
par=t(par);
x=vb[1:nrow(vb),2]*0//0//0;
rzero=t(vb[1:nrow(vb),2]*0);
out=par//rzero//cov;
out=x||out;
create out1 from out;
append from out;
quit;

%end;

%if &yreg=survAFT_exp %then %do;
proc lifereg data = data1;
   model &yvar*&cens(1) = &avar &mvar int &cvar/covb DISTRIBUTION=exponential ;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
			run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[1:nrow(vb)-1,1:nrow(vb)-1];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-2,2];
par=t(par);
out=par//cov;
create out1 from out;
append from out;
quit;
%end;

%if &yreg=survAFT_weibull %then %do;
proc lifereg data = data1;
   model &yvar*&cens(1) = &avar &mvar int &cvar/covb DISTRIBUTION=weibull ;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
			run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[1:nrow(vb)-1,1:nrow(vb)-1];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-2,2];
par=t(par);
out=par//cov;
create out1 from out;
append from out;
quit;
%end;

		%end;

		%if &interaction=true & &cvar= %then %do;
			%if &yreg=logistic %then %do;
proc logistic  data=data1 descending covout
outest=out1(drop=_link_ _status_ _LNLIKE_ _type_ _name_ )  ;
model  &yvar=&avar &mvar int ;
run;
			%end;
			%if &yreg=loglinear %then %do;
proc genmod data=data1 descending;
model &yvar=&avar &mvar int /dist=binomial link=log covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[,];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-1,2];
par=t(par);
out=par//cov;
create out1 from out;
append from out;
quit;
			%end;
			%if &yreg=poisson %then %do;
proc genmod data=data1;
model &yvar=&avar &mvar  int /dist=poisson covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[,];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-1,2];
par=t(par);
out=par//cov;
create out1 from out;
append from out;
quit;
			%end;
			%if &yreg=negbin %then %do;
proc genmod data=data1 ;
model &yvar=&avar &mvar  int /dist=negbin covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[1:nrow(vb)-1,1:ncol(vb)-1];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-1,2];
par=t(par);
out=par//cov;
create out1 from out;
append from out;
quit;
			%end;
*/add survival*/;
*/cen and non cens case*/;
*/note COX does not have intercept!!!! I add a column of zeros*/;
	%if &yreg=survCox %then %do;
	proc phreg data = data1;
model &yvar*&cens(1) = &avar &mvar int /covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
			run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[,];
use gmparms;
read all into vb;
par=vb[1:nrow(vb),2];
par=t(par);
x=vb[1:nrow(vb),2]*0//0//0;
rzero=t(vb[1:nrow(vb),2]*0);
out=par//rzero//cov;
out=x||out;
create out1 from out;
append from out;
quit;

%end;

%if &yreg=survAFT_exp %then %do;
proc lifereg data = data1;
   model &yvar*&cens(1) = &avar &mvar int /covb DISTRIBUTION=exponential ;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
			run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[1:nrow(vb)-1,1:nrow(vb)-1];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-2,2];
par=t(par);
out=par//cov;
create out1 from out;
append from out;
quit;
%end;

%if &yreg=survAFT_weibull %then %do;
proc lifereg data = data1;
   model &yvar*&cens(1) = &avar &mvar int /covb DISTRIBUTION=weibull ;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
			run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[1:nrow(vb)-1,1:nrow(vb)-1];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-2,2];
par=t(par);
out=par//cov;
create out1 from out;
append from out;
quit;
%end;
		%end;

		%if &interaction=false & &cvar^= %then %do;
			%if &yreg=logistic %then %do;
proc logistic  data=data1 descending covout
outest=out1(drop=_link_ _status_ _LNLIKE_ _type_ _name_ )  ;
model  &yvar=&avar &mvar  &cvar ;
run;
			%end;
			%if &yreg=loglinear %then %do;
proc genmod data=data1 descending;
model &yvar=&avar &mvar  &cvar/dist=binomial link=log covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[,];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-1,2];
par=t(par);
out=par//cov;
create out1 from out;
append from out;
quit;
			%end;
			%if &yreg=poisson %then %do;
proc genmod data=data1 ;
model &yvar=&avar &mvar  &cvar/dist=poisson covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[,];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-1,2];
par=t(par);
out=par//cov;
create out1 from out;
append from out;
quit;
			%end;
			%if &yreg=negbin %then %do;
proc genmod data=data1 ;
model &yvar=&avar &mvar  &cvar/dist=negbin covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[1:nrow(vb)-1,1:ncol(vb)-1];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-1,2];
par=t(par);
out=par//cov;
create out1 from out;
append from out;
quit;
			%end;
*/add survival*/;
*/cen and non cens case*/;
*/note COX does not have intercept!!!! I add a column of zeros*/;
	%if &yreg=survCox %then %do;
	proc phreg data = data1;
model &yvar*&cens(1) = &avar &mvar &cvar/covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
			run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[,];
use gmparms;
read all into vb;
par=vb[1:nrow(vb),2];
par=t(par);
x=vb[1:nrow(vb),2]*0//0//0;
rzero=t(vb[1:nrow(vb),2]*0);
out=par//rzero//cov;
out=x||out;
create out1 from out;
append from out;
quit;

%end;

%if &yreg=survAFT_exp %then %do;
proc lifereg data = data1;
   model &yvar*&cens(1) = &avar &mvar &cvar/covb DISTRIBUTION=exponential ;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
			run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[1:nrow(vb)-1,1:nrow(vb)-1];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-2,2];
par=t(par);
out=par//cov;
create out1 from out;
append from out;
quit;
%end;

%if &yreg=survAFT_weibull %then %do;
proc lifereg data = data1;
   model &yvar*&cens(1) = &avar &mvar &cvar/covb DISTRIBUTION=weibull ;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
			run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[1:nrow(vb)-1,1:nrow(vb)-1];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-2,2];
par=t(par);
out=par//cov;
create out1 from out;
append from out;
quit;
%end;
		%end;

		%if &interaction=false & &cvar= %then %do;
			%if &yreg=logistic %then %do;
proc logistic  data=data1 descending covout
outest=out1(drop=_link_ _status_ _LNLIKE_ _type_ _name_ )  ;
model  &yvar=&avar &mvar  ;
run;
			%end;
			%if &yreg=loglinear %then %do;
proc genmod data=data1 descending;
model &yvar=&avar &mvar /dist=binomial link=log covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[,];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-1,2];
par=t(par);
out=par//cov;
create out1 from out;
append from out;
quit;
			%end;
			%if &yreg=poisson %then %do;
proc genmod data=data1 ;
model &yvar=&avar &mvar /dist=poisson covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[,];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-1,2];
par=t(par);
out=par//cov;
create out1 from out;
append from out;
quit;
			%end;
			%if &yreg=negbin %then %do;
proc genmod data=data1 ;
model &yvar=&avar &mvar /dist=negbin covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[1:nrow(vb)-1,1:ncol(vb)-1];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-1,2];
par=t(par);
out=par//cov;
create out1 from out;
append from out;
quit;
			%end;
*/add survival*/;
*/cen and non cens case*/;
*/note COX does not have intercept!!!! I add a column of zeros*/;
	%if &yreg=survCox %then %do;
	proc phreg data = data1;
model &yvar*&cens(1) = &avar &mvar/covb;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
			run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[,];
use gmparms;
read all into vb;
par=vb[1:nrow(vb),2];
par=t(par);
x=vb[1:nrow(vb),2]*0//0//0;
rzero=t(vb[1:nrow(vb),2]*0);
out=par//rzero//cov;
out=x||out;
create out1 from out;
append from out;
quit;

%end;

%if &yreg=survAFT_exp %then %do;
proc lifereg data = data1;
   model &yvar*&cens(1) = &avar &mvar /covb DISTRIBUTION=exponential ;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
			run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[1:nrow(vb)-1,1:nrow(vb)-1];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-2,2];
par=t(par);
out=par//cov;
create out1 from out;
append from out;
quit;
%end;

%if &yreg=survAFT_weibull %then %do;
proc lifereg data = data1;
   model &yvar*&cens(1) = &avar &mvar /covb DISTRIBUTION=weibull ;
ods output ParameterEstimates=gmparms
			CovB=gmcovb;
			run;
proc iml;
use gmcovb;
read all into vb;
cov=vb[1:nrow(vb)-1,1:nrow(vb)-1];
use gmparms;
read all into vb;
par=vb[1:nrow(vb)-2,2];
par=t(par);
out=par//cov;
create out1 from out;
append from out;
quit;
%end;
		%end;
	%end;

************************************************************************************************************************;
	%if  &mreg=linear & &yreg=linear %then %do;
************************************************************************************************************************;

		%if &cvar^= %then %do;
			%if &casecontrol^=true %then %do;
proc reg data=data1 covout
outest=out2(drop=_model_ _type_ _name_ _depvar_ _rmse_ &mvar)  ;
model  &mvar=&avar &cvar ;
proc print;
run;
			%end;
			%if &casecontrol=true %then %do;
proc reg data=data1 covout
outest=out2(drop=_model_ _type_ _name_ _depvar_ _rmse_ &mvar)  ;
where &yvar=0;
model  &mvar=&avar &cvar;
proc print;
run;
			%end;
		%end;
		%if  &cvar= %then %do;
			%if &casecontrol^=true %then %do;
proc reg data=data1 covout
outest=out2(drop=_model_ _type_ _name_ _depvar_ _rmse_ &mvar)  ;
model  &mvar=&avar ;
proc print;
run;
			%end;
			%if &casecontrol=true %then %do;
proc reg data=data1 covout
outest=out2(drop=_model_ _type_ _name_ _depvar_ _rmse_ &mvar)  ;
where &yvar=0;
model  &mvar=&avar;
proc print;
run;
			%end;
		%end;
	%end;


************************************************************************************************************************;
	%if  &mreg=linear & &yreg^=linear %then %do;
************************************************************************************************************************;

		%if &cvar^= %then %do;
			%if &casecontrol^=true %then %do;
proc reg data=data1 covout
outest=out2(drop=_model_ _type_ _name_ _depvar_ &mvar)  ;
model  &mvar=&avar &cvar ;
proc print;
run;
			%end;
			%if &casecontrol=true %then %do;
proc reg data=data1 covout
outest=out2(drop=_model_ _type_ _name_ _depvar_  &mvar)  ;
where &yvar=0;
model  &mvar=&avar &cvar;
proc print;
run;
			%end;
		%end;
		%if  &cvar= %then %do;
			%if &casecontrol^=true %then %do;
proc reg data=data1 covout
outest=out2(drop=_model_ _type_ _name_ _depvar_  &mvar)  ;
model  &mvar=&avar ;
proc print;
run;
			%end;
			%if &casecontrol=true %then %do;
proc reg data=data1 covout
outest=out2(drop=_model_ _type_ _name_ _depvar_ &mvar)  ;
where &yvar=0;
model  &mvar=&avar;
proc print;
run;
			%end;
		%end;
	%end;

************************************************************************************************************************;
	%if &mreg=logistic %then %do;
************************************************************************************************************************;

		%if &cvar^= %then %do;
			%if &casecontrol^=true %then %do;
proc logistic data=data1 descending covout
outest=out2(drop=_link_ _status_ _LNLIKE_ _type_ _name_ )  ;
model  &mvar=&avar &cvar ;
run;
			%end;
			%if &casecontrol=true %then %do;
proc logistic data=data1 descending covout
outest=out2(drop=_link_ _status_ _LNLIKE_ _type_ _name_ )  ;
where &yvar=0;
model  &mvar=&avar &cvar ;
run;
			%end;
		%end;
		%if &cvar= %then %do;
			%if &casecontrol^=true %then %do;
proc logistic data=data1 descending covout
outest=out2(drop=_link_ _status_ _LNLIKE_ _type_ _name_ )  ;
model  &mvar=&avar  ;
run;
			%end;
			%if &casecontrol=true %then %do;
proc logistic data=data1 descending covout
outest=out2(drop=_link_ _status_ _LNLIKE_ _type_ _name_ )  ;
where &yvar=0;
model  &mvar=&avar  ;
run;
			%end;
		%end;
	%end;

**********************************************************************************************************************;

***************************   DELTA METHOD PROCEDURE  ***************************************************************;

*/PROBLEMS WITH COX REGRESSION!!!! NEED TO DEBUG! */;


	%if (&boot= | &boot=false) %then %do;

proc iml;

*/compute the causal effects*/;
		%if (&mreg=linear & &interaction=false ) | (&yreg=linear & &mreg=logistic & &interaction=false & &cvar=)  %then %do;
PROC IML;
USE out1;
READ ALL INTO VB;
NVB1= NROW(VB);
V1=VB[2:NVB1,];
theta1=VB[1,2];
theta2=VB[1,3];
USE out2;
READ ALL INTO VB;
			%if (&yreg=linear & &mreg=linear) | (&yreg=linear & &mreg=logistic) %then %do;
beta0= VB[1,1];
beta1=VB[1,2];
NVB2= NROW(VB);
V2=VB[2:NVB2,];
zero1=J(nrow(V1),nrow(V2),0);
zero2=J(nrow(V2),nrow(V1),0);
A= V2 || zero2;
B= zero1 || V1;
sigma= A // B;
zero=0;
one=1;		%end;
			%if (&yreg^=linear & &mreg=linear) %then %do;
beta1=VB[1,3];
NVB2= NROW(VB);
V2=VB[2:NVB2,2:ncol(vb)];
zero1=J(nrow(V1),nrow(V2),0);
zero2=J(nrow(V2),nrow(V1),0);
A= V2 || zero2;
B= zero1 || V1;
sigma= A // B;
zero=J(1,1,0);
one=J(1,1,1);

			%end;
effect=J(1,3);
			%if &cvar^= %then %do;
		z1=J(1,&nc,0);
z=zero||zero||z1||zero||one||zero;

gamma=J(3,2*&nc+5);
			%end;
			%if &cvar= %then %do;
gamma=J(3,5);
			%end;
			%if (&yreg=linear & &mreg=logistic) %then %do;
z=zero||zero||zero||one||zero;
*/cde and nde*/;
effect[,1]=(theta1)*(&a1-&a0);
gamma[1,]=z;
*/nie*/;
effect[,2]=(theta2)*(exp(beta0+beta1*&a1)/(1+exp(beta0+beta1*&a1))-exp(beta0+beta1*&a0)/(1+exp(beta0+beta1*&a0)));
D=exp(beta0+beta1*&a1);
E=(1+D);
A=exp(beta0+beta1*&a0);
B=(1+A);
x=(theta2)*((D*E-D**2)/E**2-(A*B-A**2)/B**2);
w=(theta2)*(&a1*(D*E-D**2)/E**2-&a0*(A*B-A**2)/B**2);
h=t(D/E-A/B);
gamma[2,]=x|| w|| zero||zero|| h;
*/te*/;
effect[,3]=effect[,1]+effect[,2];
A=exp(beta0+beta1*&a0);
B=(1+A);
D=exp(beta0+beta1*&a1);
E=(1+D);
/* Changed by @kaz-yos on 2020-04-01 based on VV2013 Appendix p14. */
/* The squared terms in the numerator is different from the denominator. */
/* Corrected from
x=(theta2)*((D*E-D**2)/(E**2)-(A*B-B**2)/(B**2));*/
x=(theta2)*((D*E-D**2)/(E**2)-(A*B-A**2)/(B**2));
/* Corrected from
w=(theta2)*(&a1*(D*E-D**2)/(E**2)-&a0*(A*B-B**2)/(B**2)); */
w=(theta2)*(&a1*(D*E-D**2)/(E**2)-&a0*(A*B-A**2)/(B**2));
t=t(D/E-A/B);
s=(&a1-&a0);
gamma[3,]=x||w||zero||s||t;
			%end;
			%if (&yreg=linear & &mreg=linear & &interaction=false) %then %do;
		effect[,1]=(theta1)*(&a1-&a0);
		effect[,2]=(theta2*beta1)*(&a1-&a0);
		effect[,3]=(theta1+theta2*beta1)*(&a1-&a0);
				%if &cvar^= %then %do;
z1=J(1,&nc,0);
z=zero||zero||z1||zero||one||zero;
gamma[1,]=zero||zero||z1||zero||one||zero||z1;
gamma[2,]=zero|| theta2||z1|| zero ||zero|| beta1||z1;
gamma[3,]=zero||theta2||z1||zero||one||beta1||z1;
				%end;
				%if &cvar= %then %do;
z=zero||zero||zero||one||zero;
gamma[1,]=zero||zero||zero||one||zero;
gamma[2,]=zero|| theta2|| zero ||zero|| beta1;
gamma[3,]=zero||theta2||zero||one||beta1;
				%end;
			%end;
			%if (&yreg^=linear & &mreg=linear ) %then %do;
effect[,1]=((theta1)*(&a1-&a0));
effect[,2]=((theta2*beta1)*(&a1-&a0));
effect[,3]=effect[,2]+effect[,1];
				%if &cvar^= %then %do;
z1=J(1,&nc,0);
z=zero||zero||z1||zero||one||zero;
gamma[1,]=zero||zero||z1||zero||one||zero||z1;
gamma[2,]=zero|| theta2||z1|| zero ||zero|| beta1||z1;
gamma[3,]=zero||theta2||z1||zero||one||beta1||z1;
				%end;
				%if &cvar= %then %do;
z=zero||zero||zero||one||zero;
gamma[1,]=zero||zero||zero||one||zero;
gamma[2,]=zero|| theta2|| zero ||zero|| beta1;
gamma[3,]=zero||theta2||zero||one||beta1;
				%end;
			%end;
se=J(1,3);
pvalue=J(1,3);
cil=J(1,3);
ciu=J(1,3);
				%if (&mreg=logistic & &yreg=linear) %then %do;
se[,1]=sqrt(gamma[1,]*sigma*t(gamma[1,]))*abs(&a1-&a0);
se[,2]=sqrt(gamma[2,]*sigma*t(gamma[2,]));
se[,3]=sqrt(gamma[3,]*sigma*t(gamma[3,]));
				%end;
				%if &mreg^=logistic  %then %do;
				%do j=1 %to 3;
se[,&j]=sqrt(gamma[&j,]*sigma*t(gamma[&j,]))*abs(&a1-&a0);
				%end;
				%end;
			%do j=1 %to 3;
pvalue[,&j] = 2*MIN(1-ABS(probnorm((effect[,&j])/(se[,&j]))),ABS(probnorm((effect[,&j])/(se[,&j]))));
				%if (&yreg=linear ) %then %do;
cil[,&j]=effect[,&j]-1.96*(se[,&j]);
ciu[, &j]=effect[,&j]+1.96*(se[,&j]);
				%end;
				%if (&yreg^=linear & &mreg=linear ) %then %do;

cil[,&j]=exp(effect[,&j]-1.96*(se[,&j]));
ciu[, &j]=exp(effect[,&j]+1.96*(se[,&j]));
effect[,&j]=exp(effect[,&j]);
				%end;
 			%END;
x=effect;
cname1 = { "effect1" "effect2" "effect3"};
create effect from x [colname=cname1];
append from x;
x=(se);
cname1 = { "se1" "se2" "se3"};
create se from x [ colname=cname1 ];
append from x;
x=cil;
cname1 = { "cil1" "cil2" "cil3"};
create cil from x [ colname=cname1 ];
append from x;
x=ciu;
cname1 = { "ciu1" "ciu2" "ciu3"};
create ciu from x [ colname=cname1 ];
append from x;
x=pvalue;
cname1 = { "p1" "p2" "p3"};
create pvalue from x [ colname=cname1 ];
append from x;
		%end;* noint;

***************************************************************************************;

		%if (&interaction=true) | (&mreg=logistic & &interaction=false & &cvar^=) |(&yreg^=linear & &mreg=logistic & &interaction=false) %then %do;

			%if &cvar^= %then %do;
USE data2;
read all into vb;
				%if &c= %then %do;
cmean=VB[1,1:ncol(vb)];
				%end;

				%if &c^= %then %do;
cmean=VB[1,1:ncol(vb)-&nc];
c=VB[1,ncol(vb)-&nc+1:ncol(vb)] ;
				%end;






				%if &c^= & &interaction=false %then %do;
effect=J(1,12);
gamma=J(12,&nc*2+5);
				%end;
				%if &c^= & &interaction=true %then %do;
effect=J(1,12);
gamma=J(12,&nc*2+6);

				%end;
				%if &c= & &interaction=false %then %do;
effect=J(1,6);
gamma=J(6,&nc*2+5);
				%end;
				%if &c= & &interaction=true %then %do;
effect=J(1,6);
gamma=J(6,&nc*2+6);
				%end;


USE out1;
READ ALL INTO VB;
NVB1= NROW(VB);
V1=VB[2:NVB1,];
theta1=VB[1,2];
theta2=VB[1,3];
				%if &interaction=true %then %do;
theta3=VB[1,4] ;
				%end;
USE out2;
READ ALL INTO VB;

				%if (&yreg=linear & &mreg=linear) | (&mreg=logistic) %then %do;
NVB2= NROW(VB);
V2=VB[2:NVB2,];
beta0=VB[1,1];
beta1=VB[1,2];
beta2= VB[1,3:ncol(vb)];
zero1=J(nrow(V1),nrow(V2),0);
zero2=J(nrow(V2),nrow(V1),0);
A= V2 || zero2;
B= zero1 || V1;
sigma= A // B;
zero=0;
one=1;
z1=J(1,&nc,0);
z=zero||zero||z1||zero||one||zero;
				%end;

				%if &yreg^=linear & &mreg=linear %then %do;
s2=VB[1,1];
s2=s2**2;
beta0=VB[1,2];
beta1=VB[1,3];
beta2= VB[1,4:ncol(vb)];
tsq=(theta3**2);
rm=s2;
asq=(&a1**2);
a1sq=(&a0**2);
NVB2= NROW(VB);
colvb=ncol(vb);
V2=VB[2:NVB2,2:colvb];
zero1=J(nrow(V1),nrow(V2),0);
zero2=J(nrow(V2),nrow(V1),0);
z2=J(nrow(V1),1,0);
z3=J(nrow(V2),1,0);
A= V2 || zero2 ||z3;
B= zero1 || V1||z2;
zeros=J(1,nrow(V1)+nrow(V2),0);
/* @kaz-yos on 2020-05-02 */
/* This is wrong. It should be the following based on V2015 p470.
where n = sample size, p = length(betas), s2 = sigma^2
D= zeros || ((2 * (s2**2)) / (n - p)) */
D= zeros ||s2;
sigma= A // B//D;
zero=0;
one=1;
z1=J(1,&nc,0);
z=zero||zero||z1||zero||one||zero;
				%if &c^= & &interaction=false %then %do;
gamma=J(12,&nc*2+6);
				%end;
				%if &c^= & &interaction=true %then %do;
gamma=J(12,&nc*2+7);
				%end;
				%if &c= & &interaction=false %then %do;
gamma=J(6,&nc*2+6);
				%end;
				%if &c= & &interaction=true %then %do;
gamma=J(6,&nc*2+7);
				%end;
				%end;


				%if &yreg=linear & &mreg=linear & &interaction=true %then %do;
*/MARGINAL CDE*/;
effect[,1]=(theta1+theta3*&m)*(&a1-&a0);
*/MARGINAL NDE*/;
effect[,2]=(theta1+theta3*beta0+theta3*beta1*&a0+(theta3*beta2*t(cmean)))*(&a1-&a0);
*/MARGINAL NIE*/;
effect[,3]=(theta2*beta1+theta3*beta1*&a0)*(&a1-&a0);
*/ MARGINAL TNDE*/;
effect[,4]=(theta1+theta3*beta0+theta3*beta1*&a1+(theta3*beta2*t(cmean)))*(&a1-&a0);
*/ MARGINAL TNIE*/;
effect[,5]=(theta2*beta1+theta3*beta1*&a1)*(&a1-&a0);
*/te marginal*/;
effect[,6]=(theta1+theta3*beta0+theta3*beta1*&a0+(theta3*beta2*t(cmean))+theta2*beta1+theta3*beta1*&a1)*(&a1-&a0);
z1=J(1,&nc,0);
zero=0;
one=1;
gamma[1,]=zero||zero||z1||zero||one||zero||&m || z1 ;
x1=theta3*&a0;
print x1;
w=theta3*t(cmean);
print w;
h1=beta0+beta1*&a0+(beta2)*t(cmean);
print h1;
gamma[2,]= theta3|| x1|| t(w) || zero|| one|| zero|| t(h1) ||z1;

x0=theta3*&a1;
w=theta3*t(cmean);
h0=beta0+beta1*&a1+(beta2)*t(cmean);
gamma[4,]=theta3|| x0|| t(w)|| zero|| one|| zero|| t(h0)||z1;
x0=theta2+theta3*&a1;                   /* This is the same as x1 except for a1. */
w0=beta1*&a1;
gamma[5,]=zero|| x0|| z1|| zero||zero|| beta1|| w0 || z1; /* m-tnie uses x1. */
x1=theta2+theta3*&a0;                   /* This is the same as x0 except for a0. */
w1=beta1*&a0;
gamma[3,]=zero|| x1|| z1|| zero||zero|| beta1|| w1 || z1; /* m-pnie uses x0. */
D=theta3*(cmean);
A=(theta3*&a1+theta3*&a0+theta2);
B=beta0+beta1*(&a1+&a0)+beta2*t(cmean);
gamma[6,]=theta3||A||(D)||zero||one||beta1||B||z1;

					%if &c^= %then %do;
*/CONDITIONAL CDE*/;
effect[,7]=(theta1)*(&a1-&a0)+(theta3*(&m))*(&a1-&a0);
*/CONDITIONAL NDE*/;
effect[,8]=(theta1+theta3*beta0+theta3*beta1*&a0+(theta3*beta2*t(c)))*(&a1-&a0);
*/CONDITIONAL NIE*/;
effect[,9]=(theta2*beta1+theta3*beta1*&a0)*(&a1-&a0);
*/CONDITIONAL TNDE*/;
effect[,10]=(theta1+theta3*beta0+theta3*beta1*&a1+(theta3*beta2*t(c)))*(&a1-&a0);
*/ CONDITIONAL TNIE*/;
effect[,11]=(theta2*beta1+theta3*beta1*&a1)*(&a1-&a0); /* cond tnie expression */
*/te conditional*/;
effect[,12]=(theta1+theta3*beta0+theta3*beta1*&a0+(theta3*beta2*t(c))+theta2*beta1+theta3*beta1*&a1)*(&a1-&a0);
*gamma=J(1,2*&nc+6);
gamma[7,]=zero||zero||z1||zero||one||zero||&m || z1;
x1=theta3*&a0;
w=theta3*t(c);
h1=beta0+beta1*&a0+(beta2)*t(c);
gamma[8,]= theta3|| x1|| t(w) || zero|| one|| zero|| t(h1) || z1;
x0=theta3*&a1;
h0=beta0+beta1*&a1+(beta2)*t(c);
gamma[10,]=theta3|| x0|| t(w)|| zero|| one|| zero|| t(h0)||z1;
/* The following line defining x0 was added by @kaz-yos on 2020-03-28 following V2015 p466.
This is the mreg linear, yreg linear, interaction true, cvar non-empty case.
This line was originally missing, causing gamma[11,] (Gamma for cond tnie) to refer to x0 defined
for gamma[10,] (Gamma for cond tnde). The second slot for gamma[11,] is the partial derivative of
effect[,11] expression above wrt beta1. So it should be (theta2 + theta3 * a1). */
x0=theta2+theta3*&a1;
w0=beta1*&a1;
gamma[11,]=zero|| x0|| z1|| zero||zero|| beta1|| w0 || z1; /* Gamma for cond tnie uses x0. */
x1=theta2+theta3*&a0;    /* This seems correct and the same as x1 for gamma[3,] */
w1=beta1*&a0;
gamma[9,]=zero|| x1|| z1|| zero||zero|| beta1|| w1 || z1;  /* Gamma for cond pnie uses x1.*/
D=theta3*(c);
A=(theta3*&a1+theta3*&a0+theta2);
B=beta0+beta1*(&a1+&a0)+beta2*t(c);
gamma[12,]=theta3||A||(D)||zero||one||beta1||B||z1;

					%end;
				%end;



				%if &yreg=linear & &mreg=logistic & &interaction=true  %then %do;
*/MARGINAL CDE*/;
effect[,1]=(theta1+theta3*&m)*(&a1-&a0);
*/MARGINAL NDE*/;
effect[,2]=(theta1+theta3*exp(beta0+beta1*&a0+sum(beta2*t(cmean)))/(1+exp(beta0+beta1*&a0+sum(beta2*t(cmean)))))*(&a1-&a0);
*/MARGINAL NIE*/;
effect[,3]=(theta2+theta3*&a0)*(
exp(beta0+beta1*&a1+sum(beta2*t(cmean)))/
(1+exp(beta0+beta1*&a1+sum(beta2*t(cmean))))-
exp(beta0+beta1*&a0+sum(beta2*t(cmean)))/
(1+exp(beta0+beta1*&a0+sum(beta2*t(cmean))))
);
*/ MARGINAL TNDE*/;
effect[,4]=(theta1+theta3*exp(beta0+beta1*&a1+sum(beta2*t(cmean)))/(1+exp(beta0+beta1*&a1+sum(beta2*t(cmean)))))*(&a1-&a0);
*/ MARGINAL TNIE*/;
effect[,5]=(theta2+theta3*&a1)*(exp(beta0+beta1*&a1+sum(beta2*t(cmean)))/(1+exp(beta0+beta1*&a1+sum(beta2*t(cmean))))-exp(beta0+beta1*&a0+sum(beta2*t(cmean)))/(1+exp(beta0+beta1*&a0+sum(beta2*t(cmean)))));
*/te marginal*/;
effect[,6]=effect[,2]+effect[,5];

gamma[1,]=z||&m || z1 ;
A=exp(beta0+beta1*&a0+beta2*t(cmean));
B=(1+A);
x=theta3*(A*B-A**2)/B**2;
w=theta3*&a0*(A*B-A**2)/B**2;
y=theta3*cmean*(A*B-A**2)/B**2;
h=A/B;
gamma[2,]=  x|| w || y || zero|| one|| zero|| t(h) || z1;
A=exp(beta0+beta1*&a1+beta2*t(cmean));
B=(1+A);
x=theta3*(A*B-A**2)/B**2;
w=theta3*&a1*(A*B-A**2)/B**2;
y=theta3*cmean*(A*B-A**2)/B**2;
h=A/B;
gamma[4,]=x|| w || y || zero|| one|| zero|| t(h) || z1;
D=exp(beta0+beta1*&a1+beta2*t(cmean));
E=(1+A);
A=exp(beta0+beta1*&a0+beta2*t(cmean));
B=(1+A);

x=(theta2+theta3*&a1)*((D*E-D**2)/E**2-(A*B-A**2)/B**2);
w=(theta2+theta3*&a1)*(&a1*(D*E-D**2)/E**2-&a0*(A*B-A**2)/B**2);
y=cmean*(theta2+theta3*&a1)*((D*E-D**2)/E**2-(A*B-A**2)/B**2);
h=t(D/E-A/B);
j=&a1*h;
gamma[5,]=x|| w|| y|| zero||zero|| h|| j || z1;
x=(theta2+theta3*&a0)*((D*E-D**2)/E**2-(A*B-A**2)/B**2);
w=(theta2+theta3*&a0)*(&a1*(D*E-D**2)/E**2-&a0*(A*B-A**2)/B**2);
y=cmean*(theta2+theta3*&a0)*((D*E-D**2)/E**2-(A*B-A**2)/B**2);
h=t(D/E-A/B);
j=&a0*h;
gamma[3,]=x|| w|| y||zero||zero|| h|| j|| z1;
A=exp(beta0+beta1*&a0+beta2*t(cmean));
B=(1+A);
D=exp(beta0+beta1*&a1+beta2*t(cmean));
E=(1+D);
/* Gamma_te corrected by @kaz-yos on 2020-04-01 based on VV2013 Appendix p14. */
/* The squared terms in the numerator should be different from the denominator. */
/* Corrected from
x=theta3*(&a1-&a0)*(A*B-B**2)/(B**2)+(theta2+theta3*&a1)*((D*E-D**2)/(E**2)-(A*B-B**2)/(B**2)); */
x=theta3*(&a1-&a0)*(A*B-A**2)/(B**2)+(theta2+theta3*&a1)*((D*E-D**2)/(E**2)-(A*B-A**2)/(B**2));
/* Corrected from
w=&a0*theta3*(&a1-&a0)*(A*B-B**2)/(B**2)+(theta2+theta3*&a1)*(&a1*((D*E-D**2)/(E**2))-&a0*((A*B-B**2)/(B**2))); */
w=&a0*theta3*(&a1-&a0)*(A*B-A**2)/(B**2)+(theta2+theta3*&a1)*(&a1*((D*E-D**2)/(E**2))-&a0*((A*B-A**2)/(B**2)));
/* Corrected from
y=theta3*cmean*(&a1-&a0)*((A*B-B**2)/(B**2))+(theta2+theta3*&a1)*(((D*E-D**2)/(E**2))-((A*B-B**2)/(B**2))); */
/* Also cmean was added to the second term of d3. VV2013 Appendix p14 lacks this, but this */
/* should be there as it is present in the Gamma_tnie d3 in p13. */
y=theta3*cmean*(&a1-&a0)*((A*B-A**2)/(B**2))+(theta2+theta3*&a1)*cmean*(((D*E-D**2)/(E**2))-((A*B-A**2)/(B**2)));
s=(&a1-&a0);
t=t(D/E-A/B);
r=(&a1-&a0)*t(A/B)+&a1*t;
gamma[6,]=x||w||y||zero||s||t||r||z1;
					%if &c^= %then %do;
*/CONDITIONAL CDE*/;
effect[,7]=(theta1)*(&a1-&a0)+(theta3*(&m))*(&a1-&a0);
*/CONDITIONAL NDE*/;
effect[,8]=(theta1+theta3*exp(beta0+beta1*&a0+sum(beta2*t(c)))/(1+exp(beta0+beta1*&a0+sum(beta2*t(c)))))*(&a1-&a0);
*/CONDITIONAL NIE*/;
effect[,9]=(theta2+theta3*&a0)*(
exp(beta0+beta1*&a1+sum(beta2*t(c)))/
(1+exp(beta0+beta1*&a1+sum(beta2*t(c))))-
exp(beta0+beta1*&a0+sum(beta2*t(c)))/
(1+exp(beta0+beta1*&a0+sum(beta2*t(c))))
);
*/CONDITIONAL TNDE*/;
effect[,10]=(theta1+theta3*exp(beta0+beta1*&a1+sum(beta2*t(c)))/(1+exp(beta0+beta1*&a1+sum(beta2*t(c)))))*(&a1-&a0);
*/ CONDITIONAL TNIE*/;
effect[,11]=(theta2+theta3*&a1)*(
exp(beta0+beta1*&a1+sum(beta2*t(c)))/
(1+exp(beta0+beta1*&a1+sum(beta2*t(c))))-
exp(beta0+beta1*&a0+sum(beta2*t(c)))/
(1+exp(beta0+beta1*&a0+sum(beta2*t(c))))
);
*/te conditional*/;
effect[,12]=effect[,8]+effect[,11];
gamma[7,]=z||&m || z1;
B=exp(beta0+beta1*&a0+beta2*t(c));
A=(1+B);
d1=theta3*(A*B-B**2)/A**2;
d2=theta3*&a0*(A*B-B**2)/A**2;
d3=theta3*c*(A*B-B**2)/A**2;
d4=0;
d5=1;
d6=0;
d7=B/A;
d8=z1;
gamma[8,]= d1|| d2 || d3 || d4|| d5|| d6|| t(d7) || d8;
B=exp(beta0+beta1*&a1+beta2*t(c));
A=(1+B);
d1=theta3*(A*B-B**2)/A**2;
d2=theta3*&a1*(A*B-B**2)/A**2;
d3=theta3*c*(A*B-B**2)/A**2;
d4=0;
d5=1;
d6=0;
d7=t(B/A);
d8=z1;
gamma[10,]=d1|| d2 || d3 || d4|| d5|| d6|| d7 || d8;
D=exp(beta0+beta1*&a1+beta2*t(c));
X=(1+D);
B=exp(beta0+beta1*&a0+beta2*t(c));
A=(1+B);
d1=(theta2+theta3*&a1)*((D*X-D**2)/X**2-(A*B-B**2)/A**2);
d2=(theta2+theta3*&a1)*(&a1*(D*X-D**2)/X**2-&a0*(A*B-B**2)/A**2);
d3=c*(theta2+theta3*&a1)*((D*X-D**2)/X**2-(A*B-B**2)/A**2);
d4=0;
d5=0;
d6=t(D/X-B/A);
d7=&a1*d6;
d8=z1;
gamma[11,]=d1|| d2 || d3 || d4|| d5|| d6|| d7 || d8;
d1=(theta2+theta3*&a0)*((D*X-D**2)/X**2-(A*B-B**2)/A**2);
d2=(theta2+theta3*&a0)*(&a1*(D*X-D**2)/X**2-&a0*(A*B-B**2)/A**2);
d3=c*(theta2+theta3*&a0)*((D*X-D**2)/X**2-(A*B-B**2)/A**2);
d4=0;
d5=0;
d6=t(D/X-B/A);
d7=&a0*d6;
d8=z1;
gamma[9,]=d1|| d2 || d3 || d4|| d5|| d6|| d7 || d8;
/* Several corrections were made by @kaz-yos on 2020-03-28 based on VV2013 Appendix p14.
This is mreg logistic, yreg linear, interaction true, cvar non-empty case.
These terms are defined for gamma[12,] for cond te. */
A=exp(beta0+beta1*&a0+beta2*t(c));
B=(1+A);
D=exp(beta0+beta1*&a1+beta2*t(c));
E=(1+D);
/* x: Two (A*B-B**2)'s were corrected to (A*B-A**2). */
/* Changed from
x=theta3*(&a1-&a0)*(A*B-B**2)/(B**2)+(theta2+theta3*&a1)*(((D*E-D**2)/(E**2))-((A*B-B**2)/(B**2))); */
x=theta3*(&a1-&a0)*(A*B-A**2)/(B**2)+(theta2+theta3*&a1)*(((D*E-D**2)/(E**2))-((A*B-A**2)/(B**2)));
/* W: Two (A*B-B**2)'s were corrected to (A*B-A**2). */
/* Changed from
w=&a0*theta3*(&a1-&a0)*(A*B-B**2)/(B**2)+(theta2+theta3*&a1)*(&a1*((D*E-D**2)/(E**2))-&a0*((A*B-B**2)/(B**2)));*/
w=&a0*theta3*(&a1-&a0)*(A*B-A**2)/(B**2)+(theta2+theta3*&a1)*(&a1*((D*E-D**2)/(E**2))-&a0*((A*B-A**2)/(B**2)));
/* y: The second c was added by @kaz-yos on 2020-03-28 based on VV2013 Appendix p14.
Note that VV2013 Appendix p14 omits this, leaving a (vector + scalar) operation.
This c should exist as the second term is the contribution from gamma[11,] for cond tnie.
That is gamma[12,] (Gamma for cond te; p14) is gamma[11,] (Gamma for cond tnie; p13)
+ (a1-a0) * gamma[8,] (Gamma for cond pnde; p12).
Also two (A*B-B**2)'s were corrected to (A*B-A**2). */
/* Changed from
y=theta3*c*(&a1-&a0)*((A*B-B**2)/(B**2))+(theta2+theta3*&a1)*(((D*E-D**2)/(E**2))-((A*B-B**2)/(B**2))); */
y=theta3*c*(&a1-&a0)*((A*B-A**2)/(B**2))+(theta2+theta3*&a1)*c*(((D*E-D**2)/(E**2))-((A*B-A**2)/(B**2)));
s=(&a1-&a0);
t=t(D/E-A/B);
r=(&a1-&a0)*t(A/B)+&a1*t;
gamma[12,]=x||w||y||zero||s||t||r||z1;
					%end;
				%end;
				%if &yreg=linear & &mreg=logistic & &interaction=false %then %do;
*/MARGINAL CDE*/;
effect[,1]=(theta1)*(&a1-&a0);
*/MARGINAL NDE*/;
effect[,2]=(theta1)*(&a1-&a0);
*/MARGINAL NIE*/;
effect[,3]=(theta2)*(exp(beta0+beta1*&a1+sum(beta2*t(cmean)))/
(1+exp(beta0+beta1*&a1+sum(beta2*t(cmean))))-exp(beta0+beta1*&a0+sum(beta2*t(cmean)))/
(1+exp(beta0+beta1*&a0+sum(beta2*t(cmean)))));
*/ MARGINAL TNDE*/;
effect[,4]=(theta1)*(&a1-&a0);
*/ MARGINAL TNIE*/;
effect[,5]=(theta2)*(exp(beta0+beta1*&a1+sum(beta2*t(cmean)))/(1+exp(beta0+beta1*&a1+sum(beta2*t(cmean))))-exp(beta0+beta1*&a0+sum(beta2*t(cmean)))/(1+exp(beta0+beta1*&a0+sum(beta2*t(cmean)))));
*/te marginal*/;
effect[,6]=effect[,2]+effect[,5];
gamma[1,]=z|| z1 ;
A=exp(beta0+beta1*&a0+beta2*t(cmean));
B=(1+A);
x=0;
w=0;
gamma[2,]=  x|| w || z1 || zero|| one|| zero||  z1;
A=exp(beta0+beta1*&a1+beta2*t(cmean));
B=(1+A);
gamma[4,]=x|| w || z1 || zero|| one|| zero||  z1;
D=exp(beta0+beta1*&a1+beta2*t(cmean));
E=(1+D);
A=exp(beta0+beta1*&a0+beta2*t(cmean));
B=(1+A);
x=(theta2)*((D*E-D**2)/E**2-(A*B-A**2)/B**2);
w=(theta2)*(&a1*(D*E-D**2)/E**2-&a0*(A*B-A**2)/B**2);
y=cmean*(theta2)*((D*E-D**2)/E**2-(A*B-A**2)/B**2);
h=t(D/E-A/B);
gamma[5,]=x|| w|| y|| zero||zero|| h||  z1;
gamma[3,]=x|| w|| y||zero||zero|| h||  z1;

/* Fixed Gamma_te by @kaz-yos on 2020-04-01 following VV2013 Appendix p14.
This is the mreg logistic, yreg linear, interaction f, cvar empty case. */
A=exp(beta0+beta1*&a0+beta2*t(Cmean));
B=(1+A);
D=exp(beta0+beta1*&a1+beta2*t(Cmean));
E=(1+D);
/* Corrected from
x=(theta2)*((D*E-E**2)/(E**2)-(A*B-B**2)/(B**2)); */
x=(theta2)*((D*E-D**2)/(E**2)-(A*B-A**2)/(B**2));
/* The squared term in the numerator is different from the denominator. */
/* Corrected from
w=((theta2)*(&a1*(D*E-E**2)/(E**2)-&a0*(A*B-B**2)/(B**2))); */
w=((theta2)*(&a1*(D*E-D**2)/(E**2)-&a0*(A*B-A**2)/(B**2)));
/* Corrected from
y=(theta2)*Cmean*((D*E-E**2)/(E**2)-(A*B-B**2)/(B**2)); */
y=(theta2)*Cmean*((D*E-D**2)/(E**2)-(A*B-A**2)/(B**2));
t=t(D/E-A/B);
s=(&a1-&a0);
gamma[6,]=x||w||y||zero||s||t||z1;

					%if &c^= %then %do;
*/CONDITIONAL CDE*/;
effect[,7]=(theta1)*(&a1-&a0);
*/CONDITIONAL NDE*/;
effect[,8]=(theta1)*(&a1-&a0);
*/CONDITIONAL NIE*/;
effect[,9]=(theta2)*(exp(beta0+beta1*&a1+sum(beta2*t(c)))/(1+exp(beta0+beta1*&a1+sum(beta2*t(c))))-exp(beta0+beta1*&a0+sum(beta2*t(c)))/(1+exp(beta0+beta1*&a0+sum(beta2*t(c)))));
*/CONDITIONAL TNDE*/;
effect[,10]=(theta1)*(&a1-&a0);
*/ CONDITIONAL TNIE*/;
effect[,11]=(theta2)*(exp(beta0+beta1*&a1+sum(beta2*t(c)))/(1+exp(beta0+beta1*&a1+sum(beta2*t(c))))-exp(beta0+beta1*&a0+sum(beta2*t(c)))/(1+exp(beta0+beta1*&a0+sum(beta2*t(c)))));
*/te conditional*/;
effect[,12]=effect[,8]+effect[,11];
gamma[7,]=z||z1;
A=exp(beta0+beta1*&a0+beta2*t(c));
B=(1+A);
x=0;
w=0;
gamma[8,]= x|| w || z1 || zero|| one|| zero|| z1;
A=exp(beta0+beta1*&a1+beta2*t(c));
B=(1+A);
x=0;
w=0;
gamma[10,]=x|| w || z1 || zero|| one|| zero|| z1;
D=exp(beta0+beta1*&a1+beta2*t(c));
E=(1+D);
A=exp(beta0+beta1*&a0+beta2*t(c));
B=(1+A);
x=(theta2)*((D*E-D**2)/E**2-(A*B-A**2)/B**2);
w=(theta2)*(&a1*(D*E-D**2)/E**2-&a0*(A*B-A**2)/B**2);
y=c*(theta2)*((D*E-D**2)/E**2-(A*B-A**2)/B**2);
h=t(D/E-A/B);
gamma[11,]=x|| w|| y|| zero||zero|| h|| z1;
x=(theta2)*((D*E-D**2)/E**2-(A*B-A**2)/B**2);
w=(theta2)*(&a1*(D*E-D**2)/E**2-&a0*(A*B-A**2)/B**2);
y=c*(theta2)*((D*E-D**2)/E**2-(A*B-A**2)/B**2);
h=t(D/E-A/B);
gamma[9,]=x|| w|| y|| zero||zero|| h|| z1;
A=exp(beta0+beta1*&a0+beta2*t(c));
B=(1+A);
D=exp(beta0+beta1*&a1+beta2*t(c));
E=(1+D);
/* Several corrections by @kaz-yos on 2020-03-28 based on VV2013 Appendix p14
This is the mreg logistic, yreg linear, interaction f, cvar non-empty case.
Note ?**2 terms must be different from the one in the denominator. */
/* Corrected from
x=(theta2)*((D*E-E**2)/(E**2)-(A*B-B**2)/(B**2)); */
x=(theta2)*((D*E-D**2)/(E**2)-(A*B-A**2)/(B**2));
/* Corrected from
w=((theta2)*(&a1*(D*E-E**2)/(E**2)-&a0*(A*B-B**2)/(B**2))); */
w=((theta2)*(&a1*(D*E-D**2)/(E**2)-&a0*(A*B-A**2)/(B**2)));
/* Corrected from
y=(theta2)*c*((D*E-E**2)/(E**2)-(A*B-B**2)/(B**2)); */
y=(theta2)*c*((D*E-D**2)/(E**2)-(A*B-A**2)/(B**2));
t=t(D/E-A/B);
s=(&a1-&a0);
gamma[12,]=x||w||y||zero||s||t||z1;

					%end;
				%end;
				%if &yreg^=linear & &mreg=linear & &interaction=true %then %do;
  	  */MARGINAL CDE*/;
	  x6=(theta1+theta3*&m)*(&a1-&a0);
effect[,1]=exp(x6);
	   */MARGINAL NDE*/;
	  x7=(theta1+theta3*beta0+theta3*beta1*&a0+sum(theta3*beta2*t(cmean))+theta3*theta2*rm)*(&a1-&a0)+1/2*tsq*rm*(asq-a1sq);
effect[,2]=exp(x7);
	  */MARGINAL NIE*/;
	  x8=(theta2*beta1+theta3*beta1*&a0)*(&a1-&a0);
effect[,3]=exp(x8);
	  */ MARGINAL TNDE*/;
	  x9=(theta1+theta3*beta0+theta3*beta1*&a1+sum(theta3*beta2*t(cmean))+theta3*theta2*rm)*(&a1-&a0)+1/2*tsq*rm*(asq-a1sq);
effect[,4]=exp(x9);
	  */ MARGINAL TNIE*/;
	  x10=(theta2*beta1+theta3*beta1*&a1)*(&a1-&a0);
effect[,5]=exp(x10);
	  */te*/;
effect[,6]=effect[,2]*effect[,5];
gamma[1,]=z||&m || z1||zero ;
x=theta3*&a0;
w=theta3*t(Cmean);
h=beta0+beta1*&a0+(beta2)*t(Cmean)+theta2*s2+theta3*s2*(&a1+&a0);
ts=s2*theta3;
f=theta3*theta2+0.5*(theta3**2)*(&a1+&a0);
gamma[2,]= theta3|| x|| t(w) || zero|| one||ts|| h || z1||f;
x=theta3*&a1;
w=theta3*t(Cmean);
h=beta0+beta1*&a1+(beta2)*t(Cmean)+theta2*s2+theta3*s2*(&a1+&a0);
ts=s2*theta3;
f=theta3*theta2+0.5*theta3**2*(&a1+&a0);
gamma[4,]=theta3|| x|| t(w)|| zero|| one|| ts|| h||z1||f;
x=theta2+theta3*&a1;
w=beta1*&a1;
gamma[5,]=zero|| x|| z1|| zero||zero|| beta1|| w || z1 || zero;
x=theta2+theta3*&a0;
w=beta1*&a0;
gamma[3,]=zero|| x|| z1|| zero||zero|| beta1|| w || z1 ||zero;
d2pnde=theta3*&a0;
d3pnde=theta3*(Cmean);
d7pnde=beta0+beta1*&a0+(beta2)*t(Cmean)+theta2*s2+theta3*s2*(&a1+&a0);
d6pnde=s2*theta3;
d9pnde=theta3*theta2+0.5*(theta3**2)*(&a1+&a0);
d2tnie=theta2+theta3*&a1;
d7tnie=beta1*&a1;
d2=d2pnde+d2tnie;
d3=d3pnde;
d6=d6pnde+beta1;
d7=d7pnde+d7tnie;
d9=d9pnde;
gamma[6,]=theta3||d2||d3||zero||one||d6||d7||z1||d9;
					%if &c^= %then %do;
*/CONDITIONAL CDE*/;
effect[,7]=exp((theta1+theta3*&m)*(&a1-&a0));
      */CONDITIONAL NDE*/;
effect[,8]=exp((theta1+theta3*beta0+theta3*beta1*&a0+sum(theta3*beta2*t(c))+theta3*theta2*rm)*(&a1-&a0)+(1/2)*tsq*rm*(asq-a1sq)
);
      */CONDITIONAL NIE*/;
	  x3=(theta2*beta1+theta3*beta1*&a0)*(&a1-&a0);
effect[,9]=exp(x3);
      */CONDITIONAL TNDE*/;
	  x4=(theta1+theta3*beta0+theta3*beta1*&a1+sum(theta3*beta2*t(c))+theta3*theta2*rm)*(&a1-&a0)+(1/2)*tsq*rm*(asq-a1sq);
effect[,10]=exp(x4);
      */ CONDITIONAL TNIE*/;
	  x5=(theta2*beta1+theta3*beta1*&a1)*(&a1-&a0);
effect[,11]=exp(x5);
	  */te*/;
effect[,12]=effect[,8]*effect[,11];
gamma[7,]=z||&m || z1||zero;
x=theta3*&a0;
w=theta3*t(c);
h=beta0+beta1*&a0+(beta2)*t(c)+theta2*s2+theta3*s2*(&a1+&a0);
ts=s2*theta3;
f=theta3*theta2+0.5*(theta3**2)*(&a1+&a0);
gamma[8,]= theta3|| x|| t(w) || zero|| one||ts|| h ||z1||f;
x=theta3*&a1;
w=theta3*t(c);
h=beta0+beta1*&a1+(beta2)*t(c)+theta2*s2+theta3*s2*(&a1+&a0);
ts=s2*theta3;
f=theta3*theta2+0.5*theta3**2*(&a1+&a0);
gamma[10,]=theta3|| x|| t(w)|| zero|| one|| ts|| h||z1||f;
x=theta2+theta3*&a1;
w=beta1*&a1;
gamma[11,]=zero|| x|| z1|| zero||zero|| beta1|| w || z1 || zero;
x=theta2+theta3*&a0;
w=beta1*&a0;
gamma[9,]=zero|| x|| z1|| zero||zero|| beta1|| w || z1 ||zero;
d2pnde=theta3*&a0;
d3pnde=theta3*(c);
d7pnde=beta0+beta1*&a0+(beta2)*t(c)+theta2*s2+theta3*s2*(&a1+&a0);
d6pnde=s2*theta3;
d9pnde=theta3*theta2+0.5*(theta3**2)*(&a1+&a0);
d2tnie=theta2+theta3*&a1;
d7tnie=beta1*&a1;
d2=d2pnde+d2tnie;
d3=d3pnde;
d6=d6pnde+beta1;
d7=d7pnde+d7tnie;
d9=d9pnde;
gamma[12,]=theta3||d2||d3||zero||one||d6||d7||z1||d9;

					%end;
				%end;
				%if &yreg^=linear & &mreg=logistic & &interaction=false %then %do;
*/MARGINAL CDE*/;
	  x6=(theta1)*(&a1-&a0);
effect[,1]=exp(x6);
	   */MARGINAL NDE*/;
effect[,2]=exp((theta1)*(&a1-&a0))*
(1+exp(theta2+beta0+beta1*&a0+sum(beta2*t(cmean))))/
(1+exp(theta2+beta0+beta1*&a0+sum(beta2*t(cmean))));
	  */MARGINAL NIE*/;
effect[,3]=(
(1+exp(beta0+beta1*&a0+sum(beta2*t(cmean))))*
(1+exp(theta2+beta0+beta1*&a1+sum(beta2*t(cmean))))
)/
(
(1+exp(beta0+beta1*&a1+sum(beta2*t(cmean))))*
(1+exp(theta2+beta0+beta1*&a0+sum(beta2*t(cmean))))
);
	  */ MARGINAL TNDE*/;
effect[,4]=exp((theta1)*(&a1-&a0))*
(1+exp(theta2+beta0+beta1*&a1+sum(beta2*t(cmean))))/
(1+exp(theta2+beta0+beta1*&a1+sum(beta2*t(cmean))));
	  */ MARGINAL TNIE*/;
effect[,5]=(
(1+exp(beta0+beta1*&a0+sum(beta2*t(cmean))))*
(1+exp(theta2+beta0+beta1*&a1+sum(beta2*t(cmean))))
)/
(
(1+exp(beta0+beta1*&a1+sum(beta2*t(cmean))))*
(1+exp(theta2+beta0+beta1*&a0+sum(beta2*t(cmean))))
);
effect[,6]=effect[,2]*effect[,5];
/* Modified by @kaz-yos on 2020-04-01 based on V2015 p474 on Gamma_cde. */
/* This is the mreg logistic, yreg non-linear, int f case.*/
/* z is defined as z=zero||zero||z1||zero||one||zero; above without common factor (a1-a0).*/
/* In the case of &yreg^=linear & &mreg=logistic no common factor multiplication is done. */
/* Confirm looking for "%if (&mreg=logistic & &yreg^=linear) %then %do;" */
/* So it must be included in Gamma_cde. */
/* Changed from
gamma[1,]=z|| z1 ; */
gamma[1,]=(z|| z1) * (&a1-&a0);
A=exp(theta2+beta0+beta1*&a0+beta2*t(Cmean));
B=(1+exp(theta2+beta0+beta1*&a0+beta2*t(Cmean)));
D=exp(theta2+beta0+beta1*&a0+beta2*t(Cmean));
E=(1+exp(theta2+beta0+beta1*&a0+beta2*t(Cmean)));
d1nde=A/B-D/E;
d2nde=&a0*(d1nde);
d3nde=(d1nde)*(Cmean);
d4nde=0;
d5nde=(&a1-&a0);
d6nde=d1nde;
d7nde=z1;
gamma[2,]= d1nde|| d2nde|| d3nde || d4nde||d5nde||d6nde||d7nde;
A=exp(theta2+beta0+beta1*&a1+beta2*t(Cmean));
B=(1+exp(theta2+beta0+beta1*&a1+beta2*t(Cmean)));
D=exp(theta2+beta0+beta1*&a1+beta2*t(Cmean));
E=(1+exp(theta2+beta0+beta1*&a1+beta2*t(Cmean)));
s=A/B-D/E;
x=&a1*(s);
w=(s)*(Cmean);
t=(&a1-&a0);
gamma[4,]=s|| x|| w || zero||t||s||z1;
A=exp(theta2+beta0+beta1*&a1+beta2*t(Cmean));
B=(1+exp(theta2+beta0+beta1*&a1+beta2*t(Cmean)));
D=exp(theta2+beta0+beta1*&a0+beta2*t(Cmean));
E=(1+exp(theta2+beta0+beta1*&a0+beta2*t(Cmean)));
F=exp(beta0+beta1*&a0+beta2*t(Cmean));
G=(1+exp(beta0+beta1*&a0+beta2*t(Cmean)));
H=exp(beta0+beta1*&a1+beta2*t(Cmean));
I=(1+exp(beta0+beta1*&a1+beta2*t(Cmean)));
d1nie=F/G-H/I+A/B-D/E;
d2nie=&a0*F/G-&a1*H/I+&a1*A/B-&a0*D/E;
d3nie=Cmean*(d1nie);
d4nie=0;
d5nie=0;
d6nie=(A/B-D/E);
d7nie=z1;
gamma[5,]=d1nie|| d2nie|| d3nie || d4nie||d5nie||d6nie||d7nie;
A=exp(theta2+beta0+beta1*&a1+beta2*t(Cmean));
B=(1+exp(theta2+beta0+beta1*&a1+beta2*t(Cmean)));
D=exp(theta2+beta0+beta1*&a0+beta2*t(Cmean));
E=(1+exp(theta2+beta0+beta1*&a0+beta2*t(Cmean)));
F=exp(beta0+beta1*&a0+beta2*t(Cmean));
G=(1+exp(beta0+beta1*&a0+beta2*t(Cmean)));
H=exp(beta0+beta1*&a1+beta2*t(Cmean));
I=(1+exp(beta0+beta1*&a1+beta2*t(Cmean)));
/* Added by @kaz-yos on 2020-04-01 based on V2015 p474. */
/* Q' = A/B; B' = D/E; K' = H/I; D' = F/G */
/* (D'+Q') - (K'+B') = (F/G + A/B) - (H/I + D/E)  */
/* Changed from
s=F/G-H/I+A/B-D/E; (cosmetic change only) */
s=(F/G + A/B) - (H/I + D/E);
/* a0(D'-B') + a(Q'-K') = &a0 * (F/G - D/E) + &a1 * (A/B - H/I) */
/* This change needs verification. VV2013 and V2015 only cover Gamma_tnie. */
/* Reasoning by @kaz-yos on 2020-04-02. */
/* Based on Pearl's decomposition (V2015 p465), the treatment indexing the outcome model
 * sum_m E[Y|a1 <- this one ,m,c] {P(m|a1,c) - P(m|a0,c)} is the "a1" that makes it
 * _total_ NIE. Thus, this a1 in the outcome model is associated with theta1 and theta3.
 * Thus, in the OR^TNIE expression in V2015 p473, a1 associated with theta3 are the ones
 * that need to change to a0 to give the OR^PNIE expression. In the Gamma_tnie to Gamma_pnie
 * change, the a1 -> a0 change should thus only occur for the terms involving theta3.
 * In the d2 expression (V2015 p474), These are only in Q and B terms. a0 and a1 in front of
 * [D-B] and [Q-K], respectively, come from the terms involving beta1 (d2 is a partial wrt beta1).
 * Thus, these should not be changed and remain the same as the corresponding term in Gamma_tnie.
 */
/* Changed from
x=&a1*F/G-&a0*H/I+&a0*A/B-&a1*D/E; */
x=&a0 * (F/G - D/E) + &a1 * (A/B - H/I);
w=Cmean*(s);
k=(A/B-D/E);
gamma[3,]=s|| x|| w || zero||zero||k||z1;
d1=(d1nie+d1nde);
d2=(d2nie+d2nde);
d3=(d3nie+d3nde);
d4=zero;
d5=(d5nie+d5nde);
d6=(d6nie+d6nde);
gamma[6,]=d1||d2||d3||zero||d5||d6||z1;

					%if &c^= %then %do;
*/CONDITIONAL CDE*/;
	  x1=exp(theta1*(&a1-&a0));
effect[,7]=x1;
      */CONDITIONAL NDE*/;
effect[,8]=exp((theta1)*(&a1-&a0))*
(1+exp(theta2+beta0+beta1*&a0+sum(beta2*t(c))))/
(1+exp(theta2+beta0+beta1*&a0+sum(beta2*t(c))));
      */CONDITIONAL NIE*/;
effect[,9]=(
(1+exp(beta0+beta1*&a0+beta2*t(c)))*
(1+exp(theta2+beta0+beta1*&a1+sum(beta2*t(c))))
)/
(
(1+exp(beta0+beta1*&a1+sum(beta2*t(c))))*
(1+exp(theta2+beta0+beta1*&a0+sum(beta2*t(c))))
);
      */CONDITIONAL TNDE*/;
effect[,10]=exp((theta1)*(&a1-&a0))*
(1+exp(theta2+beta0+beta1*&a1+sum(beta2*t(c))))/
(1+exp(theta2+beta0+beta1*&a1+sum(beta2*t(c))));
      */ CONDITIONAL TNIE*/;
effect[,11]=(
(1+exp(beta0+beta1*&a0+sum(beta2*t(c))))*
(1+exp(theta2+beta0+beta1*&a1+sum(beta2*t(c))))
)/
(
(1+exp(beta0+beta1*&a1+sum(beta2*t(c))))*
(1+exp(theta2+beta0+beta1*&a0+sum(beta2*t(c))))
);
effect[,12]=effect[,8]*effect[,11];
/* Modified by @kaz-yos on 2020-04-01 based on V2015 p474 on Gamma_cde. */
/* This is the mreg logistic, yreg non-linear, int f, cvar non-empty case.*/
/* z is defined as z=zero||zero||z1||zero||one||zero; above without common factor (a1-a0).*/
/* In the case of &yreg^=linear & &mreg=logistic no common factor multiplication is done. */
/* Confirm looking for "%if (&mreg=logistic & &yreg^=linear) %then %do;" */
/* So it must be included in Gamma_cde. */
/* Changed from
gamma[7,]=z||z1; */
gamma[7,]=(z||z1) * (&a1-&a0);
A=exp(theta2+beta0+beta1*&a0+beta2*t(c));
B=(1+exp(theta2+beta0+beta1*&a0+beta2*t(c)));
D=exp(theta2+beta0+beta1*&a0+beta2*t(c));
E=(1+exp(theta2+beta0+beta1*&a0+beta2*t(c)));
d1cnde=A/B-D/E;
d2cnde=&a0*(d1cnde);
d3cnde=(d1cnde)*(c);
d4cnde=0;
d5cnde=(&a1-&a0);
d6cnde=d1cnde;
d7cnde=z1;
gamma[8,]= d1cnde|| d2cnde||d3cnde|| d4cnde||d5cnde||d6cnde||d7cnde;
A=exp(theta2+beta0+beta1*&a1+beta2*t(c));
B=(1+exp(theta2+beta0+beta1*&a1+beta2*t(c)));
D=exp(theta2+beta0+beta1*&a1+beta2*t(c));
E=(1+exp(theta2+beta0+beta1*&a1+beta2*t(c)));
s=A/B-D/E;
x=&a1*(s);
w=(s)*(c);
t=(&a1-&a0);
gamma[10,]=s|| x|| w || zero||t||s|| z1;
A=exp(theta2+beta0+beta1*&a1+beta2*t(c));
B=(1+exp(theta2+beta0+beta1*&a1+beta2*t(c)));
D=exp(theta2+beta0+beta1*&a0+beta2*t(c));
E=(1+exp(theta2+beta0+beta1*&a0+beta2*t(c)));
F=exp(beta0+beta1*&a0+beta2*t(c));
G=(1+exp(beta0+beta1*&a0+beta2*t(c)));
H=exp(beta0+beta1*&a1+beta2*t(c));
I=(1+exp(beta0+beta1*&a1+beta2*t(c)));
d1cnie=F/G-H/I+A/B-D/E;
d2cnie=&a0*F/G-&a1*H/I+&a1*A/B-&a0*D/E;
d3cnie=c*(d1cnie);
d4cnie=0;
d5cnie=0;
d6cnie=(A/B-D/E);
d7cnie=z1;
gamma[11,]=d1cnie|| d2cnie||d3cnie || d4cnie||d5cnie||d6cnie||d7cnie;
A=exp(theta2+beta0+beta1*&a1+beta2*t(c));
B=(1+exp(theta2+beta0+beta1*&a1+beta2*t(c)));
D=exp(theta2+beta0+beta1*&a0+beta2*t(c));
E=(1+exp(theta2+beta0+beta1*&a0+beta2*t(c)));
F=exp(beta0+beta1*&a0+beta2*t(c));
G=(1+exp(beta0+beta1*&a0+beta2*t(c)));
H=exp(beta0+beta1*&a1+beta2*t(c));
I=(1+exp(beta0+beta1*&a1+beta2*t(c)));
/* Added by @kaz-yos on 2020-04-01 based on V2015 p474. */
/* Q' = A/B; B' = D/E; K' = H/I; D' = F/G */
/* (D'+Q') - (K'+B') = (F/G + A/B) - (H/I + D/E)  */
/* Changed from
s=F/G-H/I+A/B-D/E; (cosmetic change only) */
s=(F/G + A/B) - (H/I + D/E);
/* a0(D'-B') + a(Q'-K') = &a0 * (F/G - D/E) + &a1 * (A/B - H/I) */
/* This change needs verification. VV2013 and V2015 only cover Gamma_tnie. */
/* Reasoning by @kaz-yos on 2020-04-02. */
/* Based on Pearl's decomposition (V2015 p465), the treatment indexing the outcome model
 * sum_m E[Y|a1 <- this one ,m,c] {P(m|a1,c) - P(m|a0,c)} is the "a1" that makes it
 * _total_ NIE. Thus, this a1 in the outcome model is associated with theta1 and theta3.
 * Thus, in the OR^TNIE expression in V2015 p473, a1 associated with theta3 are the ones
 * that need to change to a0 to give the OR^PNIE expression. In the Gamma_tnie to Gamma_pnie
 * change, the a1 -> a0 change should thus only occur for the terms involving theta3.
 * In the d2 expression (V2015 p474), These are only in Q and B terms. a0 and a1 in front of
 * [D-B] and [Q-K], respectively, come from the terms involving beta1 (d2 is a partial wrt beta1).
 * Thus, these should not be changed and remain the same as the corresponding term in Gamma_tnie.
 */
/* Changed from
x=&a1*F/G-&a0*H/I+&a0*A/B-&a1*D/E; */
x=&a0 * (F/G - D/E) + &a1 * (A/B - H/I);
w=c*(s);
k=(A/B-D/E);
gamma[9,]=s|| x|| w || zero||zero||k||z1;
d1=(d1cnie+d1cnde);
d2=(d2cnie+d2cnde);
d3=(d3cnie+d3cnde);
d4=zero;
d5=(d5cnie+d5cnde);
d6=(d6cnie+d6cnde);
gamma[12,]=d1||d2||d3||zero||d5||d6||z1;
					%end;
				%end;
				%if &yreg^=linear & &mreg=logistic & &interaction=true %then %do;
*/MARGINAL CDE*/;
x6=(theta1+theta3*&m)*(&a1-&a0);
effect[,1]=exp(x6);
*/MARGINAL NDE*/;
effect[,2]=exp(theta1*(&a1-&a0))*
(1+exp(theta2+theta3*&a1+beta0+beta1*&a0+sum(beta2*t(cmean))))/
(1+exp(theta2+theta3*&a0+beta0+beta1*&a0+sum(beta2*t(cmean))));
*/MARGINAL NIE*/;
effect[,3]=(
(1+exp(beta0+beta1*&a0+sum(beta2*t(cmean))))*
(1+exp(theta2+theta3*&a0+beta0+beta1*&a1+sum(beta2*t(cmean))))
)/
(
(1+exp(beta0+beta1*&a1+sum(beta2*t(cmean))))*
(1+exp(theta2+theta3*&a0+beta0+beta1*&a0+sum(beta2*t(cmean))))
);
*/ MARGINAL TNDE*/;
effect[,4]=exp(theta1*(&a1-&a0))*
(1+exp(theta2+theta3*&a1+beta0+beta1*&a1+sum(beta2*t(cmean))))/
(1+exp(theta2+theta3*&a0+beta0+beta1*&a1+sum(beta2*t(cmean))));
*/ MARGINAL TNIE*/;
effect[,5]=(
(1+exp(beta0+beta1*&a0+sum(beta2*t(cmean))))*
(1+exp(theta2+theta3*&a1+beta0+beta1*&a1+sum(beta2*t(cmean))))
)/
(
(1+exp(beta0+beta1*&a1+sum(beta2*t(cmean))))*
(1+exp(theta2+theta3*&a1+beta0+beta1*&a0+sum(beta2*t(cmean))))
);
effect[,6]=effect[,2]*effect[,5];
/* Modified by @kaz-yos on 2020-04-01 based on V2015 p474 on Gamma_cde. */
/* This is the mreg logistic, yreg non-linear, int t case.*/
/* z is defined as z=zero||zero||z1||zero||one||zero; above without common factor (a1-a0).*/
/* In the case of &yreg^=linear & &mreg=logistic no common factor multiplication is done. */
/* Confirm looking for "%if (&mreg=logistic & &yreg^=linear) %then %do;" */
/* So it must be included in Gamma_cde. */
/* Changed from
gamma[1,]=z||&m || z1 ; */
gamma[1,]=(z||&m || z1) * (&a1-&a0);
A=exp(theta2+theta3*&a1+beta0+beta1*&a0+beta2*t(cmean));
B=(1+exp(theta2+theta3*&a1+beta0+beta1*&a0+beta2*t(cmean)));
D=exp(theta2+theta3*&a0+beta0+beta1*&a0+beta2*t(cmean));
E=(1+exp(theta2+theta3*&a0+beta0+beta1*&a0+beta2*t(cmean)));
d1nde=A/B-D/E;
d2nde=&a0*(d1nde);
d3nde=(d1nde)*(Cmean);
d4nde=0;
d5nde=(&a1-&a0);
d6nde=d1nde;
d7nde=&a1*A/B-&a0*D/E;
d8nde=z1;
gamma[2,]= d1nde|| d2nde|| d3nde || d4nde||d5nde||d6nde|| d7nde ||d8nde;
A=exp(theta2+theta3*&a1+beta0+beta1*&a1+beta2*t(cmean));
B=(1+exp(theta2+theta3*&a1+beta0+beta1*&a1+beta2*t(Cmean)));
D=exp(theta2+theta3*&a0+beta0+beta1*&a1+beta2*t(Cmean));
E=(1+exp(theta2+theta3*&a0+beta0+beta1*&a1+beta2*t(Cmean)));
s=A/B-D/E;
x=&a1*(s);
w=(s)*(Cmean);
t=(&a1-&a0);
h=&a1*A/B-&a0*D/E;
gamma[4,]=s|| x|| w || zero||t||s|| h ||z1;
A=exp(theta2+theta3*&a1+beta0+beta1*&a1+beta2*t(Cmean));
B=(1+exp(theta2+theta3*&a1+beta0+beta1*&a1+beta2*t(Cmean)));
D=exp(theta2+theta3*&a1+beta0+beta1*&a0+beta2*t(Cmean));
E=(1+exp(theta2+theta3*&a1+beta0+beta1*&a0+beta2*t(Cmean)));
F=exp(beta0+beta1*&a0+beta2*t(Cmean));
G=(1+exp(beta0+beta1*&a0+beta2*t(Cmean)));
H=exp(beta0+beta1*&a1+beta2*t(Cmean));
I=(1+exp(beta0+beta1*&a1+beta2*t(Cmean)));
d1nie=F/G-H/I+A/B-D/E;
d2nie=&a0*F/G-&a1*H/I+&a1*A/B-&a0*D/E;
d3nie=Cmean*(d1nie);
d4nie=0;
d5nie=0;
d6nie=A/B-D/E;
d7nie=&a1*(A/B-D/E);
d8nie=z1;
gamma[5,]=d1nie|| d2nie|| d3nie || d4nie||d5nie||d6nie|| d7nie ||d8nie;

A=exp(theta2+theta3*&a0+beta0+beta1*&a1+beta2*t(Cmean));
B=(1+exp(theta2+theta3*&a0+beta0+beta1*&a1+beta2*t(Cmean)));
D=exp(theta2+theta3*&a0+beta0+beta1*&a0+beta2*t(Cmean));
E=(1+exp(theta2+theta3*&a0+beta0+beta1*&a0+beta2*t(Cmean)));
F=exp(beta0+beta1*&a0+beta2*t(Cmean));
G=(1+exp(beta0+beta1*&a0+beta2*t(Cmean)));
H=exp(beta0+beta1*&a1+beta2*t(Cmean));
I=(1+exp(beta0+beta1*&a1+beta2*t(Cmean)));
/* Added by @kaz-yos on 2020-04-01 based on V2015 p474. */
/* Q' = A/B; B' = D/E; K' = H/I; D' = F/G */
/* (D'+Q') - (K'+B') = (F/G + A/B) - (H/I + D/E)  */
/* Changed from
s=F/G-H/I+A/B-D/E; (cosmetic change only) */
s=F/G-H/I+A/B-D/E;
/* a0(D'-B') + a(Q'-K') = &a0 * (F/G - D/E) + &a1 * (A/B - H/I) */
/* This change needs verification. VV2013 and V2015 only cover Gamma_tnie. */
/* Reasoning by @kaz-yos on 2020-04-02. */
/* Based on Pearl's decomposition (V2015 p465), the treatment indexing the outcome model
 * sum_m E[Y|a1 <- this one ,m,c] {P(m|a1,c) - P(m|a0,c)} is the "a1" that makes it
 * _total_ NIE. Thus, this a1 in the outcome model is associated with theta1 and theta3.
 * Thus, in the OR^TNIE expression in V2015 p473, a1 associated with theta3 are the ones
 * that need to change to a0 to give the OR^PNIE expression. In the Gamma_tnie to Gamma_pnie
 * change, the a1 -> a0 change should thus only occur for the terms involving theta3.
 * In the d2 expression (V2015 p474), These are only in Q and B terms. a0 and a1 in front of
 * [D-B] and [Q-K], respectively, come from the terms involving beta1 (d2 is a partial wrt beta1).
 * Thus, these should not be changed and remain the same as the corresponding term in Gamma_tnie.
 */
/* Changed from
x=&a1*F/G-&a0*H/I+&a0*A/B-&a1*D/E; */
x=&a0 * (F/G - D/E) + &a1 * (A/B - H/I);
w=Cmean*(s);
l=A/B-D/E;
k=&a0*(A/B-D/E);
gamma[3,]=s|| x|| w || zero||zero||l|| k ||z1;

d1=((d1nie)+(d1nde));
d2=((d2nie)+(d2nde));
d3=((d3nie)+(d3nde));
d4=((d4nie)+(d4nde));
d5=((d5nie)+(d5nde));
d6=((d6nie)+(d6nde));
d7=((d7nie)+(d7nde));
gamma[6,]=d1||d2||d3||d4||d5||d6||d7||z1;
					%if &c^= %then %do;
*/CONDITIONAL CDE*/;
x1=exp((theta1+theta3*&m)*(&a1-&a0));
effect[,7]=x1;
*/CONDITIONAL NDE*/;
effect[,8]=exp(theta1*(&a1-&a0))*
(1+exp(theta2+theta3*&a1+beta0+beta1*&a0+sum(beta2*t(c))))/
(1+exp(theta2+theta3*&a0+beta0+beta1*&a0+sum(beta2*t(c))));
*/CONDITIONAL NIE*/;
effect[,9]=(
(1+exp(beta0+beta1*&a0+sum(beta2*t(c))))*
(1+exp(theta2+theta3*&a0+beta0+beta1*&a1+sum(beta2*t(c))))
)/
(
(1+exp(beta0+beta1*&a1+sum(beta2*t(c))))*
(1+exp(theta2+theta3*&a0+beta0+beta1*&a0+sum(beta2*t(c))))
);
*/CONDITIONAL TNDE*/;
effect[,10]=exp(theta1*(&a1-&a0))*
(1+exp(theta2+theta3*&a1+beta0+beta1*&a1+sum(beta2*t(c))))/
(1+exp(theta2+theta3*&a0+beta0+beta1*&a1+sum(beta2*t(c))));
*/ CONDITIONAL TNIE*/;
effect[,11]=(
(1+exp(beta0+beta1*&a0+sum(beta2*t(c))))*
(1+exp(theta2+theta3*&a1+beta0+beta1*&a1+sum(beta2*t(c))))
)/
(
(1+exp(beta0+beta1*&a1+sum(beta2*t(c))))*
(1+exp(theta2+theta3*&a1+beta0+beta1*&a0+sum(beta2*t(c))))
);
effect[,12]=effect[,8]*effect[,11];
/* Modified by @kaz-yos on 2020-04-01 based on V2015 p474 on Gamma_cde. */
/* This is the mreg logistic, yreg non-linear, int t, cvar non-empty case.*/
/* z is defined as z=zero||zero||z1||zero||one||zero; above without common factor (a1-a0).*/
/* In the case of &yreg^=linear & &mreg=logistic no common factor multiplication is done. */
/* Confirm looking for "%if (&mreg=logistic & &yreg^=linear) %then %do;" */
/* So it must be included in Gamma_cde. */
/* Changed from
gamma[7,]=z||&m || z1; */
gamma[7,]=(z||&m || z1) * (&a1-&a0);
A=exp(theta2+theta3*&a1+beta0+beta1*&a0+beta2*t(c));
B=(1+exp(theta2+theta3*&a1+beta0+beta1*&a0+beta2*t(c)));
D=exp(theta2+theta3*&a0+beta0+beta1*&a0+beta2*t(c));
E=(1+exp(theta2+theta3*&a0+beta0+beta1*&a0+beta2*t(c)));
d1cnde=A/B-D/E;
d2cnde=&a0*(d1cnde);
d3cnde=(d1cnde)*(c);
d4cnde=0;
d5cnde=(&a1-&a0);
d6cnde=A/B-D/E;
d7cnde=&a1*A/B-&a0*D/E;
d8cnde=z1;
gamma[8,]= d1cnde|| d2cnde|| d3cnde || d4cnde||d5cnde||d6cnde|| d7cnde ||d8cnde;
A=exp(theta2+theta3*&a1+beta0+beta1*&a1+beta2*t(c));
B=(1+exp(theta2+theta3*&a1+beta0+beta1*&a1+beta2*t(c)));
D=exp(theta2+theta3*&a0+beta0+beta1*&a1+beta2*t(c));
E=(1+exp(theta2+theta3*&a0+beta0+beta1*&a1+beta2*t(c)));
d1=A/B-D/E;
d2=&a1*(d1);
d3=(d1)*(c);
d4=0;
d5=(&a1-&a0);
d6=A/B-D/E;
d7=&a1*A/B-&a0*D/E;
d8=z1;
gamma[10,]=d1|| d2|| d3 || d4||d5||d6|| d7 ||d8;
A=exp(theta2+theta3*&a1+beta0+beta1*&a1+beta2*t(c));
B=(1+exp(theta2+theta3*&a1+beta0+beta1*&a1+beta2*t(c)));
D=exp(theta2+theta3*&a1+beta0+beta1*&a0+beta2*t(c));
E=(1+exp(theta2+theta3*&a1+beta0+beta1*&a0+beta2*t(c)));
F=exp(beta0+beta1*&a0+beta2*t(c));
G=(1+exp(beta0+beta1*&a0+beta2*t(c)));
H=exp(beta0+beta1*&a1+beta2*t(c));
I=(1+exp(beta0+beta1*&a1+beta2*t(c)));
d1cnie=F/G-H/I+A/B-D/E;
d2cnie=&a0*F/G-&a1*H/I+&a1*A/B-&a0*D/E;
d3cnie=c*(d1cnie);
d4cnie=0;
d5cnie=0;
d6cnie=A/B-D/E;
d7cnie=&a1*(A/B-D/E);
d8cnie=z1;
gamma[11,]=d1cnie|| d2cnie|| d3cnie || d4cnie||d5cnie||d6cnie|| d7cnie ||d8cnie;
A=exp(theta2+theta3*&a0+beta0+beta1*&a1+beta2*t(c));
B=(1+exp(theta2+theta3*&a0+beta0+beta1*&a1+beta2*t(c)));
D=exp(theta2+theta3*&a0+beta0+beta1*&a0+beta2*t(c));
E=(1+exp(theta2+theta3*&a0+beta0+beta1*&a0+beta2*t(c)));
F=exp(beta0+beta1*&a0+beta2*t(c));
G=(1+exp(beta0+beta1*&a0+beta2*t(c)));
H=exp(beta0+beta1*&a1+beta2*t(c));
I=(1+exp(beta0+beta1*&a1+beta2*t(c)));
/* Added by @kaz-yos on 2020-04-01 based on V2015 p474. */
/* Q' = A/B; B' = D/E; K' = H/I; D' = F/G */
/* (D'+Q') - (K'+B') = (F/G + A/B) - (H/I + D/E)  */
/* Changed from
d1=F/G-H/I+A/B-D/E; (cosmetic change only) */
d1=(F/G + A/B) - (H/I + D/E);
/* a0(D'-B') + a(Q'-K') = &a0 * (F/G - D/E) + &a1 * (A/B - H/I) */
/* This change needs verification. VV2013 and V2015 only cover Gamma_tnie. */
/* Reasoning by @kaz-yos on 2020-04-02. */
/* Based on Pearl's decomposition (V2015 p465), the treatment indexing the outcome model
 * sum_m E[Y|a1 <- this one ,m,c] {P(m|a1,c) - P(m|a0,c)} is the "a1" that makes it
 * _total_ NIE. Thus, this a1 in the outcome model is associated with theta1 and theta3.
 * Thus, in the OR^TNIE expression in V2015 p473, a1 associated with theta3 are the ones
 * that need to change to a0 to give the OR^PNIE expression. In the Gamma_tnie to Gamma_pnie
 * change, the a1 -> a0 change should thus only occur for the terms involving theta3.
 * In the d2 expression (V2015 p474), These are only in Q and B terms. a0 and a1 in front of
 * [D-B] and [Q-K], respectively, come from the terms involving beta1 (d2 is a partial wrt beta1).
 * Thus, these should not be changed and remain the same as the corresponding term in Gamma_tnie.
 */
/* Changed from
d2=&a1*F/G-&a0*H/I+&a0*A/B-&a1*D/E; */
d2=&a0 * (F/G - D/E) + &a1 * (A/B - H/I);
d3=c*(d1);
d4=0;
d5=0;
d6=A/B-D/E;
d7=&a0*(A/B-D/E);
d8=z1;
gamma[9,]=d1|| d2|| d3 || d4||d5||d6|| d7 ||d8;
d1=((d1cnie)+(d1cnde));
d2=((d2cnie)+(d2cnde));
d3=((d3cnie)+(d3cnde));
d4=((d4cnie)+(d4cnde));
d5=((d5cnie)+(d5cnde));
d6=((d6cnie)+(d6cnde));
d7=((d7cnie)+(d7cnde));
gamma[12,]=d1||d2||d3||d4||d5||d6||d7||z1;
					%end;
				%end;
			%end;

			%if &cvar= %then %do;
				%if  &interaction=false %then %do;
effect=J(1,6);
gamma=J(6,5);
				%end;
				%if  &interaction=true %then %do;
effect=J(1,6);
gamma=J(6,6);
				%end;
USE out1;
READ ALL INTO VB;
NVB1= NROW(VB);
V1=VB[2:NVB1,];
theta1=VB[1,2];
theta2=VB[1,3];
				%if &interaction=true %then %do;
theta3=VB[1,4] ;
				%end;
USE out2;
READ ALL INTO VB;
NVB2= NROW(VB);
V2=VB[2:NVB2,];
				%if (&yreg=linear & &mreg=linear) | (&mreg=logistic) %then %do;
beta0=VB[1,1];
beta1=VB[1,2];
zero1=J(nrow(V1),nrow(V2),0);
zero2=J(nrow(V2),nrow(V1),0);
A= V2 || zero2;
B= zero1 || V1;
sigma= A // B;
zero=0;
one=1;
z=zero||zero||zero||one||zero;
				%end;
				%if &yreg^=linear & &mreg=linear %then %do;
s2=VB[1,1];
s2=s2**2;
beta0=VB[1,2];
beta1=VB[1,3];
tsq=(theta3**2);
rm=s2;
asq=(&a1**2);
a1sq=(&a0**2);
NVB2= NROW(VB);
colvb=ncol(vb);
V2=VB[2:NVB2,2:colvb];
zero1=J(nrow(V1),nrow(V2),0);
zero2=J(nrow(V2),nrow(V1),0);
z2=J(nrow(V1),1,0);
z3=J(nrow(V2),1,0);
A= V2 || zero2 ||z3;
B= zero1 || V1||z2;
zeros=J(1,nrow(V1)+nrow(V2),0);
/* @kaz-yos on 2020-05-02 */
/* This is wrong. It should be the following based on V2015 p470.
where n = sample size, p = length(betas), s2 = sigma^2
D= zeros || ((2 * (s2**2)) / (n - p)) */
D= zeros ||s2;
sigma= A // B//D;
zero=0;
one=1;
z=zero||zero||zero||one||zero;
%if  &interaction=false %then %do;

gamma=J(6,6);
				%end;
				%if  &interaction=true %then %do;

gamma=J(6,7);
				%end;
				%end;
				%if &yreg=linear & &mreg=linear & &interaction=true %then %do;
*/CONDITIONAL=MARGINAL CDE*/;
effect[,1]=(theta1)*(&a1-&a0)+(theta3*(&m))*(&a1-&a0);
*/CONDITIONAL=MARGINAL NDE*/;
effect[,2]=(theta1+theta3*beta0+theta3*beta1*&a0)*(&a1-&a0);
*/CONDITIONAL=MARGINAL NIE*/;
effect[,3]=(theta2*beta1+theta3*beta1*&a0)*(&a1-&a0);
*/CONDITIONAL=MARGINAL TNDE*/;
effect[,4]=(theta1+theta3*beta0+theta3*beta1*&a1)*(&a1-&a0);
*/ CONDITIONAL=MARGINAL TNIE*/;
effect[,5]=(theta2*beta1+theta3*beta1*&a1)*(&a1-&a0);
*/te*/;
effect[,6]=(theta1+theta3*beta0+theta3*beta1*&a0+theta2*beta1+theta3*beta1*&a1)*(&a1-&a0);
gamma[1,]=z||&m;
x1=theta3*&a0;
h1=beta0+beta1*&a0;
gamma[2,]= theta3|| x1||  zero|| one|| zero|| t(h1);
x0=theta3*&a1;
h0=beta0+beta1*&a1;
gamma[4,]=theta3|| x0||zero|| one|| zero|| t(h0);
/* gamma[5,] (Gamma for marg tnie)
mreg linear yreg linear interaction t case.
x0 expression was missing in the following, making gamma[5,] use x0 defined for gamma[4,]
(Gamma for marg tnde). gamma[5,] is the partial derivative of effect[,5] wrt beta1, thus,
it should be (theta2 + theta3 * a1). Note that (a1-a0) is factored out. */
x0=theta2+theta3*&a1; /* Added by @kaz-yos on 2020-04-01 See V2015 p466 Gamma_tnie */
w0=beta1*&a1;
gamma[5,]=zero|| x0|| zero||zero|| beta1|| w0 ;
/* gamma[3,] (Gamma for marg pnie)
mreg linear yreg linear interaction t case.
x1 expression was missing in the following, making gamma[3,] use x1 defined for gamma[2,]
(Gamma for marg pnde). gamma[3,] is the partial derivative of effect[,3] wrt beta1, thus,
it should be (theta2 + theta3 * a0). Note that (a1-a0) is factored out. */
x1=theta2+theta3*&a0; /* Added by @kaz-yos on 2020-04-01 Same as x0 except for the a1 -> a0 change. */
w1=beta1*&a0;
gamma[3,]=zero|| x1|| zero||zero|| beta1|| w1 ;
A=(theta3*&a1+theta3*&a0+theta2);
B=beta0+beta1*(&a0+&a1);
gamma[6,]=theta3||A||zero||one||beta1||B;
				%end;
				%if &yreg=linear & &mreg=logistic & &interaction=true %then %do;
effect[,1]=(theta1+theta3*&m)*(&a1-&a0);
*/CONDITIONAL=MARGINAL NDE*/;
effect[,2]=(theta1+theta3*exp(beta0+beta1*&a0)/(1+exp(beta0+beta1*&a0)))*(&a1-&a0);
*/ CONDITIONAL=MARGINAL TNIE*/;
effect[,3]=(theta2+theta3*&a0)*(exp(beta0+beta1*&a1)/(1+exp(beta0+beta1*&a1))-exp(beta0+beta1*&a0)/(1+exp(beta0+beta1*&a0)));
*/CONDITIONAL=MARGINAL TNDE*/;
effect[,4]=(theta1+theta3*exp(beta0+beta1*&a1)/(1+exp(beta0+beta1*&a1)))*(&a1-&a0);
*/ CONDITIONAL=MARGINAL TNIE*/;
effect[,5]=(theta2+theta3*&a1)*(exp(beta0+beta1*&a1)/(1+exp(beta0+beta1*&a1))-exp(beta0+beta1*&a0)/(1+exp(beta0+beta1*&a0)));
*/te*/;
effect[,6]=effect[,2]+effect[,5];
gamma[1,]=z||&m ;
A=exp(beta0+beta1*&a0);
B=(1+A);
x=theta3*(A*B-A**2)/B**2;
w=theta3*&a0*(A*B-A**2)/B**2;
h=A/B;
gamma[2,]= x|| w || zero|| one|| zero|| t(h);
A=exp(beta0+beta1*&a1);
B=(1+A);
x=theta3*(A*B-A**2)/B**2;
w=theta3*&a1*(A*B-A**2)/B**2;
h=A/B;
gamma[4,]=x|| w || zero|| one|| zero|| t(h);
D=exp(beta0+beta1*&a1);
E=(1+D);
A=exp(beta0+beta1*&a0);
B=(1+A);
x=(theta2+theta3*&a1)*((D*E-D**2)/E**2-(A*B-A**2)/B**2);
w=(theta2+theta3*&a1)*(&a1*(D*E-D**2)/E**2-&a0*(A*B-A**2)/B**2);
h=t(D/E-A/B);
j=&a1*h;
gamma[5,]=x|| w|| zero||zero|| h|| j ;
x=(theta2+theta3*&a0)*((D*E-D**2)/E**2-(A*B-A**2)/B**2);
w=(theta2+theta3*&a0)*(&a1*(D*E-D**2)/E**2-&a0*(A*B-A**2)/B**2);
h=t(D/E-A/B);
j=&a0*h;
gamma[3,]=x|| w|| zero||zero|| h|| j;
A=exp(beta0+beta1*&a0);
B=(1+A);
D=exp(beta0+beta1*&a1);
E=(1+D);
/* Fixed by @kaz-yos. Gamma_te expression based on VV2013 Appendix p14. */
/* Corrected from
x=theta3*(&a1-&a0)*(A*B-B**2)/(B**2)+(theta2+theta3*&a1)*(((D*E-D**2)/(E**2))-((A*B-B**2)/(B**2))); */
x=theta3*(&a1-&a0)*(A*B-A**2)/(B**2)+(theta2+theta3*&a1)*(((D*E-D**2)/(E**2))-((A*B-A**2)/(B**2)));
/* Corrected from
w=&a0*theta3*(&a1-&a0)*(A*B-B**2)/(B**2)+(theta2+theta3*&a1)*(&a0*((D*E-D**2)/(E**2))-&a0*((A*B-B**2)/(B**2))); */
w=&a0*theta3*(&a1-&a0)*(A*B-A**2)/(B**2)+(theta2+theta3*&a1)*(&a1*((D*E-D**2)/(E**2))-&a0*((A*B-A**2)/(B**2)));
s=(&a1-&a0);
t=t(D/E-A/B);
r=(&a1-&a0)*t(A/B)+&a1*t;
gamma[6,]=x||w||zero||s||t||r;
				%end;
				%if &yreg^=linear & &mreg=linear & &interaction=true %then %do;
			 */MARGINAL=CONDITIONAL CDE*/;
x1=(theta1+theta3*&m)*(&a1-&a0);
effect[,1]=exp(x1);
      */MARGINAL=CONDITIONAL NDE*/;
x2=(theta1+theta3*beta0+theta3*beta1*&a0+theta3*theta2*rm)*(&a1-&a0)+(1/2)*tsq*rm*(asq-a1sq);
effect[,2]=exp(x2);
      */MARGINAL=CONDITIONAL NIE*/;
x3=(theta2*beta1+theta3*beta1*&a0)*(&a1-&a0);
effect[,3]=exp(x3);
      */MARGINAL=CONDITIONAL TNDE*/;
x4=(theta1+theta3*beta0+theta3*beta1*&a1+theta3*theta2*rm)*(&a1-&a0)+(1/2)*tsq*rm*(asq-a1sq);
      effect[,4]=exp(x4);
      */ MARGINAL=CONDITIONAL TNIE*/;
x5=(theta2*beta1+theta3*beta1*&a1)*(&a1-&a0);
effect[,5]=exp(x5);
 */ MARGINAL=CONDITIONAL TE*/;
effect[,6]=effect[,2]*effect[,5];
gamma[1,]=z||&m ||zero;
x=theta3*&a0;
h=beta0+beta1*&a0+theta2*s2+theta3*s2*(&a1+&a0);
ts=s2*theta3;
f=theta3*theta2+0.5*(theta3**2)*(&a1+&a0);
gamma[2,]= theta3|| x||zero|| one||ts|| h ||f;
x=theta3*&a1;
h=beta0+beta1*&a1+theta2*s2+theta3*s2*(&a1+&a0);
ts=s2*theta3;
f=theta3*theta2+0.5*theta3**2*(&a1+&a0);
gamma[4,]=theta3|| x|| zero|| one|| ts|| h||f;
x=theta2+theta3*&a1;
w=beta1*&a1;
gamma[5,]=zero|| x||zero||zero|| beta1|| w || zero;
x=theta2+theta3*&a0;
w=beta1*&a0;
gamma[3,]=zero|| x|| zero||zero|| beta1|| w ||zero;
d2pnde=theta3*&a0;
d7pnde=beta0+beta1*&a0+theta2*s2+theta3*s2*(&a1+&a0);
d6pnde=s2*theta3;
d9pnde=theta3*theta2+0.5*(theta3**2)*(&a1+&a0);
d2tnie=theta2+theta3*&a1;
d7tnie=beta1*&a1;
d2=d2pnde+d2tnie;
d6=d6pnde+beta1;
d7=d7pnde+d7tnie;
d9=d9pnde;
gamma[6,]=theta3||d2||zero||one||d6||d7||d9;
				%end;
				%if &yreg^=linear & &mreg=logistic & &interaction=false %then %do;
	  */MARGINAL=CONDITIONAL CDE*/;
	  x1=exp(theta1*(&a1-&a0));
effect[,1]=x1;
      */MARGINAL=CONDITIONAL NDE*/;
effect[,2]=exp((theta1)*(&a1-&a0))*
(1+exp(theta2+beta0+beta1*&a0))/
(1+exp(theta2+beta0+beta1*&a0));
      */MARGINAL=CONDITIONAL NIE*/;
effect[,3]=(
(1+exp(beta0+beta1*&a0))*
(1+exp(theta2+beta0+beta1*&a1))
)/
(
(1+exp(beta0+beta1*&a1))*
(1+exp(theta2+beta0+beta1*&a0))
);
      */MARGINAL=CONDITIONAL TNDE*/;
effect[,4]=exp((theta1)*(&a1-&a0))*
(1+exp(theta2+beta0+beta1*&a1))/
(1+exp(theta2+beta0+beta1*&a1));
      */ MARGINAL=CONDITIONAL TNIE*/;
effect[,5]=(
(1+exp(beta0+beta1*&a0))*
(1+exp(theta2+beta0+beta1*&a1))
)/
(
(1+exp(beta0+beta1*&a1))*(1+exp(theta2+beta0+beta1*&a0))
);
effect[,6]=effect[,2]*effect[,5];
/* Modified by @kaz-yos on 2020-04-01 based on V2015 p475 on Gamma_cde. */
/* z is defined as z=zero||zero||z1||zero||one||zero; above without common factor (a1-a0).*/
/* In the case of &yreg^=linear & &mreg=logistic no common factor multiplication is done. */
/* Confirm looking for "%if (&mreg=logistic & &yreg^=linear) %then %do;" */
/* So it must be included in Gamma_cde. */
/* Changed from
gamma[1,]=z; */
gamma[1,]=z * (&a1-&a0);
A=exp(theta2+beta0+beta1*&a0);
B=(1+exp(theta2+beta0+beta1*&a0));
D=exp(theta2+beta0+beta1*&a0);
E=(1+exp(theta2+beta0+beta1*&a0));
d1nde=A/B-D/E;
d2nde=&a0*(d1nde);
d3nde=0;
d4nde=(&a1-&a0);
d5nde=d1nde;
gamma[2,]= d1nde|| d2nde||d3nde||d4nde||d5nde;
A=exp(theta2+beta0+beta1*&a1);
B=(1+exp(theta2+beta0+beta1*&a1));
D=exp(theta2+beta0+beta1*&a1);
E=(1+exp(theta2+beta0+beta1*&a1));
s=A/B-D/E;
x=&a1*(s);
t=(&a1-&a0);
gamma[4,]=s|| x|| zero||t||s;
A=exp(theta2+beta0+beta1*&a1);
B=(1+exp(theta2+beta0+beta1*&a1));
D=exp(theta2+beta0+beta1*&a0);
E=(1+exp(theta2+beta0+beta1*&a0));
F=exp(beta0+beta1*&a0);
G=(1+exp(beta0+beta1*&a0));
H=exp(beta0+beta1*&a1);
I=(1+exp(beta0+beta1*&a1));
d1nie=F/G-H/I+A/B-D/E;
d2nie=&a0*F/G-&a1*H/I+&a1*A/B-&a0*D/E;
d3nie=0;
d4nie=0;
d5nie=(A/B-D/E);
gamma[5,]=d1nie|| d2nie||d3nie||d4nie||d5nie;
A=exp(theta2+beta0+beta1*&a1);
B=(1+exp(theta2+beta0+beta1*&a1));
D=exp(theta2+beta0+beta1*&a0);
E=(1+exp(theta2+beta0+beta1*&a0));
F=exp(beta0+beta1*&a0);
G=(1+exp(beta0+beta1*&a0));
H=exp(beta0+beta1*&a1);
I=(1+exp(beta0+beta1*&a1));
s=F/G-H/I+A/B-D/E;
x=&a1*F/G-&a0*H/I+&a0*A/B-&a1*D/E;
k=(A/B-D/E);
gamma[3,]=s|| x|| zero||zero||k;
d1=((d1nie)+(d1nde));
d2=((d2nie)+(d2nde));
d3=((d3nie)+(d3nde));
d4=((d4nie)+(d4nde));
d5=((d5nie)+(d5nde));
gamma[6,]=d1||d2||d3||d4||d5;
			%end;
			%if &yreg^=linear & &mreg=logistic & &interaction=true %then %do;
*/MARGINAL CDE*/;
x6=(theta1+theta3*&m)*(&a1-&a0);
effect[,1]=exp(x6);
*/MARGINAL NDE*/;
effect[,2]=exp(theta1*(&a1-&a0))*
(1+exp(theta2+theta3*&a1+beta0+beta1*&a0))/
(1+exp(theta2+theta3*&a0+beta0+beta1*&a0));
*/MARGINAL NIE*/;
effect[,3]=(
(1+exp(beta0+beta1*&a0))*
(1+exp(theta2+theta3*&a0+beta0+beta1*&a1))
)/
(
(1+exp(beta0+beta1*&a1))*
(1+exp(theta2+theta3*&a0+beta0+beta1*&a0))
);
*/ MARGINAL TNDE*/;
effect[,4]=exp(theta1*(&a1-&a0))*
(1+exp(theta2+theta3*&a1+beta0+beta1*&a1))/
(1+exp(theta2+theta3*&a0+beta0+beta1*&a1));
*/ MARGINAL TNIE*/;
effect[,5]=(
(1+exp(beta0+beta1*&a0))*
(1+exp(theta2+theta3*&a1+beta0+beta1*&a1))
)/
(
(1+exp(beta0+beta1*&a1))*
(1+exp(theta2+theta3*&a1+beta0+beta1*&a0))
);
effect[,6]=effect[,2]*effect[,5];

/* Changed by @kaz-yos on 2020-04-01 based on V2015 p474. */
/* z is defined as z=zero||zero||z1||zero||one||zero; above without common factor (a1-a0).*/
/* In the case of &yreg^=linear & &mreg=logistic no common factor multiplication is done. */
/* Confirm looking for "%if (&mreg=logistic & &yreg^=linear) %then %do;" */
/* So it must be included in Gamma_cde. */
/* Changed from
gamma[1,]=z||&m; */
gamma[1,]=(z||&m) * (&a1-&a0);
A=exp(theta2+theta3*&a1+beta0+beta1*&a0);
B=(1+exp(theta2+theta3*&a1+beta0+beta1*&a0));
D=exp(theta2+theta3*&a0+beta0+beta1*&a0);
E=(1+exp(theta2+theta3*&a0+beta0+beta1*&a0));
d1nde=A/B-D/E;
d2nde=&a0*(d1nde);
d3nde=0;
d4nde=(&a1-&a0);
d5nde=d1nde;
d6nde=&a1*A/B-&a0*D/E;
gamma[2,]= d1nde|| d2nde|| d3nde||d4nde||d5nde||d6nde;
A=exp(theta2+theta3*&a1+beta0+beta1*&a1);
B=(1+exp(theta2+theta3*&a1+beta0+beta1*&a1));
D=exp(theta2+theta3*&a0+beta0+beta1*&a1);
E=(1+exp(theta2+theta3*&a0+beta0+beta1*&a1));
s=A/B-D/E;
x=&a1*(s);
t=(&a1-&a0);
h=&a1*A/B-&a0*D/E;
gamma[4,]=s|| x|| zero||t||s|| h;
A=exp(theta2+theta3*&a1+beta0+beta1*&a1);
B=(1+exp(theta2+theta3*&a1+beta0+beta1*&a1));
D=exp(theta2+theta3*&a1+beta0+beta1*&a0);
E=(1+exp(theta2+theta3*&a1+beta0+beta1*&a0));
F=exp(beta0+beta1*&a0);
G=(1+exp(beta0+beta1*&a0));
H=exp(beta0+beta1*&a1);
I=(1+exp(beta0+beta1*&a1));
d1nie=F/G-H/I+A/B-D/E;
d2nie=&a0*F/G-&a1*H/I+&a1*A/B-&a0*D/E;
d3nie=0;
d4nie=0;
d5nie=A/B-D/E;
d6nie=&a1*(A/B-D/E);
gamma[5,]=d1nie|| d2nie|| d3nie||d4nie||d5nie||d6nie;
A=exp(theta2+theta3*&a0+beta0+beta1*&a1);
B=(1+exp(theta2+theta3*&a0+beta0+beta1*&a1));
D=exp(theta2+theta3*&a0+beta0+beta1*&a0);
E=(1+exp(theta2+theta3*&a0+beta0+beta1*&a0));
F=exp(beta0+beta1*&a0);
G=(1+exp(beta0+beta1*&a0));
H=exp(beta0+beta1*&a1);
I=(1+exp(beta0+beta1*&a1));
s=F/G-H/I+A/B-D/E;
x=&a1*F/G-&a0*H/I+&a0*A/B-&a1*D/E;
l=A/B-D/E;
k=&a0*(A/B-D/E);
gamma[3,]=s|| x||zero||zero||l|| k;
d1=(d1nie+d1nde);
d2=(d2nie+d2nde);
d3=(d3nie+d3nde);
d4=(d4nie+d4nde);
d5=(d5nie+d5nde);
d6=(d6nie+d6nde);
gamma[6,]=d1||d2||d3||d4||d5||d6;

				%end;
			%END;*cvar=;
			%if &c^= %then %do;
se=J(1,12);
pvalue=J(1,12);
cil=J(1,12);
ciu=J(1,12);


				%if (&mreg=logistic & &yreg^=linear) %then %do;
					%do j=1 %to 12;
se[,&j]=sqrt(gamma[&j,]*sigma*t(gamma[&j,]));
					%end;
				%end;

				%if &mreg=linear %then %do;
					%do j=1 %to 12;
se[,&j]=sqrt(gamma[&j,]*sigma*t(gamma[&j,]))*abs(&a1-&a0);
					%end;
				%end;
				%if (&mreg=logistic & &yreg=linear) %then %do;
se[,1]=sqrt(gamma[1,]*sigma*t(gamma[1,]))*abs(&a1-&a0);
se[,2]=sqrt(gamma[2,]*sigma*t(gamma[2,]))*abs(&a1-&a0);
se[,3]=sqrt(gamma[3,]*sigma*t(gamma[3,]));
se[,4]=sqrt(gamma[4,]*sigma*t(gamma[4,]))*abs(&a1-&a0);
se[,5]=sqrt(gamma[5,]*sigma*t(gamma[5,]));
se[,6]=sqrt(gamma[6,]*sigma*t(gamma[6,]));
se[,7]=sqrt(gamma[7,]*sigma*t(gamma[7,]))*abs(&a1-&a0);
se[,8]=sqrt(gamma[8,]*sigma*t(gamma[8,]))*abs(&a1-&a0);
se[,9]=sqrt(gamma[9,]*sigma*t(gamma[9,]));
se[,10]=sqrt(gamma[10,]*sigma*t(gamma[10,]))*abs(&a1-&a0);
se[,11]=sqrt(gamma[11,]*sigma*t(gamma[11,]));
se[,12]=sqrt(gamma[12,]*sigma*t(gamma[12,]));
				%end;
				%do j=1 %to 12;

					%if (&yreg=linear ) %then %do;
cil[,&j]=effect[,&j]-1.96*(se[,&j]);
ciu[, &j]=effect[,&j]+1.96*(se[,&j]);
pvalue[,&j] = 2*MIN(1-ABS(probnorm(((effect[,&j]))/(se[,&j]))),ABS(probnorm(((effect[,&j]))/(se[,&j]))));

					%end;
					%if (&yreg^=linear ) %then %do;
					pvalue[,&j] = 2*MIN(1-ABS(probnorm((log(effect[,&j]))/(se[,&j]))),ABS(probnorm((log(effect[,&j]))/(se[,&j]))));

cil[,&j]=exp(log(effect[,&j])-1.96*(se[,&j]));
ciu[, &j]=exp(log(effect[,&j])+1.96*(se[,&j]));
					%end;
				%end;
x=effect;
cname1 = { "effect1" "effect2" "effect3" "effect4" "effect5" "effect6" "effect7" "effect8" "effect9" "effect10" "effect11" "effect12"};
create effect from x [ colname=cname1 ];
append from x;
x=se;
cname1 = { "se1" "se2" "se3" "se4" "se5" "se6" "se7" "se8" "se9" "se10" "se11" "se12"};
create se from x [ colname=cname1 ];
append from x;
x=cil;
cname1 = { "cil1" "cil2" "cil3" "cil4" "cil5" "cil6" "cil7" "cil8" "cil9" "cil10" "cil11" "cil12"};
create cil from x [ colname=cname1 ];
append from x;
x=ciu;
cname1 = { "ciu1" "ciu2" "ciu3" "ciu4" "ciu5" "ciu6" "ciu7" "ciu8" "ciu9" "ciu10" "ciu11" "ciu12"};
create ciu from x [ colname=cname1 ];
append from x;
x=pvalue;
cname1 = { "p1" "p2" "p3" "p4" "p5" "p6" "p7" "p8" "p9" "p10" "p11" "p12"};
create pvalue from x [ colname=cname1 ];
append from x;
			%end;
			%if &c= %then %do;
se=J(1,6);
pvalue=J(1,6);
cil=J(1,6);
ciu=J(1,6);



				%if (&mreg=logistic & &yreg^=linear) %then %do;
					%do j=1 %to 6;
se[,&j]=sqrt(gamma[&j,]*sigma*t(gamma[&j,]));
					%end;
				%end;
				%if (&mreg=linear & &yreg=linear) | (&mreg=linear & &yreg^=linear) %then %do;
					%do j=1 %to 6;
se[,&j]=sqrt(gamma[&j,]*sigma*t(gamma[&j,]))*abs(&a1-&a0);
					%end;
				%end;
				%if (&mreg=logistic & &yreg=linear) %then %do;
se[,1]=sqrt(gamma[1,]*sigma*t(gamma[1,]))*abs(&a1-&a0);
se[,2]=sqrt(gamma[2,]*sigma*t(gamma[2,]))*abs(&a1-&a0);
se[,3]=sqrt(gamma[3,]*sigma*t(gamma[3,]));
se[,4]=sqrt(gamma[4,]*sigma*t(gamma[4,]))*abs(&a1-&a0);
se[,5]=sqrt(gamma[5,]*sigma*t(gamma[5,]));
se[,6]=sqrt(gamma[6,]*sigma*t(gamma[6,]));
				%end;

 				%do j=1 %to 6;
					%if (&yreg=linear ) %then %do;
cil[,&j]=effect[,&j]-1.96*(se[,&j]);
ciu[, &j]=effect[,&j]+1.96*(se[,&j]);
pvalue[,&j] = 2*MIN(1-ABS(probnorm((effect[,&j])/(se[,&j]))),ABS(probnorm((effect[,&j])/(se[,&j]))));
					%end;
					%if (&yreg^=linear ) %then %do;
cil[,&j]=exp(log(effect[,&j])-1.96*(se[,&j]));
ciu[, &j]=exp(log(effect[,&j])+1.96*(se[,&j]));
pvalue[,&j] = 2*MIN(1-ABS(probnorm(log(effect[,&j])/(se[,&j]))),ABS(probnorm(log(effect[,&j])/(se[,&j]))));
					%end;
				%end;
x=effect;
cname1 = { "effect1" "effect2" "effect3" "effect4" "effect5" "effect6"};
create effect from x [ colname=cname1 ];
append from x;
x=se;
cname1 = { "se1" "se2" "se3" "se4" "se5" "se6"};
create se from x [ colname=cname1 ];
append from x;
x=cil;
cname1 = { "cil1" "cil2" "cil3" "cil4" "cil5" "cil6"};
create cil from x [ colname=cname1 ];
append from x;
x=ciu;
cname1 = { "ciu1" "ciu2" "ciu3" "ciu4" "ciu5" "ciu6"};
create ciu from x [ colname=cname1 ];
append from x;
x=pvalue;
cname1 = { "p1" "p2" "p3" "p4" "p5" "p6"};
create pvalue from x [ colname=cname1 ];
append from x;
			%end;
		%end;*other;
	%END;* deltam
***************** causal effects for delta END  *************************;


***************************   DELTA METHOD PROCEDURE -END-  ***************************************************************;



***************************   OUTPUT    *************************** ;


*/LINEAR LINEAR no interaction ;
	%if (&mreg=linear & &interaction=false ) | (&yreg=linear & &mreg=logistic & &interaction=false & &cvar=)  %then %do;
proc iml;
use effect;
read all into effect;
use se;
read all into se;
		%if (&boot= | &boot=false) %then %do;
use pvalue;
read all into pvalue;
use cil;
read all into cil;
use ciu;
read all into ciu;

		%end;
		%if (&boot^= & &boot^=false) %then %do;
use ci;
read all into ci;

x= t(effect)||(se)||ci ;
cname1 = { "Estimate"  "s.e." "95% CI lower" "95% CI upper" };
		%END;
		%if (&boot= | &boot=false) %then %do;
		%if &yreg=linear %then %do;
x= t(effect)||t(se)|| t(pvalue)||t(cil) || t(ciu) ;
cname1 = { "Estimate"  "s.e." "p-value" "95% CI lower" "95% CI upper"  };
		%end;
		%if &yreg^=linear %then %do;
x= t(effect)|| t(pvalue)||t(cil) || t(ciu) ;
cname1 = { "Estimate" "p-value" "95% CI lower" "95% CI upper"  };
		%end;
		%END;
create x3 from x [ colname=cname1 ];
append from x;
name='cde=nde' || 'nie'||'total effect' ;
name=t(name);
cname2= {"Effect"};
create x4 from name [colname=cname2];
append from name;
quit;
	%end;


*/LINEAR LINEAR interaction ;
	%if (&interaction=true) | (&mreg=logistic & &interaction=false & &cvar^=) |(&yreg^=linear & &mreg=logistic & &interaction=false) %then %do;

proc iml;
use effect;
read all into effect;
use se;
read all into se;

		%if &output=full & &c^= & &cvar^= %then %do;
name='marginal cde' || 'marginal pnde'||'marginal pnie'||'marginal tnde'||'marginal tnie'||'marginal total effect'||'conditional cde' || 'conditional pnde'||'conditional pnie'||'conditional tnde'||'conditional tnie'||'conditional total effect' ;
name=t(name);
cname2= {"effect"};
create x4 from name [colname=cname2];
append from name;
			%if (&boot^= & &boot^=false)  %then %do;
use ci;
read all into ci;
x= t(effect)||(se)||ci ;
cname1 = { "Estimate"  "s.e." "95% CI lower" "95% CI upper" };
			%END;
			%if (&boot= | &boot=false)  %then %do;
use cil;
read all into cil;
use ciu;
read all into ciu;
use pvalue;
read all into pvalue;
			%if &yreg=linear %then %do;
x= t(effect)||t(se)|| t(pvalue) ||t(cil) || t(ciu) ;
cname1 = { "Estimate"  "s.e." "p-value" "95% CI lower" "95% CI upper"  };
			%end;
			%if &yreg^=linear %then %do;
x= t(effect)|| t(pvalue) ||t(cil) || t(ciu) ;
cname1 = { "Estimate" "p-value" "95% CI lower" "95% CI upper"  };
			%end;
			%END;
create x3 from x [ colname=cname1 ];
append from x;
quit;
		%end;
		%if (&output=full & &cvar=) | (&output=full & &cvar^= & &c=)  %then %do;
name='cde'|| 'pnde'||'pnie'||'tnde'||'tnie'||'total effect';
name=t(name);
cname2= {"effect"};
create x4 from name [colname=cname2];
append from name;
			%if (&boot^= & &boot^=false) %then %do;
use ci;
read all into ci;
x= t(effect)||(se)||ci ;
cname1 = { "Estimate"  "s.e." "95% CI lower" "95% CI upper" };
			%END;
			%if (&boot= | &boot=false) %then %do;
use cil;
read all into cil;
use ciu;
read all into ciu;
use pvalue;
read all into pvalue;
use pvalue;
read all into pvalue;
%if &yreg=linear %then %do;
x= t(effect)||t(se)|| t(pvalue) ||t(cil) || t(ciu) ;
cname1 = { "Estimate"  "s.e." "p-value" "95% CI lower" "95% CI upper"  };
			%end;
			%if &yreg^=linear %then %do;
x= t(effect)|| t(pvalue) ||t(cil) || t(ciu) ;
cname1 = { "Estimate" "p-value" "95% CI lower" "95% CI upper"  };
			%end;
			%END;
create x3 from x [ colname=cname1 ];
append from x;
quit;
		%end;
		%if &output^=full %then %do;
use effect;
read all into effect;
effect1=effect[,1:2];
effect2=effect[,5:6];
effect=effect1||effect2;

			%if (&boot= | &boot=false)  %then %do;
use se;
read all into se;
se1=se[,1:2];
se2=se[,5:6];
se=se1||se2;
use pvalue;
read all into pvalue;
pval1=pvalue[,1:2];
pval2=pvalue[,5:6];
pvalue=pval1||pval2;
use cil;
read all into cil;
cil1=cil[,1:2];
cil2=cil[,5:6];
cil=cil1||cil2;
use ciu;
read all into ciu;
ciu1=ciu[,1:2];
ciu2=ciu[,5:6];
ciu=ciu1||ciu2;
			%end;
			%if (&boot^= & &boot^=false)  %then %do;
use se;
read all into se;
se1=se[1:2,];
se2=se[5:6,];
se=se1//se2;
use ci;
read all into ci;
ci1=ci[1:2,];
ci2=ci[5:6,];
ci=ci1//ci2;
x= t(effect)||(se)||ci ;
cname1 = { "Estimate"  "s.e." "95% CI lower" "95% CI upper" };
			%END;
			%if (&boot= | &boot=false) %then %do;
%if &yreg=linear %then %do;
x= t(effect)||t(se)|| t(pvalue) ||t(cil) || t(ciu) ;
cname1 = { "Estimate"  "s.e." "p-value" "95% CI lower" "95% CI upper"  };
			%end;
			%if &yreg^=linear %then %do;
x= t(effect)|| t(pvalue) ||t(cil) || t(ciu) ;
cname1 = { "Estimate" "p-value" "95% CI lower" "95% CI upper"  };
			%end;
			%END;
create x3 from x [ colname=cname1 ];
append from x;
run;
name='cde'|| 'nde'||'nie'||'total effect';
name=t(name);
cname2= {"Effect"};
create x4 from name [colname=cname2];
append from name;
quit;
		%END;
	%end;*end interaction;




DATA x5;
MERGE x4 x3;
proc print data=x5;
quit;


%if &output^=full %then %do;
proc iml;
name= 'proportion mediated';
cname2= {"Effect"};
create xbis from name [colname=cname2];
append from name;
use x3;
read all into VB;
	%if &interaction=true %then %do;
nde=VB[2,1];
nie=VB[3,1];
	%end;
	%if &interaction=false %then %do;
nde=VB[1,1];
nie=VB[2,1];
	%end;

%if  (&mreg=logistic & &interaction=false & &cvar^=) |(&yreg^=linear & &mreg=logistic & &interaction=false) %then %do;

nde=VB[2,1];
nie=VB[3,1];
%end;

	%if &yreg=linear %then %do;
pm=nie/(nde+nie);
	%end;
	%if &yreg^=linear %then %do;
pm=(nde*(nie-1))/(nde*nie-1);
	%end;
cname1 = { "Estimate"};
create xtris from pm [ colname=cname1 ];
append from pm;
run;
quit;
DATA x6;
MERGE xbis xtris;
proc print data=x6;
quit;
%end;


***************************   OUTPUT -END-   *************************** ;


proc datasets library=work nolist;
delete data1 data2 data3 out1 out2 x3 x4 x5;
run;

%if &output^=full %then %do;
proc datasets library=work nolist;
delete x6 xbis xtris;
run;
%end;

%if &boot= %then %do;
proc datasets library=work nolist;
delete cil ciu effect pvalue se;
run;
%end;

%if &cvar^= %then %do;
 %do i = 1 %to &nc;
 proc datasets library=work nolist;
delete data2&i data2new&i;
run;
 %end;
 %end;

%if &yreg=loglinear | &yreg=poisson | &yreg=negbin | &yreg=survAFT_exp | &yreg=survAFT_weibull | &yreg=survCox %then %do;
proc datasets library=work nolist;
delete gmcovb gmparms;
run;
%end;

%if (&boot^= & &boot^=false) %then %do;
	%do t = 1 %to &n;
 proc datasets library=work nolist;
delete data1&t out1&t out2&t pmethod&t;
run;
%end;
proc datasets library=work nolist;
delete effect se ci Bootdata;
run;
%end;



%mend;
