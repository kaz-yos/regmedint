# regmedint (developmental repo) <img src="hex.png" align="right" height="139"/>
[![Travis-CI Build Status](https://travis-ci.org/kaz-yos/regmedint.svg?branch=master)](https://travis-ci.org/kaz-yos/regmedint)
[![](http://www.r-pkg.org/badges/version/regmedint)](http://www.r-pkg.org/pkg/regmedint)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/regmedint)](http://www.r-pkg.org/pkg/regmedint)

This is an R reimplementation of the regression-based causal mediation analysis methods as implemented in the SAS macro by Valeri and VanderWeele (2013 and 2015). The original is found at

# Status

The eAppendix for SAS macro for causal mediation analysis with survival data (Valeri & VanderWeele 2015) describes the following grid of models.

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">


<colgroup>
<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />
</colgroup>
<thead>
<tr>
<th scope="col" class="org-left">yreg \\ mreg</th>
<th scope="col" class="org-left">linear</th>
<th scope="col" class="org-left">logistic</th>
</tr>
</thead>

<tbody>
<tr>
<td class="org-left">linear</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">logistic</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">loglinear</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">poisson</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">negbin</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">survCox</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">survAFT exp</td>
<td class="org-left">1</td>
<td class="org-left">1</td>
</tr>


<tr>
<td class="org-left">survAFT weibull</td>
<td class="org-left">1</td>
<td class="org-left">1</td>
</tr>
</tbody>
</table>

The current implementation status:

1.  None (left blank)
2.  SAS reference results created
3.  R implementation ongoing
4.  R implementation complete


<a id="org30ff663"></a>

# Formulas

-   V2015: VanderWeele (2015) Explanation in Causal Inference.
-   VV2013A: Valeri & VanderWeele (2013) Appendix
-   VV2015A: Valeri & VanderWeele (2015) Appendix

As seen below and <https:/github.com/harvard-P01/causalMediation/blob/master/R/NDE.R/>, there are only four formulas.


<a id="orgff32dd5"></a>

## Effect formulas

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">


<colgroup>
<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />
</colgroup>
<thead>
<tr>
<th scope="col" class="org-left">yreg \\ mreg</th>
<th scope="col" class="org-left">linear</th>
<th scope="col" class="org-left">logistic</th>
</tr>
</thead>

<tbody>
<tr>
<td class="org-left">linear</td>
<td class="org-left">V2015 p466 Proposition 2.3</td>
<td class="org-left">V2015 p471 Proposition 2.5</td>
</tr>


<tr>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">logistic</td>
<td class="org-left">V2015 p468 Proposition 2.4</td>
<td class="org-left">V2015 p473 Proposition 2.6</td>
</tr>


<tr>
<td class="org-left">loglinear</td>
<td class="org-left">VV2013A p8 Use Proposition 2.4</td>
<td class="org-left">VV2013A p8 Use Proposition 2.6</td>
</tr>


<tr>
<td class="org-left">poisson</td>
<td class="org-left">VV2013A p8 Use Proposition 2.4</td>
<td class="org-left">VV2013A p8 Use Proposition 2.6</td>
</tr>


<tr>
<td class="org-left">negbin</td>
<td class="org-left">VV2013A p8 Use Proposition 2.4</td>
<td class="org-left">VV2013A p8 Use Proposition 2.6</td>
</tr>


<tr>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">survCox</td>
<td class="org-left">V2015 p496 Proposition 4.4 (Use 2.4)</td>
<td class="org-left">V2015 p499 Proposition 4.6 (Use 2.6)</td>
</tr>


<tr>
<td class="org-left">survAFT exp</td>
<td class="org-left">V2015 p494 Proposition 4.1 (Use 2.4)</td>
<td class="org-left">V2015 p495 Proposition 4.3 (Use 2.6)</td>
</tr>


<tr>
<td class="org-left">survAFT weibull</td>
<td class="org-left">V2015 p494 Proposition 4.1 (Use 2.4)</td>
<td class="org-left">V2015 p495 Proposition 4.3 (Use 2.6)</td>
</tr>
</tbody>
</table>

Note the loglinear outcome model means log link with binomial error.


<a id="org1aaad20"></a>

## Standard error formulas

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">


<colgroup>
<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />
</colgroup>
<thead>
<tr>
<th scope="col" class="org-left">yreg \\ mreg</th>
<th scope="col" class="org-left">linear</th>
<th scope="col" class="org-left">logistic</th>
</tr>
</thead>

<tbody>
<tr>
<td class="org-left">linear</td>
<td class="org-left">V2015 p466 Proposition 2.3</td>
<td class="org-left">V2015 p471 Proposition 2.5</td>
</tr>


<tr>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">logistic</td>
<td class="org-left">V2015 p468 Proposition 2.4</td>
<td class="org-left">V2015 p473 Proposition 2.6</td>
</tr>


<tr>
<td class="org-left">loglinear</td>
<td class="org-left">VV2013A p8 Use Proposition 2.4</td>
<td class="org-left">VV2013A p8 Use Proposition 2.6</td>
</tr>


<tr>
<td class="org-left">poisson</td>
<td class="org-left">VV2013A p8 Use Proposition 2.4</td>
<td class="org-left">VV2013A p8 Use Proposition 2.6</td>
</tr>


<tr>
<td class="org-left">negbin</td>
<td class="org-left">VV2013A p8 Use Proposition 2.4</td>
<td class="org-left">VV2013A p8 Use Proposition 2.6</td>
</tr>


<tr>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">survCox</td>
<td class="org-left">V2015 p496 Use Proposition 2.4</td>
<td class="org-left">V2015 p499 Use Proposition 2.6</td>
</tr>


<tr>
<td class="org-left">survAFT exp</td>
<td class="org-left">V2015 p494 Use Proposition 2.4</td>
<td class="org-left">V2015 p495 Use Proposition 2.6</td>
</tr>


<tr>
<td class="org-left">survAFT weibull</td>
<td class="org-left">V2015 p494 Use Proposition 2.4</td>
<td class="org-left">V2015 p495 Use Proposition 2.6</td>
</tr>
</tbody>
</table>


<a id="org49f67de"></a>

# Implementation progress

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">


<colgroup>
<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />
</colgroup>
<thead>
<tr>
<th scope="col" class="org-left">mreg</th>
<th scope="col" class="org-left">yreg</th>
<th scope="col" class="org-left">est</th>
<th scope="col" class="org-left">se</th>
</tr>
</thead>

<tbody>
<tr>
<td class="org-left">linear</td>
<td class="org-left">linear</td>
<td class="org-left">calc<sub>myreg</sub><sub>mreg</sub><sub>linear</sub><sub>yreg</sub><sub>linear</sub><sub>est</sub></td>
<td class="org-left">calc<sub>myreg</sub><sub>mreg</sub><sub>linear</sub><sub>yreg</sub><sub>linear</sub><sub>se</sub></td>
</tr>


<tr>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">implemented, need cvar=NULL handling</td>
<td class="org-left">implemented, need cvar=NULL handling</td>
</tr>


<tr>
<td class="org-left">linear</td>
<td class="org-left">logistic</td>
<td class="org-left">calc<sub>myreg</sub><sub>mreg</sub><sub>linear</sub><sub>yreg</sub><sub>logistic</sub><sub>est</sub></td>
<td class="org-left">calc<sub>myreg</sub><sub>mreg</sub><sub>linear</sub><sub>yreg</sub><sub>logistic</sub><sub>se</sub></td>
</tr>


<tr>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">implemented, need cvar=NULL handling</td>
<td class="org-left">se<sub>pm</sub> missing, need cvar=NULL handling</td>
</tr>


<tr>
<td class="org-left">logistic</td>
<td class="org-left">linear</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">logistic</td>
<td class="org-left">logistic</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>
</tbody>
</table>


<a id="org17e78e0"></a>

# TODOs


<a id="orge2bce86"></a>

## 2020-03-15


<a id="org86422e2"></a>

### TODO Complete implementation of calc<sub>myreg</sub><sub>mreg</sub><sub>linear</sub><sub>yreg</sub><sub>logistic</sub>

1.  DONE Expand tests to include interactions (software check)

2.  DONE Test changes in estiamtes as m<sub>cde</sub> and c<sub>cond</sub> change (math check)

3.  TODO Complete input UI -> output UI


<a id="orgc913806"></a>

### TODO Copy test for calc<sub>myreg</sub><sub>mreg</sub><sub>linear</sub><sub>yreg</sub><sub>linear</sub>


<a id="orgba95ac0"></a>

# Design


<a id="org48490a3"></a>

## regmedint UI function


<a id="org7c6c2bf"></a>

### new<sub>regmedint</sub> internal constructor

1.  fit<sub>mreg</sub>

2.  fit<sub>yreg</sub>

3.  calc<sub>myreg</sub> calls a specialized worker function

    1.  eg calc<sub>myreg</sub><sub>mreg</sub><sub>linear</sub><sub>yreg</sub><sub>linear</sub>

        1.  calc<sub>myreg</sub><sub>mreg</sub><sub>linear</sub><sub>yreg</sub><sub>linear</sub><sub>est</sub> returns an effect calculator

        2.  calc<sub>myreg</sub><sub>mreg</sub><sub>linear</sub><sub>yreg</sub><sub>linear</sub><sub>se</sub> returns an se calculator


<a id="org604fe8d"></a>

## regmedint object


<a id="orga90f9cf"></a>

### mreg<sub>fit</sub> mediator regression model object as is


<a id="orgcf08e09"></a>

### yreg<sub>fit</sub> outcome regression model object as is


<a id="orgb7cad40"></a>

### myreg<sub>funs</sub> list

1.  est<sub>fun</sub>: (a0,a1,m<sub>cde,c</sub><sub>cond</sub>) -> (cde,pnde,tnie,tnde,pnie,te,pm)

2.  se<sub>fun</sub>: (a0,a1,m<sub>cde,c</sub><sub>cond</sub>) -> se for (cde,pnde,tnie,tnde,pnie,te,pm)


<a id="orgc33f581"></a>

### args preserves arguments given to the UI


<a id="org937fed6"></a>

## methods for regmedint


<a id="org480b94f"></a>

### print.regmedint


<a id="orge2d951d"></a>

### summary.regmedint: regmedint -> regmedint<sub>summary</sub>

1.  coef.regmedint<sub>summary</sub>: regmedint<sub>summary</sub> -> data.frame of (coef, se, ci, p)


<a id="org7c91457"></a>

### coef.regmedint: regmedint -> vector (cde,pnde,tnie,tnde,pnie,te,pm)


<a id="orge68fa47"></a>

### confint.regmedint: regmedint -> matrix of (lower,upper)


<a id="org44ae55c"></a>

# Similar or related R projects

-   mediation: <https:/cran.r-project.org/web/packages/mediation/index.html/>
-   medflex: <https:/cran.r-project.org/web/packages/medflex/index.html/>
-   intmed: <https:/cran.r-project.org/web/packages/intmed/index.html/>
-   mediator: <https:/github.com/GerkeLab/mediator/>
