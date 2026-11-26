{smcl}
{* *! version 1.0.0  26nov2026}{...}
{vieweralsosee "stcoxem" "help stcoxem"}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "[ST] stcox postestimation" "help stcox_postestimation"}{...}
{viewerjumpto "Postestimation commands" "stcoxem_postestimation##description"}{...}
{viewerjumpto "predict" "stcoxem_postestimation##predict"}{...}
{viewerjumpto "Remarks" "stcoxem_postestimation##remarks"}{...}
{viewerjumpto "Examples" "stcoxem_postestimation##examples"}{...}
{title:Title}

{phang}
{bf:stcoxem postestimation} {hline 2} Postestimation tools for stcoxem


{marker description}{...}
{title:Postestimation commands}

{pstd}
{cmd:stcoxem} stores its results in the layout of
{cmd:stcox, shared() frailty(gamma)}. The postestimation commands are
therefore {cmd:stcox}'s own; they run unchanged and are documented in
{manhelp stcox_postestimation ST:stcox postestimation}. The following are
available:

{synoptset 20}{...}
{p2coldent:Command}Description{p_end}
{synoptline}
{synopt:{helpb stcox postestimation##estatcon:estat concordance}}compute
the concordance probability{p_end}
{synopt:{helpb stcox diagnostics:estat phtest}}test
proportional-hazards assumption based on Schoenfeld residuals{p_end}
{synopt:{helpb estat gofplot}}goodness-of-fit
plot{p_end}
{synopt:{helpb stcurve}}plot the survivor, hazard, and cumulative hazard
functions{p_end}
{synopt:{helpb estat ic}}Akaike's and Schwarz's Bayesian information criteria
(AIC and BIC){p_end}
{synopt:{helpb estat summarize}}summary statistics for the estimation
sample{p_end}
{synopt:{helpb estat vce}}variance-covariance matrix of the estimators
(VCE){p_end}
{synopt:{helpb estimates}}cataloging estimation results{p_end}
{synopt:{helpb lincom}}point estimates, standard errors, testing, and
inference for linear combinations of coefficients{p_end}
{synopt:{helpb lrtest}}likelihood-ratio test{p_end}
{synopt:{helpb margins}}marginal means, predictive margins, marginal effects,
and average marginal effects{p_end}
{synopt:{helpb marginsplot}}graph the results from margins{p_end}
{synopt:{helpb nlcom}}point estimates, standard errors, testing, and
inference for nonlinear combinations of coefficients{p_end}
{synopt:{helpb stcoxem postestimation##predict:predict}}predictions,
residuals, influence statistics, and other diagnostic measures{p_end}
{synopt:{helpb predictnl}}point estimates, standard errors, testing, and
inference for generalized predictions{p_end}
{synopt:{helpb test}}Wald tests of simple and composite linear
hypotheses{p_end}
{synopt:{helpb testnl}}Wald tests of nonlinear hypotheses{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
{cmd:estat concordance}, as after {cmd:stcox}, is not available for data with
delayed entry or multiple records per subject.{p_end}
{p 4 6 2}
{cmd:lrtest} and {cmd:estat ic} use {cmd:e(ll)}, which is on {cmd:stcox}'s
scale; see {help stcoxem##rem_ll:stcoxem}.{p_end}
{p 4 6 2}
Commands that need {cmd:e(V)} ({cmd:test}, {cmd:lincom}, {cmd:estat vce},
{cmd:estat phtest}, standard errors from {cmd:margins}, ...) are not available
after {cmd:vce(none)}.{p_end}


{marker predict}{...}
{title:Syntax for predict}

{p 8 19 2}
{cmd:predict} {dtype} {newvar} {ifin} [{cmd:,} {it:sv_statistic}
{opt atfrail:ty}[{cmd:(}{varname}{c |}{it:#}{cmd:)}]
{opt nooff:set} {opt part:ial}]

{p 8 19 2}
{cmd:predict} {dtype} {c -(}{it:stub}{cmd:*}|{it:newvarlist}{c )-} {ifin}
{cmd:,} {it:mv_statistic} [{opt part:ial}]

{synoptset 20 tabbed}{...}
{synopthdr:sv_statistic}
{synoptline}
{syntab:Main}
{synopt:{opt hr}}predicted hazard ratio, also known as the relative hazard;
the default{p_end}
{synopt:{opt xb}}linear prediction xb{p_end}
{synopt:{opt stdp}}standard error of the linear prediction{p_end}
{synopt:{opt bases:urv}}baseline survivor function{p_end}
{synopt:{opt basec:hazard}}baseline cumulative hazard function{p_end}
{synopt:{opt basehc}}baseline hazard contributions{p_end}
{synopt:{opt mg:ale}}martingale-like residuals{p_end}
{synopt:{opt csn:ell}}Cox-Snell residuals{p_end}
{synopt:{opt dev:iance}}deviance residuals{p_end}
{synopt:{opt ld:isplace}}likelihood displacement values{p_end}
{synopt:{opt lm:ax}}LMAX measures of influence{p_end}
{synopt:{opt eff:ects}}log frailties{p_end}
{synoptline}

{synopthdr:mv_statistic}
{synoptline}
{syntab:Main}
{synopt:{opt sco:res}}efficient score residuals{p_end}
{synopt:{opt esr}}synonym for {opt scores}{p_end}
{synopt:{opt dfb:eta}}DFBETA measures of influence{p_end}
{synopt:{opt sch:oenfeld}}Schoenfeld residuals{p_end}
{synopt:{opt sca:ledsch}}scaled Schoenfeld residuals{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
{opt atfrailty} and {opt atfrailty()} are allowed with {opt basesurv},
{opt basechazard}, and {opt basehc}, as after {cmd:stcox}.{p_end}

{pstd}
The statistics, their defaults, and the options are those of {cmd:predict}
after {cmd:stcox, shared()}; see
{help stcox postestimation##predict:predict after stcox}. In particular, as
after {cmd:stcox, shared()}: {opt hr}, {opt xb}, {opt stdp}, and the baseline
functions are computed with the log frailty set to 0 (u=1) and, for the
baseline functions, x=0; the residuals and influence statistics use the
estimated frailties; {opt mgale}, {opt csnell}, and {opt deviance} are
cumulative per subject (stored on the subject's last record) unless
{opt partial} is specified; {opt effects} is the estimated log frailty
ln E(u | data) of each record's group.


{marker remarks}{...}
{title:Remarks}

{pstd}
{bf:How predict works.} For fits without weights, {cmd:e(predict)} is
{cmd:stcox_p}, {cmd:stcox}'s own program. Statistics that need the baseline
hazard or residuals ({opt basechazard}, {opt mgale}, {opt effects},
{opt scores}, {opt schoenfeld}, ...) are computed by {cmd:stcox_p} refitting
{cmd:stcox}'s penalized frailty model at {cmd:e(b)} and {cmd:e(theta)}. At
these values the refit lands on {cmd:stcoxem}'s solution: the frailties, the
baseline cumulative hazard, and the residuals agree with {cmd:stcoxem}'s own
({opt ufrailty()}, {opt generate()}) to about 1e-10, and with those of
{cmd:stcox, shared()} to the convergence tolerance.

{pstd}
The refit has one parameter per {opt shared()} group. With very many groups
(hundreds of thousands), these statistics may therefore be slow or exceed
memory, even when {cmd:stcoxem} itself is fast. {opt hr}, {opt xb}, and
{opt stdp} need no refit. The frailties and the baseline are also available
directly from {cmd:stcoxem}'s {opt ufrailty()} and {opt generate()} options.

{pstd}
{bf:Weighted fits.} With {cmd:pweight}s, {cmd:e(predict)} is
{cmd:stcoxem_p}, because {cmd:stcox}'s weighted frailty refit is a different
estimator. {opt effects} and {opt atfrailty} use {cmd:stcoxem}'s frailties
(recomputed at {cmd:e(b)} and {cmd:e(theta)}). The other refit statistics
are computed by {cmd:stcox_p} for the Cox model with ln E(u) as an offset,
evaluated at {cmd:e(b)} without iteration; this is the same baseline hazard
and the same residuals. The results equal those of an unweighted fit on the
data expanded by the weights.

{pstd}
{bf:Statistics that use e(V).} {opt stdp}, {opt dfbeta}, {opt ldisplace},
{opt lmax}, {cmd:estat phtest}, and all Wald-type inference use
{cmd:stcoxem}'s {cmd:e(V)}. Under {cmd:vce(oim)} this differs by a few
percent from the variance of {cmd:stcox, shared()}, so these results differ
from {cmd:stcox}'s by the same order. As after {cmd:stcox}: {opt ldisplace}
and {opt lmax} are not available after {cmd:vce(robust)} or
{cmd:vce(cluster)}; {opt dfbeta} is not available after {cmd:vce(none)} or
{cmd:vce(bootstrap)}; {opt stdp} is not available after {cmd:vce(none)}.

{pstd}
{bf:Data changes.} As after {cmd:stcox}, {cmd:predict} refuses (error 459)
if the variables used in the estimation have changed since the fit (checked
with {cmd:e(datasignature)}). Adding new variables is allowed.

{pstd}
{bf:Replay.} Typing {cmd:stcoxem} redisplays the results in {cmd:stcoxem}'s
format; typing {cmd:stcox} redisplays them in {cmd:stcox}'s format.
{cmd:estimates store}, {cmd:estimates restore}, and {cmd:estimates table}
work as usual.


{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. webuse catheter}{p_end}
{phang2}{cmd:. stset time, failure(infect)}{p_end}
{phang2}{cmd:. stcoxem age female, shared(patient)}{p_end}

{pstd}Log frailties and martingale residuals{p_end}
{phang2}{cmd:. predict nu, effects}{p_end}
{phang2}{cmd:. predict mg, mgale}{p_end}

{pstd}Baseline cumulative hazard at the estimated frailty{p_end}
{phang2}{cmd:. predict H, basechazard atfrailty}{p_end}

{pstd}Test of proportional hazards{p_end}
{phang2}{cmd:. estat phtest, detail}{p_end}

{pstd}Survivor function for women and men at mean age{p_end}
{phang2}{cmd:. stcurve, survival at(female=0) at(female=1)}{p_end}

{pstd}LR test of a covariate{p_end}
{phang2}{cmd:. estimates store full}{p_end}
{phang2}{cmd:. stcoxem age, shared(patient)}{p_end}
{phang2}{cmd:. lrtest full .}{p_end}

{pstd}Predictive margins{p_end}
{phang2}{cmd:. stcoxem age i.female, shared(patient)}{p_end}
{phang2}{cmd:. margins female}{p_end}


{title:Also see}

{psee}
Help: {help stcoxem}, {manhelp stcox_postestimation ST:stcox postestimation}
{p_end}
