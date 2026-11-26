{smcl}
{* *! version 1.0.0  26nov2026}{...}
{vieweralsosee "stcoxem postestimation" "help stcoxem_postestimation"}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "[ST] stcox" "help stcox"}{...}
{vieweralsosee "[ST] stcox postestimation" "help stcox_postestimation"}{...}
{vieweralsosee "[ST] stset" "help stset"}{...}
{viewerjumpto "Syntax" "stcoxem##syntax"}{...}
{viewerjumpto "Description" "stcoxem##description"}{...}
{viewerjumpto "Options" "stcoxem##options"}{...}
{viewerjumpto "Weights" "stcoxem##weights"}{...}
{viewerjumpto "Remarks" "stcoxem##remarks"}{...}
{viewerjumpto "Examples" "stcoxem##examples"}{...}
{viewerjumpto "Stored results" "stcoxem##results"}{...}
{viewerjumpto "Methods and formulas" "stcoxem##methods"}{...}
{viewerjumpto "References" "stcoxem##references"}{...}
{viewerjumpto "Author" "stcoxem##author"}{...}
{title:Title}

{phang}
{bf:stcoxem} {hline 2} Cox proportional hazards model with shared gamma frailty,
fit by the EM algorithm


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:stcoxem} [{varlist}] {ifin}{cmd:,} {opth sh:ared(varname)} [{it:options}]

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{p2coldent:* {opth sh:ared(varname)}}variable identifying the shared-frailty
groups{p_end}
{synopt:{opth str:ata(varname)}}stratify on {it:varname}: a separate baseline
hazard within each stratum{p_end}
{synopt:{opt forceshared}}allow delayed entry or gaps{p_end}

{syntab:SE/Robust}
{synopt:{opth vce(vcetype)}}{it:vcetype} may be {opt oim} (default), {opt r:obust},
{opt cl:uster} {it:clustvar}, {opt boot:strap}, or {opt none}{p_end}
{synopt:{opt reps(#)}}number of bootstrap replications; default is
{cmd:reps(50)}{p_end}
{synopt:{opt bsseed(#)}}random-number seed for {cmd:vce(bootstrap)}{p_end}
{synopt:{opt bsw:armstart}}start each bootstrap replication from the
estimates{p_end}
{synopt:{opt nopl:se}}report the observed-information standard error of
theta instead of the profile-likelihood one{p_end}
{synopt:{opt pl:ci}}report a profile-likelihood confidence interval for
theta{p_end}
{synopt:{opt nolr:test}}do not perform the likelihood-ratio test of
theta=0{p_end}

{syntab:Reporting}
{synopt:{opt level(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synopt:{opt nohr}}report coefficients, not hazard ratios{p_end}
{synopt:{opt nohead:er}}suppress the header{p_end}
{synopt:{opt nocoef}}suppress the coefficient table{p_end}
{synopt:{it:{help stcoxem##display_options:display_options}}}control columns
and column formats, row spacing, line width, and factor-variable
labeling{p_end}

{syntab:Maximization}
{synopt:{opt max:iter(#)}}maximum number of EM iterations; default is
{cmd:maxiter(}{ccl maxiter}{cmd:)}{p_end}
{synopt:{opt tol:erance(#)}}tolerance for the parameters; default is
{cmd:tolerance(1e-8)}{p_end}
{synopt:{opt ltol:erance(#)}}tolerance for the log likelihood; default is
{cmd:ltolerance(1e-7)}{p_end}
{synopt:{opt cox:iter(#)}}Newton steps for the coefficients per EM
iteration; default is {cmd:coxiter(1)}{p_end}
{synopt:{opt accel:erate(method)}}EM acceleration: {opt squarem} (default),
{opt aitken}, or {opt none}{p_end}
{synopt:{opt noai:tken}}synonym for {cmd:accelerate(none)}{p_end}
{synopt:{opt nolo:g}}suppress the iteration log{p_end}
{synopt:{opt dots}}display one dot per EM iteration instead of the
log{p_end}

{syntab:Generated variables}
{synopt:{opth gen:erate(newvar)}}create {it:newvar} containing the baseline
cumulative hazard contribution of each record{p_end}
{synopt:{opth uf:railty(newvar)}}create {it:newvar} containing the estimated
frailty of each record's group{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
* {opt shared()} is required.{p_end}
{p 4 6 2}
You must {cmd:stset} your data before using {cmd:stcoxem}; see
{manhelp stset ST}.{p_end}
{p 4 6 2}
{it:varlist} may contain factor variables and time-series operators; see
{help fvvarlist} and {help tsvarlist}. {it:varlist} may be omitted, which fits
a model with frailty only.{p_end}
{p 4 6 2}
{cmd:pweight}s are specified with {cmd:stset}, not with {cmd:stcoxem}; see
{help stcoxem##weights:Weights}.{p_end}
{p 4 6 2}
See {help stcoxem postestimation} for features available after
estimation.{p_end}

{p 4 4 2}
Typing {cmd:stcoxem} without arguments redisplays the last estimates;
{opt level()}, {opt nohr}, {opt noheader}, {opt nocoef}, and the
{it:display_options} may be specified.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:stcoxem} fits a Cox proportional hazards model with a shared gamma
frailty,

{p 8 8 2}
h(t | x, u) = u h0(t) exp(xb),{space 4}u ~ Gamma(mean 1, variance theta),

{pstd}
in which all records of a group defined by {opt shared()} share one
unobserved frailty u, and the baseline hazard h0(t) is left unspecified.
This is the model fit by {cmd:stcox, shared() frailty(gamma)}. {cmd:stcoxem}
fits it by the expectation-maximization (EM) algorithm of
{help stcoxem##klein1992:Klein (1992)}, implemented in Mata, instead of the
penalized partial likelihood used by {cmd:stcox}. The estimates of b, theta,
the frailties, and the baseline hazard agree with those of {cmd:stcox}; the
EM algorithm is written for large data (hundreds of thousands of groups),
where each iteration costs one pass over the risk sets.

{pstd}
{cmd:stcoxem} also allows what {cmd:stcox} does not allow together with
{opt shared()}: {opt strata()}, {opt vce(robust)}, {opt vce(cluster)}, and
{cmd:pweight}s.

{pstd}
The results are stored in the layout of {cmd:stcox, shared()}, so
{cmd:predict}, {cmd:estat}, {cmd:stcurve}, {cmd:margins}, {cmd:test},
{cmd:lrtest}, and the other postestimation commands of {cmd:stcox} work after
{cmd:stcoxem}; see {help stcoxem postestimation}.


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{opth shared(varname)} is required and specifies the variable that defines
the frailty groups: records with the same value of {it:varname} share one
frailty. {it:varname} may be numeric or string. Typically it is the subject
identifier (one frailty per subject, recurrent events) or a cluster such as a
family or a hospital. Every record of a subject (as set by {cmd:stset}'s
{opt id()}) should belong to the same group; otherwise a warning is shown, as
by {cmd:stcox}.

{phang}
{opth strata(varname)} specifies a stratified model: a separate baseline
hazard is estimated for each level of {it:varname}, while b, theta, and the
frailties are common. {it:varname} may be numeric or string. Groups may span
strata. ({cmd:stcox} does not allow {opt strata()} with {opt shared()}.)

{phang}
{opt forceshared} allows the model to be fit when the data contain delayed
entry or gaps, which are otherwise refused (error 498), as by {cmd:stcox}; see
{help st_forceshared}. The results are then consistent only under the
assumption that the frailty distribution is independent of the covariates and
the truncation points. Multiple records per subject without gaps (for example,
after {cmd:stsplit}) are not delayed entry.

{dlgtab:SE/Robust}

{phang}
{opt vce(vcetype)} specifies how the standard errors of the coefficients are
computed.

{phang2}
{cmd:vce(oim)}, the default, uses the observed information of b and theta
together, obtained from the EM algorithm by the method of
{help stcoxem##louis1982:Louis (1982)}. The standard errors of the
coefficients agree closely with those of {cmd:stcox, shared()} (within a few
percent).

{phang2}
{cmd:vce(robust)} uses the robust (sandwich) estimator of
{help stcoxem##linwei1989:Lin and Wei (1989)} for the coefficients,
clustered on the {opt shared()} groups, at the fitted frailties, with the
small-sample factor G/(G-1). It equals the cluster-robust variance of
{cmd:stcox, vce(cluster} {it:shared_var}{cmd:)} with ln(u) as an offset.

{phang2}
{cmd:vce(cluster} {it:clustvar}{cmd:)} is the same with clusters defined by
{it:clustvar}. Every {opt shared()} group must lie within one cluster.
{cmd:vce(robust)} is {cmd:vce(cluster} {it:shared_var}{cmd:)}.

{phang2}
{cmd:vce(bootstrap)} resamples whole {opt shared()} groups with replacement
({opt reps()} times; each copy of a group is a new group) and refits the
model; the standard errors of b and theta are the standard deviations of the
replicates.

{phang2}
{cmd:vce(none)} computes no standard errors: no {cmd:e(V)}, no Wald test, no
standard error of theta. This is the fastest fit and is meant for runs that
need only the point estimates (b, theta, {cmd:e(H0mean)}).

{pmore}
Under {cmd:vce(robust)} and {cmd:vce(cluster)}, theta's standard error is
model based (see {opt noplse}); the model-based variance of b is stored in
{cmd:e(V_modelbased)}.

{phang}
{opt reps(#)} specifies the number of bootstrap replications; the default is
{cmd:reps(50)}, the minimum 2. Replications that do not converge are dropped
and counted in {cmd:e(boot_reps_maxiter)}.

{phang}
{opt bsseed(#)} sets the Mata random-number seed for {cmd:vce(bootstrap)}, so
that the replications can be reproduced.

{phang}
{opt bswarmstart} starts each bootstrap replication from the estimates
instead of from b=0, theta=0.5. It is faster, but when theta is near 0 a
replication may stop at a spurious solution with a large theta, inflating
the standard error of theta; the default cold start avoids this.

{phang}
{opt noplse} specifies that the standard error of theta be the
observed-information (Louis) value. By default, under {cmd:vce(oim)},
{cmd:vce(robust)}, and {cmd:vce(cluster)}, the standard error of theta is
obtained from the curvature of the profile log likelihood of theta, which
requires five additional fits with theta held fixed. The observed-information
standard error of theta is too small, typically by 10% to 20%, because it
treats the baseline hazard as known (see
{help stcoxem##methods:Methods and formulas}); the profile-likelihood value
agrees with {cmd:stcox, shared()}. Both are stored, in {cmd:e(se_theta)} and
{cmd:e(se_theta_oim)}; {cmd:e(se_theta_source)} records which one
{cmd:e(se_theta)} holds. The standard errors of b are not affected. Under
{cmd:vce(bootstrap)} the standard error of theta is the bootstrap one, and
under {cmd:vce(none)} there is none.

{phang}
{opt plci} adds a profile-likelihood confidence interval for theta, at the
{opt level()} of the fit, found by refitting the model at fixed values of
theta and inverting the likelihood-ratio test. It is shown in the theta row
of the table and stored in {cmd:e(theta_pl_lower)} and
{cmd:e(theta_pl_upper)}. It requires many refits.

{phang}
{opt nolrtest} suppresses the likelihood-ratio test of theta=0, which
otherwise requires fitting the Cox model without frailty.

{dlgtab:Reporting}

{phang}
{opt level(#)}; see {helpb estimation options##level():[R] Estimation options}.

{phang}
{opt nohr} reports the coefficients b rather than the hazard ratios exp(b).
This option affects only how results are displayed. Theta is never
exponentiated.

{phang}
{opt noheader} suppresses the header above the coefficient table.

{phang}
{opt nocoef} suppresses the coefficient table.

{marker display_options}{...}
{phang}
{it:display_options}:
{opt noci},
{opt nopv:alues},
{opt noomit:ted},
{opt vsquish},
{opt noempty:cells},
{opt base:levels},
{opt allbase:levels},
{opt nofvlab:el},
{opt fvwrap(#)},
{opt fvwrapon(style)},
{opth cformat(%fmt)},
{opt pformat(%fmt)},
{opt sformat(%fmt)}, and
{opt nolstretch};
see {helpb estimation options##display_options:[R] Estimation options}.

{dlgtab:Maximization}

{phang}
{opt maxiter(#)} sets the maximum number of EM iterations; the default is
the current value of {cmd:set maxiter} ({cmd:c(maxiter)}, 300 unless reset),
as for {cmd:stcox}'s {opt iterate()}; see {helpb set maxiter}. If it is
reached, {cmd:e(converged)} is 0 and a warning is shown. EM converges
linearly, so slowly converging data (for example, delayed entry with a small
theta) may need a larger value.

{phang}
{opt tolerance(#)} and {opt ltolerance(#)} set the convergence criteria. The
EM iterations stop when the largest relative change in the parameters (b and
theta), |new-old|/(|old|+1), is below {opt tolerance()} and the relative
change in the log likelihood, |new-old|/(|old|+1), is below
{opt ltolerance()}. The defaults are {cmd:tolerance(1e-8)} and
{cmd:ltolerance(1e-7)}. The parameter tolerance is tighter than {cmd:stcox}'s
because EM converges linearly: with a looser tolerance it stops measurably
short of the maximum.

{phang}
{opt coxiter(#)} sets the number of Newton steps for b within each EM
iteration; the default is {cmd:coxiter(1)}. The solution does not depend on
{opt coxiter()}; more steps may reduce the number of EM iterations.

{phang}
{opt accelerate(method)} specifies the EM acceleration:
{opt squarem} ({help stcoxem##varadhan2008:Varadhan and Roland 2008}; the
default), {opt aitken} (Aitken's delta-squared extrapolation), or {opt none}
(plain EM). Acceleration changes the number of iterations, not the solution.

{phang}
{opt noaitken} is a synonym for {cmd:accelerate(none)}.

{phang}
{opt nolog} suppresses the iteration log, which shows the log likelihood and
theta at each EM iteration.

{phang}
{opt dots} displays one dot per EM iteration instead of the iteration log.
{opt dots} may not be combined with {opt nolog}.

{dlgtab:Generated variables}

{phang}
{opth generate(newvar)} creates {it:newvar} containing, for each record in
the estimation sample, the increase of the baseline cumulative hazard over
the record's interval, H0(_t) - H0(_t0), at x=0 and u=1. Summed over a
subject's records it is the subject's baseline cumulative hazard; see also
{cmd:predict, basechazard}.

{phang}
{opth ufrailty(newvar)} creates {it:newvar} containing the estimated frailty
of each record's group, the posterior mean E(u | data). Its logarithm is
what {cmd:predict, effects} returns after {cmd:stcoxem}.


{marker weights}{...}
{title:Weights}

{pstd}
{cmd:pweight}s are specified when the data are {cmd:stset},

{phang2}{cmd:. stset time [pweight=wvar], failure(died) id(id)}{p_end}

{pstd}
and {cmd:stcoxem} uses them; other weight types are not allowed.

{pstd}
A weight k is applied to the group: a group whose records have weight k
contributes as k independent copies of the group, each with its own frailty.
Thus integer weights give the same estimates as the data expanded k times
with a new group identifier for each copy. If the weight varies within a
group, the group's mean weight is used for the frailty part of the
likelihood. ({cmd:stcox} does not allow {cmd:pweight}s with
{opt shared()}.)

{pstd}
Under {cmd:vce(oim)}, the standard errors are those of the expanded data
(copies treated as independent groups). For sampling weights, specify
{cmd:vce(robust)} or {cmd:vce(cluster)}, which give the usual sampling-weight
sandwich estimator.


{marker remarks}{...}
{title:Remarks}

{pstd}
Remarks are presented under the following headings:

{phang2}{help stcoxem##rem_stcox:Relation to stcox, shared()}{p_end}
{phang2}{help stcoxem##rem_speed:Large data}{p_end}
{phang2}{help stcoxem##rem_theta:Standard error of theta}{p_end}
{phang2}{help stcoxem##rem_ll:Log likelihood and tests}{p_end}
{phang2}{help stcoxem##rem_h0:Baseline hazard and e(H0mean)}{p_end}
{phang2}{help stcoxem##rem_other:Other remarks}{p_end}

{marker rem_stcox}{...}
{title:Relation to stcox, shared()}

{pstd}
{cmd:stcoxem} and {cmd:stcox, shared(}{it:id}{cmd:)} (which uses Breslow ties
by default) estimate the same model and give the same estimates up to the
convergence tolerance: b and theta to 5 or 6 significant digits, the
frailties and the baseline cumulative hazard to about 1e-5, and the log
likelihood {cmd:e(ll)}. Ties are always handled by the Breslow method.

{pstd}
The results are stored as {cmd:stcox, shared()} stores them ({cmd:e(cmd)} is
{cmd:stcox_fr}, {cmd:e(cmd2)} is {cmd:stcox}). {cmd:e(estimator)} is
{cmd:em} and {cmd:e(cmdline)} begins with {cmd:stcoxem}, which identifies
{cmd:stcoxem}'s results. Typing {cmd:stcoxem} redisplays them in
{cmd:stcoxem}'s format; typing {cmd:stcox} redisplays them in {cmd:stcox}'s.

{marker rem_speed}{...}
{title:Large data}

{pstd}
For runs that need only the estimates, specify {cmd:vce(none)}: the
standard errors, the five profile-likelihood fits for theta, and (with
{opt nolrtest}) the comparison Cox model are then skipped.

{phang2}{cmd:. stcoxem i.drug age, shared(id) vce(none) nolrtest nolog}{p_end}

{pstd}
Postestimation statistics that refit the model (for example,
{cmd:predict, mgale} or {cmd:basechazard}) use {cmd:stcox}'s own code, which
for unweighted fits re-solves {cmd:stcox}'s penalized likelihood with one
parameter per group. With very many groups this can be slow; see
{help stcoxem postestimation}.

{marker rem_theta}{...}
{title:Standard error of theta}

{pstd}
The default standard error of theta is computed from the profile log
likelihood (see {opt noplse}) and agrees with {cmd:stcox, shared()}. The
observed-information value ({cmd:e(se_theta_oim)}, or {cmd:e(se_theta)} with
{opt noplse}) is too small. For a confidence interval for theta that does not
rely on normality, use {opt plci}; theta is bounded at 0, and its sampling
distribution is skewed when it is small.

{marker rem_ll}{...}
{title:Log likelihood and tests}

{pstd}
{cmd:e(ll)} is the marginal log likelihood on the scale of
{cmd:stcox, shared()} (with the baseline hazard profiled out, Breslow), so
{cmd:lrtest} and {cmd:estat ic} can compare {cmd:stcoxem} models with each
other and with {cmd:stcox} models fit on the same sample. The
likelihood-ratio test of theta=0 compares the model with the Cox model
without frailty ({cmd:e(ll_c)}); because theta=0 is on the boundary, the
statistic is distributed as a 50:50 mixture of chi2(0) and chi2(1)
(chibar2(01)), as in {cmd:stcox}.

{pstd}
The Wald test in the header ({cmd:e(chi2)}) tests all coefficients; it is
not computed under {cmd:vce(none)}.

{marker rem_h0}{...}
{title:Baseline hazard and e(H0mean)}

{pstd}
The baseline hazard is the Breslow estimator at x=0 and u=1, as in
{cmd:stcox}. {cmd:e(H0mean)} is the mean over subjects of the baseline
cumulative hazard accumulated over each subject's follow-up (the sum of
{opt generate()} over a subject's records), with subjects defined by
{cmd:stset}'s {opt id()} (each record is a subject if {opt id()} was not
set). The covariates are centred internally, so the fit itself is not
affected by their location, but the baseline is reported at x=0: if x=0 is
far outside the data (for example, a year of birth), the baseline can
overflow and {cmd:e(H0mean)} is missing. Centre such covariates near 0
before fitting.

{marker rem_other}{...}
{title:Other remarks}

{pstd}
Covariates that are collinear (or constant) are omitted, with a note, and
appear with the {cmd:o.} prefix in {cmd:e(b)}, as with {cmd:stcox}.

{pstd}
Multiple-record data are handled through {cmd:stset} as usual; delayed entry
and gaps require {opt forceshared}, as with {cmd:stcox}. Time-varying covariates must be created by splitting
the records (for example, with {cmd:stsplit}); {cmd:stcoxem} has no
{opt tvc()} option.

{pstd}
Starting values are b=0 and theta=0.5. When the data contain little frailty
variation, theta approaches 0 and EM converges slowly; check
{cmd:e(converged)}.


{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. webuse catheter}{p_end}
{phang2}{cmd:. stset time, failure(infect)}{p_end}

{pstd}Gamma shared frailty on patient{p_end}
{phang2}{cmd:. stcoxem age female, shared(patient)}{p_end}

{pstd}The same model by {cmd:stcox}{p_end}
{phang2}{cmd:. stcox age female, shared(patient)}{p_end}

{pstd}Coefficients instead of hazard ratios; replay{p_end}
{phang2}{cmd:. stcoxem age female, shared(patient) nolog}{p_end}
{phang2}{cmd:. stcoxem, nohr}{p_end}

{pstd}Robust standard errors and a profile-likelihood CI for theta{p_end}
{phang2}{cmd:. stcoxem age female, shared(patient) vce(robust) plci}{p_end}

{pstd}Stratified baseline hazard{p_end}
{phang2}{cmd:. stcoxem age, shared(patient) strata(female)}{p_end}

{pstd}Bootstrap standard errors{p_end}
{phang2}{cmd:. stcoxem age female, shared(patient) vce(bootstrap) reps(200) bsseed(12345)}{p_end}

{pstd}Estimates only, for large data{p_end}
{phang2}{cmd:. stcoxem age female, shared(patient) vce(none) nolrtest nolog}{p_end}

{pstd}Frailties and baseline cumulative hazard{p_end}
{phang2}{cmd:. stcoxem age female, shared(patient) ufrailty(u) generate(dH0)}{p_end}
{phang2}{cmd:. predict nu, effects}{p_end}

{pstd}Sampling weights (set with {cmd:stset}){p_end}
{phang2}{cmd:. generate w = 1 + (age > 45)}{p_end}
{phang2}{cmd:. stset time [pweight=w], failure(infect)}{p_end}
{phang2}{cmd:. stcoxem age female, shared(patient) vce(robust)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:stcoxem} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(N_sub)}}number of subjects{p_end}
{synopt:{cmd:e(N_fail)}}number of failures{p_end}
{synopt:{cmd:e(risk)}}total time at risk{p_end}
{synopt:{cmd:e(N_g)}}number of groups{p_end}
{synopt:{cmd:e(g_min)}}smallest group size (records){p_end}
{synopt:{cmd:e(g_avg)}}average group size{p_end}
{synopt:{cmd:e(g_max)}}largest group size{p_end}
{synopt:{cmd:e(N_strata)}}number of strata ({opt strata()} only){p_end}
{synopt:{cmd:e(N_clust)}}number of clusters ({cmd:vce(robust)},
{cmd:vce(cluster)}){p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom{p_end}
{synopt:{cmd:e(rank)}}rank of {cmd:e(V)}{p_end}
{synopt:{cmd:e(chi2)}}Wald chi-squared{p_end}
{synopt:{cmd:e(ll)}}log likelihood{p_end}
{synopt:{cmd:e(ll_c)}}log likelihood, Cox model without frailty{p_end}
{synopt:{cmd:e(chi2_c)}}chi-squared, LR test of theta=0{p_end}
{synopt:{cmd:e(p_c)}}p-value, LR test of theta=0{p_end}
{synopt:{cmd:e(theta)}}frailty variance theta{p_end}
{synopt:{cmd:e(se_theta)}}standard error of theta{p_end}
{synopt:{cmd:e(se_theta_oim)}}observed-information standard error of
theta{p_end}
{synopt:{cmd:e(se_theta_profile)}}profile-likelihood standard error of theta
(missing if it could not be computed){p_end}
{synopt:{cmd:e(theta_profile_nonconv)}}number of the five profile fits that
did not converge{p_end}
{synopt:{cmd:e(se_theta_naive)}}standard error of theta ignoring its
correlation with b (diagnostic){p_end}
{synopt:{cmd:e(schur_term)}}b-theta coupling term (diagnostic){p_end}
{synopt:{cmd:e(theta_pl_lower)}}lower profile-likelihood confidence limit for
theta ({opt plci}){p_end}
{synopt:{cmd:e(theta_pl_upper)}}upper profile-likelihood confidence limit for
theta ({opt plci}){p_end}
{synopt:{cmd:e(theta_pl_nonconv)}}number of {opt plci} fits that did not
converge{p_end}
{synopt:{cmd:e(H0mean)}}mean baseline cumulative hazard over subjects{p_end}
{synopt:{cmd:e(converged)}}{cmd:1} if converged, {cmd:0} otherwise{p_end}
{synopt:{cmd:e(iter)}}number of EM iterations{p_end}
{synopt:{cmd:e(tolerance)}}parameter tolerance{p_end}
{synopt:{cmd:e(ltolerance)}}log-likelihood tolerance{p_end}
{synopt:{cmd:e(mem_est)}}estimated memory use (MB; lower bound){p_end}
{synopt:{cmd:e(squarem_att)}, {cmd:e(squarem_acc)}}SQUAREM steps attempted
and accepted{p_end}
{synopt:{cmd:e(aitken_theta_att)}, {cmd:e(aitken_theta_acc)}}Aitken steps
for theta attempted and accepted{p_end}
{synopt:{cmd:e(aitken_beta_att)}, {cmd:e(aitken_beta_acc)}}Aitken steps for
b attempted and accepted{p_end}
{synopt:{cmd:e(N_reps)}}number of bootstrap replications used{p_end}
{synopt:{cmd:e(boot_reps)}}number of bootstrap replications requested{p_end}
{synopt:{cmd:e(boot_reps_used)}}number of replications used{p_end}
{synopt:{cmd:e(boot_reps_maxiter)}}number of replications that did not
converge{p_end}
{synopt:{cmd:e(boot_reps_badtheta)}}number of replications with invalid
theta{p_end}
{synopt:{cmd:e(boot_reps_badbeta)}}number of replications with invalid
b{p_end}
{synopt:{cmd:e(boot_avg_iter)}}average EM iterations per replication{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:stcox_fr}{p_end}
{synopt:{cmd:e(cmd2)}}{cmd:stcox}{p_end}
{synopt:{cmd:e(estimator)}}{cmd:em}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}{cmd:_t}{p_end}
{synopt:{cmd:e(t0)}}{cmd:_t0}{p_end}
{synopt:{cmd:e(shared)}}frailty grouping variable{p_end}
{synopt:{cmd:e(strata)}}stratification variable{p_end}
{synopt:{cmd:e(method)}}{cmd:breslow}{p_end}
{synopt:{cmd:e(ties)}}{cmd:breslow}{p_end}
{synopt:{cmd:e(wtype)}}weight type ({cmd:pweight}){p_end}
{synopt:{cmd:e(wexp)}}weight expression{p_end}
{synopt:{cmd:e(clustvar)}}name of cluster variable{p_end}
{synopt:{cmd:e(vce)}}{it:vcetype} specified in {cmd:vce()}{p_end}
{synopt:{cmd:e(vcetype)}}title used to label Std. err.{p_end}
{synopt:{cmd:e(se_theta_source)}}source of {cmd:e(se_theta)}: {cmd:plse},
{cmd:oim ...}, or {cmd:bootstrap}{p_end}
{synopt:{cmd:e(chi2type)}}{cmd:Wald}{p_end}
{synopt:{cmd:e(accelerate)}}EM acceleration method{p_end}
{synopt:{cmd:e(boot_start)}}{cmd:cold} or {cmd:warm} ({cmd:vce(bootstrap)}){p_end}
{synopt:{cmd:e(sh_warn)}}{cmd:sh_warn} if subjects span frailty groups{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}
{synopt:{cmd:e(predict)}}program used to implement {cmd:predict}:
{cmd:stcox_p}, or {cmd:stcoxem_p} with weights{p_end}
{synopt:{cmd:e(estat_cmd)}}program used to implement {cmd:estat}{p_end}
{synopt:{cmd:e(datasignature)}}the checksum{p_end}
{synopt:{cmd:e(datasignaturevars)}}variables used in calculation of
checksum{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}
{synopt:{cmd:e(V_modelbased)}}model-based variance ({cmd:vce(robust)},
{cmd:vce(cluster)}){p_end}
{synopt:{cmd:e(boot_theta_all)}}theta of each bootstrap replication{p_end}
{synopt:{cmd:e(boot_beta_all)}}b of each bootstrap replication{p_end}
{synopt:{cmd:e(boot_iter_all)}}EM iterations of each replication{p_end}
{synopt:{cmd:e(boot_outcome_all)}}outcome of each replication (1 used, 2 not
converged, 3 invalid theta, 4 invalid b){p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}


{marker methods}{...}
{title:Methods and formulas}

{pstd}
Let group i have frailty u_i ~ Gamma(1/theta, 1/theta), D_i failures, and
cumulative hazard H_i = sum over its records of
[H0(_t) - H0(_t0)] exp(xb). Given b, theta, and H0, the posterior of u_i is
Gamma(1/theta + D_i, 1/theta + H_i).

{pstd}
E-step: E(u_i) = (1/theta + D_i)/(1/theta + H_i) and E(ln u_i) are computed
in closed form.

{pstd}
M-step: b maximizes the Cox partial likelihood with offset ln E(u_i)
({opt coxiter()} Newton steps); H0 is the Breslow estimator with the same
offset, computed in one pass over the sorted risk sets per stratum; theta
solves its own score equation, ln(theta) + digamma(1/theta) - 1 +
mean[E(u_i) - E(ln u_i)] = 0, by a safeguarded Newton step. The algorithm
is that of {help stcoxem##klein1992:Klein (1992)}; the fixed point is the
maximum of the marginal likelihood with the Breslow baseline profiled out,
which is the estimate {cmd:stcox, shared()} reports.

{pstd}
Log likelihood: with the Breslow increments dH0_k = d_k/S_k at the distinct
failure times (d_k failures, S_k the risk-set sum of E(u) exp(xb)),

{p 8 8 2}
e(ll) = sum_i [ln Gamma(1/theta + D_i) - ln Gamma(1/theta) -
(1/theta) ln(theta) - (1/theta + D_i) ln(1/theta + H_i)] +
sum over failures of xb + sum_k d_k [ln(dH0_k) - ln(d_k) + 1],

{pstd}
that is, the marginal log likelihood with the baseline profiled out, minus
the data-only constant sum_k d_k ln(d_k), plus the number of failures. This
puts it on {cmd:stcox}'s scale, where the log likelihood equals the Breslow
partial log likelihood as theta -> 0. With weights, D_i, d_k, and the sums
are weighted.

{pstd}
Variance: under {cmd:vce(oim)}, the observed information of (b, theta) is
the complete-data information minus the missing information (Louis 1982),
evaluated at the EM solution. It treats the Breslow baseline as known, so
it ignores the uncertainty in H0; for b this matters little, for theta it
makes the standard error too small
({help stcoxem##nielsen1992:Nielsen et al. 1992}). The profile-likelihood
standard error of theta is 1/sqrt(-d2 pl(theta)/d theta2), where pl(theta)
is the log likelihood maximized over b and H0 at fixed theta; the second
derivative is computed by central differences at two step sizes with
Richardson extrapolation. The robust variance of b is the
{help stcoxem##linwei1989:Lin and Wei (1989)} sandwich with score residuals
at the fitted frailties, summed within clusters, times G/(G-1).

{pstd}
Covariates are centred internally; b, theta, the log likelihood, and the
variances do not change, and the baseline is transformed back to x=0.
Collinearity is checked on the centred, scaled covariates. Tied times are
ordered deterministically, so repeated fits give identical results.


{marker references}{...}
{title:References}

{marker klein1992}{...}
{phang}
Klein, J. P. 1992. Semiparametric estimation of random effects using the
Cox model based on the EM algorithm. {it:Biometrics} 48: 795-806.

{marker linwei1989}{...}
{phang}
Lin, D. Y., and L. J. Wei. 1989. The robust inference for the Cox
proportional hazards model. {it:Journal of the American Statistical
Association} 84: 1074-1078.

{marker louis1982}{...}
{phang}
Louis, T. A. 1982. Finding the observed information matrix when using the
EM algorithm. {it:Journal of the Royal Statistical Society, Series B} 44:
226-233.

{marker nielsen1992}{...}
{phang}
Nielsen, G. G., R. D. Gill, P. K. Andersen, and T. I. A. Sorensen. 1992. A
counting process approach to maximum likelihood estimation in frailty
models. {it:Scandinavian Journal of Statistics} 19: 25-43.

{marker varadhan2008}{...}
{phang}
Varadhan, R., and C. Roland. 2008. Simple and globally convergent methods
for accelerating the convergence of any EM algorithm.
{it:Scandinavian Journal of Statistics} 35: 335-353.


{marker author}{...}
{title:Author}

{pstd}
Barna Szabo-Soderberg{break}
Karolinska Institutet, Stockholm, Sweden{break}
{browse "https://github.com/szabobarnaKI/stcoxem"}


{title:Also see}

{psee}
Help: {help stcoxem postestimation}, {manhelp stcox ST},
{manhelp stcox_postestimation ST:stcox postestimation}, {manhelp stset ST}
{p_end}
