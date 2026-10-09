{smcl}
{* *! version 1.0.0  28sep2026}{...}
{viewerjumpto "Syntax" "stcoxem##syntax"}{...}
{viewerjumpto "Description" "stcoxem##description"}{...}
{viewerjumpto "Options" "stcoxem##options"}{...}
{viewerjumpto "Weights" "stcoxem##weights"}{...}
{viewerjumpto "Remarks" "stcoxem##remarks"}{...}
{viewerjumpto "Stored results" "stcoxem##results"}{...}
{viewerjumpto "Methods and formulas" "stcoxem##methods"}{...}
{viewerjumpto "References" "stcoxem##references"}{...}
{viewerjumpto "Examples" "stcoxem##examples"}{...}
{title:Title}

{phang}
{bf:stcoxem} {hline 2} Semiparametric Cox model with shared gamma frailty, fit by the
Klein (1992) EM algorithm


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:stcoxem} {it:{help varlist:indepvars}} {ifin}{cmd:,} {cmdab:sh:ared(}{it:varname}{cmd:)}
[{it:options}]

{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt:{cmdab:sh:ared(}{it:varname}{cmd:)}}variable identifying the shared-frailty
groups (clusters); {bf:required}{p_end}
{synopt:{cmdab:str:ata(}{it:varname}{cmd:)}}estimate a separate baseline hazard
within each level of {it:varname}{p_end}

{syntab:SE/CI}
{synopt:{cmdab:vce:(}{it:vcetype}{cmd:)}}{it:vcetype} may be {cmd:none} (default),
{cmd:oim}, {cmd:robust}, or {cmd:bootstrap}{p_end}
{synopt:{cmd:plci}}report a profile-likelihood confidence interval for {cmd:theta}{p_end}
{synopt:{cmdab:nopl:se}}suppress the profile-likelihood-based SE for {cmd:theta};
by default {cmd:plse} runs under {cmd:vce(oim)}/{cmd:vce(robust)} only{p_end}
{synopt:{cmdab:rep:s(}{it:#}{cmd:)}}number of bootstrap replications; default is
{cmd:reps(50)}; requires {cmd:vce(bootstrap)}{p_end}
{synopt:{cmd:bsseed(}{it:#}{cmd:)}}Mata random-number seed for
{cmd:vce(bootstrap)}{p_end}
{synopt:{cmd:bswarmstart}}warm-start each bootstrap replicate from the main fit
instead of cold-starting{p_end}
{synopt:{cmdab:l:evel(}{it:#}{cmd:)}}set confidence level; default is
{cmd:level(95)}{p_end}

{syntab:Reporting}
{synopt:{cmdab:noh:r}}report coefficients instead of hazard ratios{p_end}
{synopt:{cmdab:nol:og}}suppress the per-iteration EM log{p_end}
{synopt:{cmd:dots}}show a progress dot per EM iteration instead of the full log{p_end}

{syntab:EM algorithm}
{synopt:{cmdab:maxi:ter(}{it:#}{cmd:)}}maximum EM iterations; default is
{cmd:maxiter(200)}{p_end}
{synopt:{cmdab:tol:erance(}{it:#}{cmd:)}}EM convergence tolerance on relative
parameter change; default is {cmd:tolerance(1e-6)}{p_end}
{synopt:{cmdab:coxi:ter(}{it:#}{cmd:)}}Newton steps for beta per EM cycle;
default is {cmd:coxiter(1)}{p_end}
{synopt:{cmdab:accel:erate(}{it:method}{cmd:)}}EM acceleration: {cmd:squarem}
(default), {cmd:aitken}, or {cmd:none}{p_end}
{synopt:{cmdab:noai:tken}}synonym for {cmd:accelerate(none)}{p_end}

{syntab:Output variables}
{synopt:{cmdab:gen:erate(}{it:newvar}{cmd:)}}record-level Breslow baseline
cumulative hazard{p_end}
{synopt:{cmd:ufrailty(}{it:newvar}{cmd:)}}record-level posterior mean frailty
E[u|data] for the record's {cmd:shared()} cluster{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
{it:indepvars} may contain factor variables; see {help fvvarlist}.{p_end}
{p 4 6 2}
{it:indepvars} may be omitted, fitting a frailty-only model with no covariates.{p_end}
{p 4 6 2}
{cmd:pweight}s are supported but are set on {cmd:stset}, not on {cmd:stcoxem};
see {help stcoxem##weights:Weights} below.{p_end}
{p 4 6 2}
The data must already be {cmd:stset}; see {help stset}.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:stcoxem} fits a semiparametric Cox proportional hazards model with a
shared multiplicative gamma frailty, using the EM algorithm of
{help stcoxem##KleinRef:Klein (1992)}. Subjects (or higher-level clusters,
identified by {cmd:shared()}) share an unobserved frailty {cmd:u} that
multiplies their hazard; {cmd:u} is assumed Gamma-distributed with mean 1 and
unknown variance {cmd:theta}, which {cmd:stcoxem} estimates jointly with the
regression coefficients.

{pstd}
The model is

{pstd}
{cmd:h(t | x, u) = u * h0(t) * exp(x'b)},  {cmd:u ~ Gamma(mean 1, variance theta)},
shared within each {cmd:shared()} group,

{pstd}
where {cmd:h0(t)} is an unspecified (nonparametric) baseline hazard, estimated
as a Breslow step function exactly as {helpb stcox} does.

{pstd}
{cmd:stcoxem} was written as a faster, more stable alternative to
{cmd:stcox, shared() frailty(gamma)} for very large registry
data: it implements the E- and M-steps as O(n) risk-set sweeps in Mata, does
one Newton step for the regression coefficients per EM cycle by default (see
{opt coxiter()}). The primary reason for developing {cmd:stcoxem} was to 
obtain point estimates for theta and the mean of the individual cumulative 
baseline hazards faster, therefore as default it computes no standard errors 
at all unless requested (see {opt vce()}).


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{opt shared(varname)} identifies the shared-frailty groups: records sharing
the same value of {it:varname} share one common (unobserved) frailty draw.
{it:varname} may be numeric or string. {cmd:shared()} is required.

{phang}
{opt strata(varname)} estimates a separate nonparametric baseline hazard
within each level of {it:varname} (as {helpb stcox}'s {opt strata()} does).
{it:varname} may be numeric or string.

{dlgtab:SE/CI}

{phang}
{opt vce(vcetype)} selects how (or whether) standard errors are computed.

{phang2}
{cmd:vce(none)}, the default, computes point estimates only. No standard
errors, z statistics, p-values, or confidence intervals are available, and
{cmd:e(V)} is not returned. This is the fastest option and is recommended
for large-scale production runs (see {help stcoxem##remarks:Remarks}). Under
{cmd:vce(none)}, {opt plse} does not run by default either, since there is
no {cmd:e(se_theta)} for it to feed.

{phang2}
{cmd:vce(oim)} reports model-based standard errors for both the regression
coefficients and {cmd:theta}, from a joint Louis (1982) EM-information
correction. Beta's SE matches {cmd:stcox}'s closely; {cmd:theta}'s analytic
SE is a known, understood ~13-20% underestimate relative to {cmd:stcox},
bootstrap, or R's {cmd:frailtyEM} (a semiparametric-information deficit --
see {help stcoxem##remarks:Remarks} and
{help stcoxem##methods:Methods and formulas}), which is why {opt plse} runs
by default alongside {cmd:vce(oim)} and is what {cmd:e(se_theta)} reports.

{phang2}
{cmd:vce(robust)} reports a cluster-robust (Lin & Wei 1989) sandwich SE for
the regression coefficients, clustered on {cmd:shared()}; {cmd:theta}'s SE
is the same joint Louis-correction value as {cmd:vce(oim)}, and {opt plse}
runs by default here too, feeding {cmd:e(se_theta)} exactly as it does under
{cmd:vce(oim)}.

{phang2}
{cmd:vce(bootstrap)} reports nonparametric cluster-bootstrap standard errors
(resampling whole {cmd:shared()} clusters with replacement, each resampled
copy assigned a fresh cluster id) for both the coefficients and {cmd:theta}.
See {opt reps()}, {opt bsseed()}, and {opt bswarmstart}.

{phang}
{opt plci} additionally reports a profile-likelihood confidence interval for
{cmd:theta}, obtained by refitting the model at a grid of fixed {cmd:theta}
values and inverting the likelihood-ratio test. This is independent of
{opt vce()} and does not rely on the Louis decomposition. It is
substantially more expensive than a single fit, since each bound search is
several full EM refits.

{phang}
{opt noplse} suppresses the profile-likelihood-based standard error for
{cmd:theta} ({cmd:plse}). By default, {cmd:plse} runs automatically under
{cmd:vce(oim)} and {cmd:vce(robust)} -- where it is what {cmd:e(se_theta)}
reports -- but not under {cmd:vce(none)} or {cmd:vce(bootstrap)}, where there
is no {cmd:e(se_theta)} for it to feed and running it would only cost time
for no displayed benefit. {cmd:plse} evaluates the profile log-likelihood in
{cmd:theta} at five points (two step sizes, Richardson-extrapolated) and is
noticeably more expensive than a plain fit; specify {opt noplse} under
{cmd:vce(oim)}/{cmd:vce(robust)} for maximum speed when {cmd:theta}'s
profile-likelihood SE is not needed (the analytic Louis/EM SE remains
available in {cmd:e(se_theta_oim)} either way). See
{help stcoxem##remarks:Remarks}.

{phang}
{opt reps(#)} sets the number of bootstrap replications for
{cmd:vce(bootstrap)}; the default is {cmd:reps(50)}, and {it:#} must be at
least 2.

{phang}
{opt bsseed(#)} sets the Mata random-number seed (via {cmd:rseed()}) used to
draw the bootstrap resamples, for reproducibility. If omitted, the ambient
Mata random-number state is used.

{phang}
{opt bswarmstart} warm-starts every bootstrap replicate's EM iterations from
the main fit's converged coefficients and {cmd:theta}, instead of the
default cold start ({cmd:beta=0}, {cmd:theta=0.5}). Warm-starting is faster
but can bias {cmd:theta}'s bootstrap SE upward if a resample's true
{cmd:theta} is near zero (the warm start can pull a resample toward a
spurious high-{cmd:theta} fixed point); compare against the cold-started
default if results look implausible.

{phang}
{opt level(#)}; see {help level:{it:[R] level}}.

{dlgtab:Reporting}

{phang}
{opt nohr} reports the estimated coefficients on the linear-predictor scale
instead of the default hazard-ratio (exponentiated) scale. {cmd:theta} is
never exponentiated, regardless of {opt nohr}.

{phang}
{opt nolog} suppresses the per-iteration EM log (the default display shows
one line per EM iteration with the log-likelihood and current {cmd:theta}).

{phang}
{opt dots} shows a single progress dot per EM iteration instead of the full
per-iteration log. {opt nolog} and {opt dots} may not be combined.

{dlgtab:EM algorithm}

{phang}
{opt maxiter(#)} sets the maximum number of EM iterations; the default is
{cmd:maxiter(200)}. {cmd:e(converged)} reports whether the algorithm
actually converged (relative parameter change below {opt tolerance()})
before hitting this limit.

{phang}
{opt tolerance(#)} sets the EM convergence criterion: iteration stops once
the largest relative change in {cmd:theta} and every regression coefficient,
between consecutive EM iterations, falls below {it:#}. The default is
{cmd:tolerance(1e-6)}.

{phang}
{opt coxiter(#)} sets the number of Newton steps taken for the regression
coefficients within each EM cycle's M-step; the default, {cmd:coxiter(1)},
takes a single step. Increasing {opt coxiter()} does not change the fitted
model -- {cmd:stcoxem} converges to the same fixed point regardless of
{opt coxiter()}, confirmed to floating-point precision when compared at
matching {opt tolerance()} -- but can reduce the number of EM iterations
needed to reach a very tight {opt tolerance()}.

{phang}
{opt accelerate(method)} selects the EM acceleration scheme. {cmd:squarem}
(Varadhan and Roland 2008), the default, and {cmd:aitken} (Aitken's
Delta-squared) both extrapolate past several raw EM steps to their implied
fixed point, guarded against numerically unstable jumps; {cmd:none} runs
plain, unaccelerated EM. Acceleration changes only how many iterations are
needed, never the converged answer.

{phang}
{opt noaitken} is a synonym for {opt accelerate(none)}, retained for
backward compatibility.

{dlgtab:Output variables}

{phang}
{opt generate(newvar)} creates {it:newvar} holding each record's Breslow
baseline cumulative hazard, {cmd:H0(t) - H0(t0)}, on the estimation sample.

{phang}
{opt ufrailty(newvar)} creates {it:newvar} holding each record's posterior
mean frailty, {cmd:E[u | data]}, for the record's {cmd:shared()} cluster
(constant within cluster). This is comparable to {cmd:stcox}'s
{cmd:predict, effects}.


{marker weights}{...}
{title:Weights}

{pstd}
{cmd:stcoxem} supports {help pweight}s, but -- unlike most Stata
estimation commands -- the weight is specified when the data are
{cmd:stset}, not in {cmd:stcoxem}'s own syntax:

{phang2}{cmd:. stset timevar [pweight=wtvar], ...}{p_end}
{phang2}{cmd:. stcoxem ..., shared(id) ...}{p_end}

{pstd}
{cmd:stcoxem} reads the weight {help stset} recorded and reports which
weight (if any) it used at the top of the output. Only {cmd:pweight} is
supported; {cmd:stset} with any other weight type raises an error.

{pstd}
A {cmd:pweight} is interpreted at the {cmd:shared()}-cluster level: a
cluster whose records carry a constant weight {cmd:k} is treated as {cmd:k}
independent replicate copies of that cluster for the frailty likelihood
(matching what you would get by literally replicating the cluster {cmd:k}
times, each copy with its own frailty draw). A weight that varies within a
cluster is handled via the cluster's mean weight.


{marker remarks}{...}
{title:Remarks}

{pstd}
{opt plse} runs automatically under {cmd:vce(oim)} and {cmd:vce(robust)},
because it is what {cmd:e(se_theta)} reports under those settings, and is
skipped by default under {cmd:vce(none)} and {cmd:vce(bootstrap)}, where
there is no {cmd:e(se_theta)} for it to feed. If you are running
{cmd:vce(oim)}/{cmd:vce(robust)} but do not need {cmd:theta}'s
profile-likelihood SE, specify {opt noplse} to skip this extra computation
and get a faster fit (the analytic Louis/EM SE remains available in
{cmd:e(se_theta_oim)}):

{phang2}{cmd:. stcoxem ..., shared(id) vce(oim) noplse}{p_end}

{pstd}
For large datasets and with the only goal of obtaining point estimates, 
{cmd:vce(none)} provides a fast fit with no SE of any kind computed and 
{opt plse} already does not run, so no extra option is needed for the fastest 
possible fit:

{phang2}{cmd:. stcoxem ..., shared(id) vce(none)}{p_end}

{pstd}
For very large numbers of records or clusters, fit progress can be tracked
with {opt dots} without the volume of a full per-iteration log, or
suppressed entirely with {opt nolog}.

{pstd}
{cmd:stcoxem} does not implement a {helpb predict} command. Per-record
quantities are instead obtained directly through {opt generate()} (baseline
cumulative hazard) and {opt ufrailty()} (posterior mean frailty).

{pstd}
Time-varying covariates and counting-process (multiple-record) data are
supported directly through {cmd:stset}'s usual enter()/exit()/id() syntax;
{cmd:stcoxem} does not split records for you.


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:stcoxem} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N_orig)}}number of records used{p_end}
{synopt:{cmd:e(N_clust)}}number of {cmd:shared()} clusters{p_end}
{synopt:{cmd:e(N_subj)}}number of subjects (from {cmd:stset}'s {cmd:id()}, or
equal to {cmd:e(N_orig)} if none was set){p_end}
{synopt:{cmd:e(n_strata)}}number of strata{p_end}
{synopt:{cmd:e(theta)}}estimated frailty variance{p_end}
{synopt:{cmd:e(se_theta)}}SE of {cmd:theta} (from {cmd:plse} if it ran and
converged reliably, else the analytic Louis/EM value); {cmd:vce(oim)} or
{cmd:vce(robust)} only{p_end}
{synopt:{cmd:e(se_theta_oim)}}the analytic Louis/EM SE of {cmd:theta},
always available alongside {cmd:e(se_theta)}; {cmd:vce(oim)}/{cmd:vce(robust)}
only{p_end}
{synopt:{cmd:e(se_theta_profile)}}profile-likelihood-based SE of
{cmd:theta}, missing if it did not compute reliably; present whenever
{opt plse} ran, which by default is under {cmd:vce(oim)}/{cmd:vce(robust)}
(unless {opt noplse} was given){p_end}
{synopt:{cmd:e(theta_profile_nonconv)}}number of the 5 internal profile
refits (out of 5) that did not converge within {opt maxiter()}; present
whenever {opt plse} ran{p_end}
{synopt:{cmd:e(se_theta_naive)}}{cmd:theta}'s SE ignoring beta/theta
coupling (diagnostic); {cmd:vce(oim)}/{cmd:vce(robust)} with covariates only{p_end}
{synopt:{cmd:e(schur_term)}}the beta/theta Schur-complement coupling term
(diagnostic); {cmd:vce(oim)}/{cmd:vce(robust)} with covariates only{p_end}
{synopt:{cmd:e(ll)}}Cox partial log-likelihood plus frailty terms, correct
up to an additive tie-pattern constant -- valid for comparisons within the
same model (e.g. across {opt plci()} refits) but not for cross-model
LR/AIC comparisons{p_end}
{synopt:{cmd:e(H0mean)}}mean cumulative baseline hazard across subjects{p_end}
{synopt:{cmd:e(iter)}}number of EM iterations run{p_end}
{synopt:{cmd:e(converged)}}1 if the EM algorithm converged, 0 if it stopped
at {opt maxiter()}{p_end}
{synopt:{cmd:e(mem_est)}}a lower-bound estimate (MB) of the fit's memory
footprint; not a measurement of true peak memory{p_end}
{synopt:{cmd:e(theta_pl_lower)}, {cmd:e(theta_pl_upper)}}profile-likelihood
CI bounds for {cmd:theta}; {opt plci} only{p_end}
{synopt:{cmd:e(theta_pl_nonconv)}}number of profile refits during the
{opt plci} bound search that did not converge{p_end}
{synopt:{cmd:e(aitken_theta_att)}, {cmd:e(aitken_theta_acc)}}Aitken
acceleration theta-jump attempted/accepted counts; {opt accelerate(aitken)}
only{p_end}
{synopt:{cmd:e(aitken_beta_att)}, {cmd:e(aitken_beta_acc)}}Aitken
acceleration beta-jump attempted/accepted counts; {opt accelerate(aitken)}
with covariates only{p_end}
{synopt:{cmd:e(squarem_att)}, {cmd:e(squarem_acc)}}SQUAREM acceleration
attempted/accepted counts; default {opt accelerate(squarem)} only{p_end}
{synopt:{cmd:e(boot_reps)}}requested number of bootstrap replications{p_end}
{synopt:{cmd:e(boot_reps_used)}}number of replications that converged to a
usable estimate{p_end}
{synopt:{cmd:e(boot_reps_maxiter)}}number of replications that hit
{opt maxiter()} without converging{p_end}
{synopt:{cmd:e(boot_reps_badtheta)}, {cmd:e(boot_reps_badbeta)}}number of
replications with a non-finite {cmd:theta} or coefficient{p_end}
{synopt:{cmd:e(boot_avg_iter)}}average EM iterations per bootstrap
replicate{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:stcoxem}{p_end}
{synopt:{cmd:e(depvar)}}{cmd:_t}{p_end}
{synopt:{cmd:e(title)}}title of the estimation command{p_end}
{synopt:{cmd:e(sharedvar)}}variable given to {opt shared()}{p_end}
{synopt:{cmd:e(stratavar)}}variable given to {opt strata()}, if any{p_end}
{synopt:{cmd:e(wexp)}}the weight variable read from {cmd:stset}, if any{p_end}
{synopt:{cmd:e(vce)}}{it:vcetype} given to {opt vce()}{p_end}
{synopt:{cmd:e(se_theta_source)}}{cmd:plse} or an {cmd:oim} fallback note,
describing what {cmd:e(se_theta)} reports; {cmd:vce(oim)}/{cmd:vce(robust)}
only{p_end}
{synopt:{cmd:e(accelerate)}}acceleration method used{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the coefficients;
{opt vce()} other than {cmd:none} only{p_end}
{synopt:{cmd:e(boot_theta_all)}, {cmd:e(boot_beta_all)}}per-replicate
{cmd:theta} and coefficient estimates; {cmd:vce(bootstrap)} only{p_end}
{synopt:{cmd:e(boot_iter_all)}, {cmd:e(boot_outcome_all)}}per-replicate
iteration counts and outcome codes (1=usable, 2=hit maxiter, 3=invalid
theta, 4=invalid beta); {cmd:vce(bootstrap)} only{p_end}
{p2colreset}{...}


{marker methods}{...}
{title:Methods and formulas}

{pstd}
{cmd:stcoxem} implements the EM algorithm of
{help stcoxem##KleinRef:Klein (1992)} for a shared gamma-frailty Cox model.
Given current parameter values, the E-step uses the closed-form
gamma-Poisson posterior moments of each cluster's frailty, {cmd:E[u|data]}
and {cmd:E[log u|data]}, computed from each cluster's cumulative hazard at
its own covariates. The M-step re-maximizes an offset Cox partial
likelihood, with the log posterior-mean frailty entered as an offset: the
regression coefficients by Newton's method ({opt coxiter()} steps per
cycle), and the Breslow baseline hazard in one O(n) risk-set sweep per
stratum; {cmd:theta} is then updated by a stable Newton step on its own
profile score. EM acceleration ({opt accelerate()}) extrapolates past
several raw iterations to their implied fixed point.

{pstd}
Standard errors under {opt vce(oim)}/{opt vce(robust)} use a joint Louis
(1982) EM-information decomposition for the regression coefficients and
{cmd:theta} together (the complete-data beta/theta cross-information is
exactly zero for this model, so only the missing-information term couples
them). {cmd:theta}'s own missing-information term is exact; the coupling
with the regression coefficients uses a first-order (delta-method)
linearization. The resulting analytic SE for {cmd:theta} is a known,
understood underestimate relative to fully efficient semiparametric
estimators (the decomposition treats the nonparametric Breslow baseline as
known -- see
{help stcoxem##NielsenRef:Nielsen, Gill, Andersen, and Sørensen (1992)}),
which is why the profile-likelihood SE ({opt plse}, which runs by default
alongside {opt vce(oim)}/{opt vce(robust)}) or {opt plci}/{cmd:vce(bootstrap)}
are the recommended sources for {cmd:theta}'s uncertainty. {opt vce(robust)}
additionally reports a
cluster-robust (Lin and Wei 1989) sandwich SE for the regression
coefficients.


{marker references}{...}
{title:References}

{marker KleinRef}{...}
{phang}
Klein, J. P. 1992. Semiparametric estimation of random effects using the
Cox model based on the EM algorithm. {it:Biometrics} 48: 795-806.

{phang}
Lin, D. Y., and L. J. Wei. 1989. The robust inference for the Cox
proportional hazards model. {it:Journal of the American Statistical
Association} 84: 1074-1078.

{phang}
Louis, T. A. 1982. Finding the observed information matrix when using the
EM algorithm. {it:Journal of the Royal Statistical Society, Series B} 44:
226-233.

{marker NielsenRef}{...}
{phang}
Nielsen, G. G., R. D. Gill, P. K. Andersen, and T. I. A. Sørensen. 1992. A
counting process approach to maximum likelihood estimation in frailty
models. {it:Scandinavian Journal of Statistics} 19: 25-43.

{phang}
Varadhan, R., and C. Roland. 2008. Simple and globally convergent methods
for accelerating expectation-maximization algorithms. {it:Scandinavian
Journal of Statistics} 35: 335-353.


{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. stset followup_years, enter(time 0) exit(time .) failure(outcome==1) id(patientid)}{p_end}

{pstd}Fit a shared-frailty model, point estimates only (fastest){p_end}
{phang2}{cmd:. stcoxem i.drug i.gender i.dm age_entry, shared(patientid) vce(none)}{p_end}

{pstd}With model-based standard errors{p_end}
{phang2}{cmd:. stcoxem i.drug i.gender i.dm age_entry, shared(patientid) vce(oim)}{p_end}

{pstd}Stratified by a third variable, with a profile-likelihood CI for theta{p_end}
{phang2}{cmd:. stcoxem i.drug i.gender age_entry, shared(patientid) strata(dm) plci}{p_end}

{pstd}With a pweight (set on {cmd:stset}, not on {cmd:stcoxem}){p_end}
{phang2}{cmd:. stset followup_years [pweight=samplewt], enter(time 0) exit(time .) failure(outcome==1) id(patientid)}{p_end}
{phang2}{cmd:. stcoxem i.drug i.gender i.dm age_entry, shared(patientid)}{p_end}

{pstd}Cluster-bootstrap standard errors, 200 replications{p_end}
{phang2}{cmd:. stcoxem i.drug i.gender i.dm age_entry, shared(patientid) vce(bootstrap) reps(200) bsseed(12345)}{p_end}

{pstd}Save per-record baseline cumulative hazard and posterior mean frailty{p_end}
{phang2}{cmd:. stcoxem i.drug i.gender i.dm age_entry, shared(patientid) generate(H0) ufrailty(uhat)}{p_end}
