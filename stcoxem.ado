*! version 1.0.0  26nov2026  Barna Szabo-Soderberg (rev 55)
*! stcoxem.ado — Semiparametric Cox model with shared gamma frailty,
*! fit by the Klein (1992) EM algorithm. Stata 19.0 / Mata.
*!
*! REV 55 (2026-11-26): accelerate(aitken) fixed. An Aitken jump is taken only
*!   when the ratio of the last two steps, r = d1/d2, satisfies |r| < 1 (the
*!   sequence contracts; Aitken's formula assumes it). From the cold start
*!   theta accelerates (0.500, 0.572, 0.661: r = 1.24) and the old code
*!   extrapolated to the wrong side (theta 0.0065), where EM barely moves:
*!   5000 iterations, theta 78% off on the registry test set (cert 100905).
*!   Applied to theta and to each beta component. Rev 54 probe code removed.
*! REV 54 (2026-11-26): maxiter() default is c(maxiter) (set maxiter; 300 by
*!   default), as stcox.
*! REV 53 (2026-11-26): delayed entry or gaps with shared() -> error 498 unless
*!   the (undocumented in stcox) option forceshared is given, exactly as stcox:
*!   stcox's own _st_err_sharedgaps is called after the sample is marked.
*! REV 52 (2026-11-26): version 1.0.0; help files (stcoxem.sthlp,
*!   stcoxem_postestimation.sthlp), package files, certification script.
*! - Model without covariates: the table shows only the theta row (rev <=51
*!   printed "(no standard errors: vce(none))" and the note on standard
*!   errors of hazard ratios).
*! - Note when e(H0mean) is missing (baseline at x=0 out of range, e.g. a
*!   covariate far from 0).
*! REV 51 (2026-11-26): postestimation (stcox's own commands).
*! - e(cmd)="stcox_fr", e(cmd2)="stcox", e(estat_cmd)="stcox_estat" (as
*!   stcox, shared()); e(estimator)="em" identifies stcoxem. Unweighted
*!   fits: e(predict)="stcox_p" -- predict, estat, stcurve, margins run
*!   stcox's code unchanged (Path A; refit statistics re-solve stcox's
*!   penalized fit at e(theta), which lands on the EM fixed point).
*!   Weighted fits: e(predict)="stcoxem_p" (stcoxem_p.ado): effects and
*!   atfrailty from stcoxem's own frailties; refit statistics through
*!   stcox_p on a plain-Cox view of the results (offset ln u, iter(0) at
*!   e(b)) -- stcox's weighted shared-frailty refit is a different
*!   estimator (Path B). ln u is rebuilt at e(b), e(theta) by the hidden
*!   subcommand -stcoxem _lnu newvar- (Mata stcoxem_lnu).
*! - Typing -stcox- now replays the results with stcox_fr's display;
*!   -stcoxem- replays with stcoxem's.
*! REV 50 (2026-11-26): e(ll) on stcox's scale; LR test of theta=0.
*! - e(ll) = previous value + sum(wt*d) (weighted number of failures). The
*!   old value was the marginal log likelihood minus sum_k d_k ln d_k
*!   (Breslow); stcox, shared() reports a scale that equals the Breslow
*!   partial likelihood as theta->0, which differs from it by exactly D.
*!   Data-only constant: valid with ties, strata, late entry, pweights.
*!   lrtest / estat ic now comparable with stcox. The EM convergence
*!   criterion (ltolerance) is now evaluated on this scale.
*! - Comparison model (plain Cox, Breslow, same sample/strata/weights/X)
*!   fitted in Mata (Newton-Raphson, stcoxem_cox_ll) -> e(ll_c),
*!   e(chi2_c) = max(0, 2(ll-ll_c)), e(p_c) = chi2tail(1,chi2_c)/2, exactly
*!   as stcox_fr; "LR test of theta=0: chibar2(01)" line under the table.
*!   Option nolrtest (stcox's) skips it. Hidden e(ll_c_conv).
*! REV 49 (2026-11-26): e() and display in the layout of stcox, shared().
*! - e(b)/e(V): full factor-variable stripe with base levels (0b.drug) and
*!   omitted terms as zero columns, posted with depname(_t) buildfvinfo.
*! - Collinear or constant covariates are omitted ("note: x omitted because
*!   of collinearity", o. in e(b)) instead of error 459 (rev 47b); same
*!   centred Mata check. Base detection covers interactions (1.a#0b.b).
*! - e() as stcox_fr: N N_sub N_fail risk (st_hc) N_g g_min g_avg g_max
*!   theta se_theta ll converged chi2 df_m rank (Wald test of all coefs,
*!   _post_vce_rank), shared method ties t0 chi2type cmdline datasignature
*!   (+vars); hidden crittype fr_title marginsprop marginsfootnote.
*!   Renamed: N_orig -> N, N_subj -> N_sub, sharedvar -> shared,
*!   stratavar -> strata, n_strata -> N_strata, wexp = "= var" + wtype;
*!   e(title) dropped; new e(estimator)="em". Under vce(bootstrap):
*!   e(vcetype)="Bootstrap", e(N_reps); e(se_theta) is the bootstrap SE
*!   (user decision: plse no longer runs by default under vce(bootstrap);
*!   rev 45-48 let plse's SE replace it).
*!   e(cmd) stayed "stcoxem" until rev 51 (postestimation, Path A).
*! - Display: stcox's header (st_hcd), _coef_table (eform, base levels,
*!   cformat() etc.), theta row as stcox_fr (+ plci CI), notes; replay by
*!   typing stcoxem [, level() nohr noheader nocoef]. vce(none): plain table.
*! - Mata results go to stcoxem_b/stcoxem_V (rev <=48 overwrote matrices
*!   named b and V); r() copied to scalars right after the Mata call.
*! - The data sort order is preserved (sortpreserve).
*! - Bug fix: `stcoxem, shared(id)` (no covariates) used every variable in
*!   the data as covariates (optional varlist without default=none = _all;
*!   present since the varlist became optional).
*! REV 48 (2026-11-25): vce(cluster clustvar) added; vce(robust) ==
*! vce(cluster <shared() var>). Every shared() group must lie within one
*! cluster (else error 498). Robust/cluster beta sandwich now multiplied by
*! Gc/(Gc-1), Gc = number of vce clusters (Stata convention; rev <=47c had
*! no factor). theta's SE stays model-based (plse / Louis). e(N_g) = number
*! of frailty groups (was e(N_clust)); under robust/cluster e(N_clust),
*! e(clustvar), e(vcetype)="Robust" and e(V_modelbased) (the oim V) are
*! posted. vce() abbreviations r/rob.., cl/clu.., boot.. accepted.
*! REV 47c (2026-11-25): default vce() is oim, as in stcox (user decision;
*! rev 46-47b defaulted to robust). vce(robust) must now be requested; plse
*! still runs by default under vce(oim)/vce(robust).
*! REV 47b (2026-11-25): collinearity check moved from _rmcoll to Mata on
*! centred, unit-length columns (_rmcoll called age+1e6 collinear with the
*! constant: its own precision limit on uncentred data, test t7);
*! e(se_theta_source) set by if/else instead of ereturn local = cond().
*! REV 47a (2026-11-25): precision and reproducibility.
*! - Ties: Mata order() breaks ties randomly (a different permutation on
*!   every call, confirmed in log), so identical calls differed by ~1e-12.
*!   The risk-set sweep now sorts with the row number as tie-breaker:
*!   identical calls on the same data give identical results.
*! - Convergence (Stata convention): parameter change measured as
*!   |new-old|/(|old|+1) (a coefficient near 0 no longer stalls), AND the
*!   log likelihood's relative change must be below ltolerance() (new,
*!   default 1e-7). Default tolerance() is now 1e-8 (was 1e-6): EM converges
*!   linearly, and at 1e-6 the estimates stopped ~3e-6 (relative) and e(ll)
*!   ~1.7e-4 short of the fixed point (confirmed in log). beta's previous
*!   value is now taken at the start of each EM cycle (matters for coxiter>1).
*! - Covariates are centred internally (column means of the estimation
*!   sample): beta, theta, e(ll) and all variances are unchanged in exact
*!   arithmetic; the baseline (generate(), e(H0mean)) is rescaled back to
*!   x=0. Reduces cancellation in the running risk-set sums at scale; the
*!   sums are also reset to exactly 0 whenever the risk set empties.
*! - Collinear or constant covariates: error 459 naming them
*!   instead of invsym() silently freezing the coefficient; a coefficient
*!   that becomes unidentified during the EM also stops with error 459.
*! - plci's chi2 critical value is computed in Mata; plse's SE is passed
*!   through a tempname scalar, not a local (full double precision).
*! - log[Gamma(a+d)/Gamma(a)]: asymptotic series only where d^2 < a/1000.
*! REV 46 (superseded by 47c): default vce() changed from none to robust (cluster-robust beta
*! sandwich over shared(); theta's SE from plse, which runs by default under
*! vce(robust)). A plain -stcoxem x, shared(id)- now reports SEs, z, p and
*! CIs. Large production batch runs must now ask for vce(none) explicitly to
*! get the fastest fit (no sensitivity pass, no plse refits).
*! REV 45: plse no longer runs by default under vce(none). It previously ran
*! (5 extra profile-likelihood refits) regardless of vce(), but e(se_theta)
*! is only ever created under vce(oim)/vce(robust) -- under vce(none) the
*! result landed solely in e(se_theta_profile), unused and undisplayed. Now
*! doplse defaults to 0 whenever vce(none) (the recommended production
*! setting); noplse still forces it off under oim/robust as before, and
*! nothing changes for vce(oim)/vce(robust)/vce(bootstrap).
*! REV 44: console-output cleanup, part 2. Removed the "se_theta source: ..."
*! block (plse-vs-oim provenance note) and the "(coefficients shown as
*! hazard ratios ...)" line from the results display; removed the
*! "Profile-likelihood se_theta: ..." printf (all 4 branches: success and
*! the 3 could-not-compute variants) from the plse block. All underlying
*! values are unchanged and still available via e() (e(se_theta),
*! e(se_theta_profile), e(se_theta_oim), e(theta_profile_nonconv)) -- only
*! the automatic printing was removed.
*! REV 43: markout ... strok added. -markout- treats ANY string variable in
*! its varlist as unconditionally missing unless strok is given, so a
*! string strata() (Test D) silently zeroed out `touse' ("no observations").
*! Would have hit a string shared() the same way; strok covers that too.
*! REV 42: strata() bug fix: numeric strata codes are now always remapped to
*! 1..S (previously a stratum coded 0, negative, fractional or with gaps was
*! silently dropped from all risk sets).
*! REV 41: Louis/oim block rescaled to the same cluster-replication semantics
*! as rev 40 (a cluster with weight Wi = Wi independent copies): Mtt uses Wi
*! (was Wi^2, rev 16, which assumed a multiplicative score weight), Mbb gets
*! /Wi (infl_clust already carries Wi), Mbt loses its extra Wi factor. Before
*! this, weighted vce(oim) gave ~0 beta SEs and a missing se_theta (Test 3).
*! vce(robust) is the standard sampling-weight sandwich (wt in the score) and
*! is intentionally NOT the same as the replicated-data model-based SE.
*! REV 40: pweight point-estimate fix. Hi/Di were accumulated as wt*H and
*! wt*d (inflating the frailty posterior by the weight; confirmed by Test 2:
*! theta 6.34 vs 3.62 for replicated clusters). They are now
*! sum(wt_r/Wi * H_r) and sum(wt_r/Wi * d_r), Wi = cluster mean weight, so
*! cluster-constant weight k reproduces k replicated clusters exactly. The
*! e(ll) cluster terms are scaled by Wi. The Louis/oim block (Mbb, Mbt, Mtt,
*! Wi^2 from rev 16) is NOT yet reconciled with this replication semantics:
*! under pweight prefer vce(robust)/vce(bootstrap) until Test 3 is resolved.
*! REV 39: pweights are now read from -stset- (stset ... [pweight=var]),
*! NOT from stcoxem's own syntax. Rev 38 declared [pweight] in -syntax-,
*! which leaves `exp' as "=var" (leading "=") and errored with "=var
*! invalid name"; and any weight given to -stset- was silently ignored.
*! Now: char _dta[st_wt]/_dta[st_wv] are read; only pweight is accepted;
*! a weight -stset- has but that cannot be read is an error, never a silent
*! fallback to unit weights. Every fit prints "(weights: ...)".
*! NOTE: Hi/Di are still accumulated as wt*H and wt*d (record-frequency
*! interpretation). Whether that is right for cluster-level pweights is
*! under test (stcoxem_test_pweight.do, Test 2: weights vs replicated clusters).
*!
*! Unlike stpexpfrailty (piecewise-exponential baseline), this estimates the
*! baseline hazard nonparametrically (Breslow step function, as in -stcox-).
*! Built for speed at scale: O(n) risk-set sweeps, Breslow ties, one Newton
*! step for beta per EM cycle by default, no SE computation unless asked.
*!
*! SYNTAX
*!   stcoxem [indepvars] [if] [in] , shared(varname)
*!        [ strata(varname) vce(oim|robust|cluster clustvar|none|bootstrap) coxiter(#)
*!          generate(newvar) maxiter(#) tolerance(#) ltolerance(#) nolog dots
*!          accelerate(squarem|aitken|none) noaitken level(#) nohr noheader nocoef
*!          reps(#) bsseed(#) bswarmstart plci noplse ufrailty(newvar) ]
*!   vce() default: oim, as in stcox (rev 47c; was robust in rev 46-47b).
*!   Weights: set with -stset [pweight=var]- (pweight only).
*!
*! HISTORY (condensed; full narrative was in rev 38's header)
*!   rev 38: vce(robust) beta sandwich at-risk term now includes the frailty
*!           offset (r0_robust = wt*exp(x'b+offset)); infl_rec untouched.
*!   rev 37: noplse fixed (-syntax- populates local `plse' for a declared
*!           noPLse option, compared to the string "noplse").
*!   rev 34: ufrailty(newvar): record-level posterior mean frailty.
*!   rev 33: e(se_theta) reports plse's profile SE when it ran reliably;
*!           analytic Louis value always kept in e(se_theta_oim).
*!   rev 32: plse on by default (noplse to skip).
*!   rev 28-31: plse (numerical 2nd derivative of profile log-lik in theta,
*!           two step sizes + Richardson, tight internal tolerance, center
*!           point also refit); rev 31 fixed the e(ll) formula (removed the
*!           spurious sum(Di*ln(w1)) term); e(ll) is correct up to an
*!           additive tie-pattern constant (fine within-model, not for
*!           cross-model LR/AIC).
*!   rev 27: factor variables: fvexpand (labels) and fvrevar default mode
*!           (columns) walked together; base levels dropped from both.
*!   rev 26: plci (profile-likelihood CI for theta).
*!   rev 25: accelerate(): SQUAREM default, aitken, none.
*!   rev 24: bootstrap cold start default; bswarmstart opts in.
*!   rev 21: resync Hi/Di/H0diff/offset after an accelerated jump.
*!   rev 17: vce(bootstrap): cluster bootstrap over shared() groups.
*!   rev 16: Mtt uses Wi^2.
*!   Known, explained: vce(oim)/vce(robust) se_theta sits ~13-20% below
*!   stcox/bootstrap/frailtyEM (semiparametric-information deficit, Louis
*!   decomposition treats the Breslow baseline as known); Ctt/Mtt proven
*!   exact at EM convergence. Use plse/plci/bootstrap for theta's SE.
*!
*! DATA: must be -stset- (_t0 _t _d _st); counting-process data and left
*! truncation supported; pre-split time-varying covariates yourself.
*!
*! MODEL:  h(t | x, u_i) = u_i * h0_strat(t) * exp(x'beta)
*!         u_i ~ Gamma(mean 1, variance theta), shared within shared()-group
*!
*! ALGORITHM (Klein 1992): E-step closed-form gamma-Poisson posterior
*! moments from Hi=sum[H0(t)-H0(t0)]exp(x'b); M-step: offset log E[u_i] in
*! an offset Cox partial likelihood (beta by Newton, Breslow baseline in one
*! O(n) sweep per stratum), theta by a stable Newton step.
*!
*! LIMITATIONS: Breslow ties only; reported SEs are approximations except
*! bootstrap; e(mem_est) is a lower bound.

program define stcoxem, eclass sortpreserve
    version 19.0
    * rev 51: hidden subcommand for stcoxem_p -- ln u at e(b), e(theta)
    if `"`1'"' == "_lnu" {
        gettoken sub 0 : 0
        LnU `0'
        exit
    }
    if replay() {
        * rev 49: replay; rev 51 will post e(cmd)=stcox_fr (Path A), so a
        * stcoxem fit is also recognised by e(estimator)=="em". A model
        * without covariates (stcoxem, shared(id) ...) is estimated, as is
        * any call with -estimate- (stcox convention).
        local zero : copy local 0
        syntax [, ESTImate SHared(passthru) *]
        if "`estimate'" == "" & `"`shared'"' == "" {
            if "`e(cmd)'" != "stcoxem" & !("`e(cmd)'" == "stcox_fr" & "`e(estimator)'" == "em") {
                error 301
            }
            Display `zero'
            exit
        }
        local 0 : copy local zero
    }
    Estimate `0'
end


program define Estimate, eclass
    version 19.0
    local cmdline : copy local 0

    * rev 49: default=none -- an optional varlist left empty otherwise
    * means _all (stcoxem, shared(id) took every variable; test log 003758)
    syntax [varlist(default=none numeric fv ts)] [if] [in] , ///
        SHared(varname) ///
        [ STRata(varname) VCE(string) COXiter(integer 1) GENerate(name) ///
          MAXiter(integer -1) TOLerance(real 1e-8) LTOLerance(real 1e-7) noLOg DOTS ///
          ACCELerate(string) noAItken LEVEL(cilevel) noHR noHEADer noCOEF ///
          REPS(integer 50) BSSEED(integer -1) BSWarmstart PLci noPLse ///
          UFrailty(name) noLRtest ESTImate FORCEShared * ]

    * rev 49: display options (cformat() etc.) as in stcox
    _get_diopts diopts options, `options'
    if `"`options'"' != "" {
        local w : word 1 of `options'
        di as err `"option {bf:`w'} not allowed"'
        exit 198
    }

    * rev 54: default maxiter = c(maxiter), as stcox
    if `maxiter' < 0 local maxiter = c(maxiter)
    if `maxiter' < 1 {
        di as err "maxiter() must be a positive integer"
        exit 198
    }

    if "`log'" != "" & "`dots'" != "" {
        di as err "may not specify both nolog and dots"
        exit 198
    }

    * ---- must already be stset -----------------------------------------
    foreach v in _t0 _t _d _st {
        capture confirm variable `v'
        if _rc {
            di as err "data not stset; run -stset- (with -id()- for" ///
                " multiple-record data) before -stcoxem-"
            exit 119
        }
    }

    local xvars `varlist'

    * rev 47c: default is oim, as in stcox (rev 46-47b: robust; through rev 45: none)
    * rev 48: vce(cluster clustvar); abbreviations r/rob.., cl/clu.., boot..
    local clustvar
    gettoken vcetype vcearg : vce
    local vcearg = strtrim(`"`vcearg'"')
    local lt = length("`vcetype'")
    if `lt' == 0 local vcetype "oim"
    else if "`vcetype'" == substr("robust", 1, `lt') local vcetype "robust"
    else if `lt' >= 2 & "`vcetype'" == substr("cluster", 1, `lt') local vcetype "cluster"
    else if `lt' >= 4 & "`vcetype'" == substr("bootstrap", 1, `lt') local vcetype "bootstrap"
    if !inlist("`vcetype'", "none", "oim", "robust", "cluster", "bootstrap") {
        di as err "vce() must be oim, robust, cluster clustvar, none, or bootstrap"
        exit 198
    }
    if "`vcetype'" == "cluster" {
        if `"`vcearg'"' == "" {
            di as err "vce(cluster) requires a cluster variable"
            exit 198
        }
        unab clustvar : `vcearg', max(1)
    }
    else if `"`vcearg'"' != "" {
        di as err "vce(`vce') invalid"
        exit 198
    }
    local vce "`vcetype'"
    local vcemode = cond("`vce'"=="none", 0, cond("`vce'"=="oim", 1, cond(inlist("`vce'", "robust", "cluster"), 2, 3)))
    if "`vce'"=="bootstrap" & `reps' < 2 {
        di as err "reps() must be at least 2 for vce(bootstrap)"
        exit 198
    }

    * ---- weights: taken from -stset- ---------------------------------------
    local wtype : char _dta[st_wt]
    local wvar  : char _dta[st_wv]
    local wfull : char _dta[st_w]
    if "`wfull'" != "" & ("`wtype'" == "" | "`wvar'" == "") {
        di as err "could not read the weight from -stset- (st_w=`wfull')"
        exit 198
    }
    if "`wtype'" != "" & "`wtype'" != "pweight" {
        di as err "only pweights are allowed (data were -stset- with `wtype')"
        exit 101
    }
    tempvar wt
    if "`wtype'" == "pweight" {
        quietly gen double `wt' = `wvar'
        di as txt "(weights: pweight=`wvar' from -stset-)"
    }
    else {
        quietly gen byte `wt' = 1
    }

    * ---- sample ----------------------------------------------------------
    marksample touse
    quietly replace `touse' = 0 if _st != 1
    * rev 43: strok -- markout otherwise treats ANY string variable in its
    * varlist as missing, which zeroed out `touse' for a string strata()/shared().
    markout `touse' `shared' `strata' `clustvar' _t0 _t _d `wt', strok
    quietly count if `touse'
    if r(N) == 0 {
        di as err "no observations"
        exit 2000
    }

    * rev 53: delayed entry / gaps with shared(): stcox's own check (498 unless
    * forceshared; with forceshared a warning)
    _st_err_sharedgaps "`shared'" "`forceshared'" "`touse'"

    * ---- factor-variable / time-series expansion ------------------------
    * xvars_full (fvexpand) and xvars_full_revar (fvrevar, DEFAULT mode) are
    * built from the same expression, one token per term in the same order.
    * rev 49: the full list (base levels included) is kept as the e(b)
    * stripe, as in stcox; `fitpos' holds the positions of the estimated
    * terms. A term is a base level if ANY of its #-components is (e.g.
    * 1.drug#0b.DM), and fvexpand's own omitted terms (o.) are not fitted.
    local K 0
    local fullnames
    local fitpos
    local xvars_expand
    local xvars_revar
    if "`xvars'" != "" {
        quietly fvexpand `xvars' if `touse'
        local xvars_full `r(varlist)'
        quietly fvrevar `xvars' if `touse'
        local xvars_full_revar `r(varlist)'
        local nfull       : word count `xvars_full'
        local nfull_revar : word count `xvars_full_revar'
        if `nfull' != `nfull_revar' {
            di as err "could not reliably match factor-variable names to" ///
                " estimated columns (this can happen with empty cells or" ///
                " unusual interactions). Please pre-generate dummy/interaction" ///
                " variables yourself (e.g. via -tabulate var, gen()- or" ///
                " -fvrevar-) and pass them as a plain numeric varlist instead."
            exit 198
        }
        local K `nfull'
        local fullnames `xvars_full'
        forvalues i = 1/`nfull' {
            local v  : word `i' of `xvars_full'
            local rv : word `i' of `xvars_full_revar'
            local isbase = ustrregexm("`v'", "(^|#)[0-9]+b\.")
            local isomit = ustrregexm("`v'", "(^|#)([0-9]+b?|c)?o[A-Za-z0-9]*\.")
            if !`isbase' & !`isomit' {
                local xvars_expand `xvars_expand' `v'
                local xvars_revar  `xvars_revar'  `rv'
                local fitpos       `fitpos' `i'
            }
        }
    }
    local p : word count `xvars_revar'

    * ---- rev 49: collinear or constant covariates are omitted (o. in
    *      e(b)/e(V), "note: ... omitted because of collinearity", as stcox)
    *      instead of error 459 (rev 47b). Checked in Mata on centred,
    *      unit-length columns (_rmcoll works on uncentred data and flagged a
    *      covariate shifted by 1e6 as collinear with the constant). Mata
    *      returns the positions (within the fitted terms) in local collinpos.
    if `p' > 0 {
        local collinpos
        mata: stcoxem_collin("`xvars_revar'", "`touse'")
        if "`collinpos'" != "" {
            local keep_e
            local keep_r
            local keep_p
            forvalues j = 1/`p' {
                local nm : word `j' of `xvars_expand'
                local rv : word `j' of `xvars_revar'
                local fp : word `j' of `fitpos'
                local isin : list j in collinpos
                if `isin' {
                    di as txt "note: `nm' omitted because of collinearity."
                    local oname
                    OmitName `nm'
                    local newfull
                    forvalues k = 1/`K' {
                        local w : word `k' of `fullnames'
                        if `k' == `fp' local w `oname'
                        local newfull `newfull' `w'
                    }
                    local fullnames `newfull'
                }
                else {
                    local keep_e `keep_e' `nm'
                    local keep_r `keep_r' `rv'
                    local keep_p `keep_p' `fp'
                }
            }
            local xvars_expand `keep_e'
            local xvars_revar  `keep_r'
            local fitpos       `keep_p'
            local p : word count `xvars_revar'
        }
    }

    * ---- numeric group ids for shared() and strata() -------------------
    capture confirm numeric variable `shared'
    if _rc {
        tempvar sharedid
        quietly egen long `sharedid' = group(`shared') if `touse'
    }
    else local sharedid `shared'

    * ---- rev 48: vce(cluster): numeric cluster ids; every shared() group
    *      must lie within one cluster (as panels within clusters in xt) ----
    local clid
    if "`vce'" == "cluster" {
        tempvar clid
        quietly egen long `clid' = group(`clustvar') if `touse'
        local nestok 0
        local nclus 0
        mata: stcoxem_nested("`sharedid'", "`clid'", "`touse'")
        if !`nestok' {
            di as err "shared() groups are not nested within clusters"
            di as err "every group of `shared' must lie within one cluster of `clustvar'"
            exit 498
        }
        if `nclus' < 2 {
            di as err "vce(cluster `clustvar') requires at least 2 clusters"
            exit 459
        }
    }

    if "`strata'" != "" {
        * rev 42: ALWAYS remap to consecutive codes 1..S (the Mata code loops
        * s=1..S and selects strat==s).
        tempvar stratid
        quietly egen long `stratid' = group(`strata') if `touse'
        quietly summarize `stratid' if `touse'
        local S = r(max)
    }
    else {
        tempvar stratid
        quietly gen byte `stratid' = 1
        local S = 1
    }

    if "`log'" == "" & "`dots'" == "" local dolog 1
    else                              local dolog 0
    local dodots     = ("`dots'" != "")

    * ---- accelerate(): squarem (default), aitken, none; noaitken = none ---
    if "`accelerate'" != "" & "`aitken'" == "noaitken" & "`accelerate'" != "none" {
        di as err "may not combine noaitken with accelerate(`accelerate'')"
        exit 198
    }
    if "`accelerate'" == "" {
        local accelerate = cond("`aitken'"=="noaitken", "none", "squarem")
    }
    if !inlist("`accelerate'", "squarem", "aitken", "none") {
        di as err "accelerate() must be squarem, aitken, or none"
        exit 198
    }
    local accelmode = cond("`accelerate'"=="none", 0, cond("`accelerate'"=="aitken", 1, 2))
    local usebswarm  = ("`bswarmstart'" != "")

    * ---- plci(): profile-likelihood CI for theta ----
    local doplci = ("`plci'" != "")

    * ---- plse: default ON under vce(oim)/robust/cluster; noplse (declared
    *      noPLse fills local `plse' with "noplse", rev 37) switches it off.
    *      rev 49 (user decision 2026-11-26): OFF under vce(none) and
    *      vce(bootstrap) -- under bootstrap e(se_theta) is the bootstrap SE
    *      (rev 45-48: plse ran and replaced it) ----
    local doplse = cond(`vcemode' == 0 | `vcemode' == 3, 0, cond("`plse'" == "noplse", 0, 1))

    * ---- rev 50: LR test of theta=0 against plain Cox (as stcox_fr; off
    *      with nolrtest, stcox_fr's option) ----
    local dolrtest = ("`lrtest'" != "nolrtest")

    * ---- subject identifier for H0mean, from -stset-'s id() ----
    local stid : char _dta[st_id]
    if "`stid'" == "" {
        tempvar subjid
        quietly gen long `subjid' = _n if `touse'
        local nosubjid 1
    }
    else {
        capture confirm numeric variable `stid'
        if _rc {
            tempvar subjid
            quietly egen long `subjid' = group(`stid') if `touse'
        }
        else local subjid `stid'
        local nosubjid 0
    }

    * ---- optional per-record baseline-hazard output variable ----------
    local dogen = ("`generate'" != "")
    if `dogen' {
        confirm new variable `generate'
        quietly gen double `generate' = .
    }
    else {
        tempvar generate
        quietly gen double `generate' = .
    }

    * ---- optional per-record posterior-mean-frailty output (rev 34) ----
    local doufrailty = ("`ufrailty'" != "")
    if `doufrailty' {
        confirm new variable `ufrailty'
        quietly gen double `ufrailty' = .
    }
    else {
        tempvar ufrailty
        quietly gen double `ufrailty' = .
    }

    * ---- rev 49: Mata writes its results to stcoxem_b/stcoxem_V/stcoxem_Vmb
    *      (rev <=48: b and V, which overwrote any user matrices so named);
    *      drop leftovers of an interrupted run first ----
    foreach m in stcoxem_b stcoxem_V stcoxem_Vmb {
        capture matrix drop `m'
    }

    * ================= hand off to Mata =================================
    mata: stcoxem_driver( ///
        "`xvars_revar'", "`sharedid'", "`stratid'", "`wt'", "`touse'", "`subjid'", ///
        `S', `maxiter', `coxiter', `tolerance', `ltolerance', `dolog', `dodots', `accelmode', ///
        `vcemode', "`generate'", `reps', `bsseed', `usebswarm', `doplci', `level', ///
        `doplse', "`ufrailty'", "`clid'", `dolrtest')

    * ---- copy r() at once: later commands overwrite it ----
    local rlist theta ll N_orig N_subj N_clust N_clust_vce H0mean iter converged ///
        mem_est_mb se_theta se_theta_profile theta_profile_nonconv se_theta_naive ///
        schur_term theta_pl_lower theta_pl_upper theta_pl_nonconv aitken_theta_att ///
        aitken_theta_acc aitken_beta_att aitken_beta_acc squarem_att squarem_acc ///
        boot_reps boot_reps_used boot_reps_maxiter boot_reps_badtheta ///
        boot_reps_badbeta boot_avg_iter
    if `dolrtest' local rlist `rlist' ll_c ll_c_conv
    foreach s of local rlist {
        tempname R_`s'
        scalar `R_`s'' = r(`s')
    }
    tempname bfit Vfit Vmbfit
    local haveb 0
    local haveV 0
    local haveVmb 0
    if `p' > 0 {
        matrix `bfit' = stcoxem_b
        local haveb 1
    }
    capture confirm matrix stcoxem_V
    if !_rc {
        matrix `Vfit' = stcoxem_V
        local haveV 1
    }
    capture confirm matrix stcoxem_Vmb
    if !_rc {
        matrix `Vmbfit' = stcoxem_Vmb
        local haveVmb 1
    }
    foreach m in stcoxem_b stcoxem_V stcoxem_Vmb {
        capture matrix drop `m'
    }

    if `dogen' {
        label variable `generate' "record-level Breslow baseline cumulative hazard, stcoxem"
    }
    if `doufrailty' {
        label variable `ufrailty' "record-level posterior mean frailty E[u|data] (shared() cluster), stcoxem"
    }

    * ---- group information (records per shared() group, as stcox_fr) ----
    tempname g_N g_min g_max g_avg
    tempvar gT
    quietly bysort `touse' `sharedid': gen long `gT' = _N if _n == 1 & `touse'
    quietly summarize `gT' if `touse', meanonly
    scalar `g_N'   = r(N)
    scalar `g_min' = r(min)
    scalar `g_max' = r(max)
    scalar `g_avg' = r(mean)
    * subjects (stset id()) spread over several frailty groups -> warning, as stcox_fr
    local shwarn 0
    if "`stid'" != "" {
        tempvar dif
        quietly bysort `touse' `stid' (`sharedid'): gen byte `dif' = ///
            (`sharedid'[1] != `sharedid'[_N]) if `touse'
        quietly count if `dif' == 1
        if r(N) > 0 local shwarn 1
    }

    * ---- rev 49: e(b)/e(V) on the full stripe (base and omitted terms as
    *      0 columns, as stcox) ----
    tempname bF VF VmbF
    if `K' > 0 {
        local bnarg = cond(`haveb', "`bfit'", "")
        local Vnarg = cond(`haveV', "`Vfit'", "")
        local Voarg = cond(`vcemode' >= 1, "`VF'", "")
        mata: stcoxem_expandbv("`bnarg'", "`Vnarg'", "`fitpos'", `K', "`bF'", "`Voarg'")
        matrix colnames `bF' = `fullnames'
        if `vcemode' >= 1 {
            matrix colnames `VF' = `fullnames'
            matrix rownames `VF' = `fullnames'
        }
        if `vcemode' == 2 {
            local Vnarg = cond(`haveVmb', "`Vmbfit'", "")
            mata: stcoxem_expandbv("", "`Vnarg'", "`fitpos'", `K', "`VmbF'", "`VmbF'")
            matrix colnames `VmbF' = `fullnames'
            matrix rownames `VmbF' = `fullnames'
        }
    }

    * ================= ereturn (layout of stcox, shared()) ================
    local Nobs = `R_N_orig'
    if `K' > 0 & `vcemode' >= 1 {
        ereturn post `bF' `VF', depname(_t) obs(`Nobs') esample(`touse') buildfvinfo ADDCONS
    }
    else if `K' > 0 {
        ereturn post `bF', depname(_t) obs(`Nobs') esample(`touse') buildfvinfo ADDCONS
    }
    else {
        ereturn post, obs(`Nobs') esample(`touse')
        ereturn local depvar "_t"
    }
    tempvar esamp
    quietly gen byte `esamp' = e(sample)

    * weights before signestimationsample (it signs the weight variable too)
    if "`wtype'" == "pweight" {
        ereturn local wtype "pweight"
        ereturn local wexp  "= `wvar'"
    }

    * -- macros as stcox_fr; e(cmd)="stcox_fr" (posted last, below) and
    *    e(cmd2)/e(predict)/e(estat_cmd) as stcox_fr (rev 51)
    ereturn local estimator "em"
    ereturn local shared    "`shared'"
    ereturn local method    "breslow"
    ereturn local ties      "breslow"
    ereturn local t0        "_t0"
    ereturn local vce       "`vce'"
    if `vcemode' == 2 ereturn local vcetype "Robust"
    if `vcemode' == 3 ereturn local vcetype "Bootstrap"
    ereturn hidden local crittype        "log likelihood"
    ereturn hidden local fr_title        "Gamma shared frailty"
    ereturn hidden local marginsprop     addcons
    ereturn hidden local marginsfootnote _multirecordcheck
    if "`strata'" != "" {
        ereturn local  strata   "`strata'"
        ereturn scalar N_strata = `S'
    }
    if `shwarn' ereturn local sh_warn sh_warn

    * -- vce clusters (robust = clusters are the shared() groups)
    if `vcemode' == 2 {
        ereturn scalar N_clust = `R_N_clust_vce'
        if "`vce'" == "cluster" ereturn local clustvar "`clustvar'"
        else                    ereturn local clustvar "`shared'"
        if `K' > 0 ereturn matrix V_modelbased = `VmbF'
    }

    * -- Wald test of all coefficients and rank, as stcox_fr
    if `K' > 0 & `vcemode' >= 1 {
        ereturn local chi2type "Wald"
        tempname bW
        matrix `bW' = e(b)
        local wnames : colfullnames `bW'
        capture test `wnames'
        if !_rc {
            ereturn scalar chi2 = r(chi2)
            ereturn scalar df_m = r(df)
        }
        else {
            ereturn scalar chi2 = 0
            ereturn scalar df_m = 0
        }
        _post_vce_rank
    }
    else {
        ereturn scalar df_m = `p'
        ereturn scalar rank = `p'
    }

    * -- subjects, failures, time at risk (stcox's own routine; weighted
    *    under pweights, as stcox)
    st_hc `esamp'

    * -- groups
    ereturn scalar N_g   = `g_N'
    ereturn scalar g_min = `g_min'
    ereturn scalar g_avg = `g_avg'
    ereturn scalar g_max = `g_max'

    * -- theta and its SE: under vce(bootstrap) the bootstrap SE (rev 49);
    *    otherwise plse's profile SE when it ran and converged, else the
    *    analytic Louis value (oim/robust/cluster)
    ereturn scalar theta = `R_theta'
    if `vcemode' == 3 {
        ereturn scalar se_theta = `R_se_theta'
        ereturn local  se_theta_source "bootstrap"
    }
    else if `vcemode' >= 1 {
        ereturn scalar se_theta_oim = `R_se_theta'
        local altsrc "oim"
        if `doplse' & !missing(`R_se_theta_profile') {
            ereturn scalar se_theta = `R_se_theta_profile'
            ereturn local  se_theta_source "plse"
        }
        else {
            ereturn scalar se_theta = `R_se_theta'
            if `doplse' ereturn local se_theta_source "`altsrc' (plse did not converge reliably)"
            else        ereturn local se_theta_source "`altsrc' (noplse)"
        }
    }
    ereturn scalar ll        = `R_ll'
    * rev 50: LR test of theta=0 (chibar2(01)), exactly as stcox_fr
    if `dolrtest' & !missing(`R_ll_c') {
        ereturn scalar ll_c   = `R_ll_c'
        ereturn scalar chi2_c = max(0, 2*(e(ll)-e(ll_c)))
        ereturn scalar p_c    = chi2tail(1, e(chi2_c))*0.5
        ereturn hidden scalar ll_c_conv = `R_ll_c_conv'
    }
    ereturn scalar converged = `R_converged'

    * -- stcoxem's own results
    ereturn scalar iter       = `R_iter'
    ereturn scalar H0mean     = `R_H0mean'
    ereturn scalar mem_est    = `R_mem_est_mb'
    ereturn scalar tolerance  = `tolerance'
    ereturn scalar ltolerance = `ltolerance'
    ereturn local  accelerate "`accelerate'"
    if `accelmode' == 1 {
        ereturn scalar aitken_theta_att = `R_aitken_theta_att'
        ereturn scalar aitken_theta_acc = `R_aitken_theta_acc'
        if `p' > 0 {
            ereturn scalar aitken_beta_att = `R_aitken_beta_att'
            ereturn scalar aitken_beta_acc = `R_aitken_beta_acc'
        }
    }
    else if `accelmode' == 2 {
        ereturn scalar squarem_att = `R_squarem_att'
        ereturn scalar squarem_acc = `R_squarem_acc'
    }
    if (`vcemode'==1 | `vcemode'==2) & `p' > 0 {
        ereturn scalar se_theta_naive = `R_se_theta_naive'
        ereturn scalar schur_term     = `R_schur_term'
    }
    if `doplci' {
        ereturn scalar theta_pl_lower   = `R_theta_pl_lower'
        ereturn scalar theta_pl_upper   = `R_theta_pl_upper'
        ereturn scalar theta_pl_nonconv = `R_theta_pl_nonconv'
        ereturn hidden scalar plci_level = `level'
    }
    if `doplse' {
        ereturn scalar se_theta_profile      = `R_se_theta_profile'
        ereturn scalar theta_profile_nonconv = `R_theta_profile_nonconv'
    }
    if `vcemode' == 3 {
        ereturn scalar N_reps             = `R_boot_reps_used'
        ereturn scalar boot_reps          = `R_boot_reps'
        ereturn scalar boot_reps_used     = `R_boot_reps_used'
        ereturn scalar boot_reps_maxiter  = `R_boot_reps_maxiter'
        ereturn scalar boot_reps_badtheta = `R_boot_reps_badtheta'
        ereturn scalar boot_reps_badbeta  = `R_boot_reps_badbeta'
        ereturn scalar boot_avg_iter      = `R_boot_avg_iter'
        if !`usebswarm' ereturn local boot_start "cold"
        else            ereturn local boot_start "warm"
        tempname boot_theta_all boot_iter_all boot_outcome_all
        matrix `boot_theta_all'   = boot_theta_all
        matrix `boot_iter_all'    = boot_iter_all
        matrix `boot_outcome_all' = boot_outcome_all
        ereturn matrix boot_theta_all   = `boot_theta_all'
        ereturn matrix boot_iter_all    = `boot_iter_all'
        ereturn matrix boot_outcome_all = `boot_outcome_all'
        if `p' > 0 {
            tempname boot_beta_all
            matrix `boot_beta_all' = boot_beta_all
            matrix colnames `boot_beta_all' = `xvars_expand'
            ereturn matrix boot_beta_all = `boot_beta_all'
        }
    }

    * -- data signature (needed by stcox's predict, rev 51); variables as
    *    stcox_fr: _t _t0 _d, the covariates without operators, shared()
    local sigvars
    if "`xvars_expand'" != "" {
        unopvarlist `xvars_expand'
        local sigvars `r(varlist)'
    }
    local sigall _t _t0 _d `sigvars' `shared' `strata' `clustvar'
    local sigall : list uniq sigall
    signestimationsample `sigall'

    ereturn local cmdline `"stcoxem `cmdline'"'

    * -- rev 51: postestimation. Unweighted: stcox's own predict (Path A).
    *    Weighted: stcoxem_p (frailties from stcoxem; refit statistics on a
    *    plain-Cox view with offset ln u -- Path B)
    ereturn local estat_cmd "stcox_estat"
    if "`wtype'" == "pweight" ereturn local predict "stcoxem_p"
    else                      ereturn local predict "stcox_p"
    ereturn local cmd2      "stcox"
    ereturn local cmd       "stcox_fr"

    if `nosubjid' {
        di as txt "(no id() in -stset-: each record is treated as its own" ///
            " subject for e(H0mean) and e(N_sub))"
    }
    Display, level(`level') `hr' `header' `coef' `diopts'
end


* ---------------------------------------------------------------------
* rev 51: ln u (posterior-mean log frailty) at e(b) and e(theta), the EM
* fixed point at the fitted values, stored in a new double variable for
* e(sample). Used by stcoxem_p (predict after weighted fits): effects,
* atfrailty, and the offset of the plain-Cox view for refit statistics.
* Called as -stcoxem _lnu newvar-; leaves e() unchanged.
* ---------------------------------------------------------------------
program define LnU
    syntax newvarname
    if "`e(estimator)'" != "em" | "`e(cmd2)'" != "stcox" error 301
    if "`e(datasignature)'" != "" checkestimationsample
    tempvar touse xb gid sid wt
    quietly gen byte `touse' = e(sample)
    tempname bb
    local K 0
    capture matrix `bb' = e(b)
    if !_rc local K = colsof(`bb')
    if `K' > 0 quietly _predict double `xb' if `touse', xb
    else       quietly gen double `xb' = 0 if `touse'
    quietly egen long `gid' = group(`e(shared)') if `touse'
    if "`e(strata)'" != "" quietly egen long `sid' = group(`e(strata)') if `touse'
    else                   quietly gen byte `sid' = 1 if `touse'
    if "`e(wtype)'" != "" quietly gen double `wt' `e(wexp)' if `touse'
    else                  quietly gen byte `wt' = 1 if `touse'
    quietly gen double `varlist' = .
    local lnu_conv 0
    local lnu_iter 0
    mata: stcoxem_lnu("`xb'", "`gid'", "`sid'", "`wt'", "`touse'", "`varlist'")
    if !`lnu_conv' {
        di as txt "(note: frailties at the estimates did not converge in" ///
            " `lnu_iter' iterations; predictions approximate)"
    }
end


* ---------------------------------------------------------------------
* rev 49: name of an omitted term: 1.drug -> 1o.drug, c.age -> co.age,
* x -> o.x, L.x -> oL.x, each #-component of an interaction alike.
* Returned in the caller's local oname.
* ---------------------------------------------------------------------
program define OmitName
    args term
    local out
    local rest `term'
    while "`rest'" != "" {
        gettoken c rest : rest, parse("#")
        if "`c'" == "#" {
            local out `out'#
            continue
        }
        if ustrregexm("`c'", "^([0-9]+)n?\.(.+)") {
            local c = ustrregexs(1) + "o." + ustrregexs(2)
        }
        else if ustrregexm("`c'", "^c\.(.+)") {
            local c = "co." + ustrregexs(1)
        }
        else if ustrregexm("`c'", "^([A-Za-z][A-Za-z0-9]*)\.(.+)") {
            local c = "o" + ustrregexs(1) + "." + ustrregexs(2)
        }
        else local c = "o.`c'"
        local out `out'`c'
    }
    c_local oname `out'
end


* ---------------------------------------------------------------------
* rev 49: results display (also replay): stcox's header (st_hcd),
* _coef_table, theta row, LR test of theta=0 (once e(ll_c) exists, rev 50)
* and notes -- layout of stcox_fr's Display.
* ---------------------------------------------------------------------
program define Display
    syntax [, Level(cilevel) noHR noHEADer noCOEF *]
    _get_diopts diopts options, `options'
    if `"`options'"' != "" {
        local w : word 1 of `options'
        di as err `"option {bf:`w'} not allowed"'
        exit 198
    }

    if "`header'" == "" {
        di
        di as txt "Cox regression with Breslow method for ties, EM estimation"
        st_hcd
    }
    di
    if "`coef'" != "" exit

    * is there a coefficient vector with a variance matrix?
    tempname b V
    local K 0
    capture matrix `b' = e(b)
    if !_rc local K = colsof(`b')
    local hasV 0
    if `K' > 0 {
        capture matrix `V' = e(V)
        if !_rc {
            if rowsof(`V') > 0 local hasV 1
        }
    }

    local showpl = !missing(e(theta_pl_lower)) & !missing(e(theta_pl_upper)) & ///
        (e(plci_level) == `level')

    if `hasV' {
        * (_coef_table prints "(Std. err. adjusted for # clusters in x)"
        * itself from e(N_clust)/e(clustvar) -- test log 002322)
        if "`hr'" != "" {
            _coef_table, level(`level') plus `diopts'
        }
        else {
            _coef_table, level(`level') eform("Haz. ratio") plus `diopts'
        }
        if `"`s(width_col1)'"' != "" {
            local w1 = s(width_col1)
            local wrest = `s(width)' - `w1' - 1
        }
        else {
            local w1 13
            local wrest 64
        }
    }
    else if `K' == 0 {
        * rev 52: no covariates -- header row and the theta row only
        local w1 13
        local wrest 64
        di as txt "{hline `w1'}{c TT}{hline `wrest'}"
        di as txt %12s "_t" " {c |}"
        di as txt "{hline `w1'}{c +}{hline `wrest'}"
    }
    else {
        local w1 13
        local wrest 64
        local lbl = cond("`hr'" == "", "Haz. ratio", "Coefficient")
        di as txt "{hline `w1'}{c TT}{hline `wrest'}"
        di as txt %12s "_t" " {c |} " %10s "`lbl'" "   (no standard errors: vce(none))"
        di as txt "{hline `w1'}{c +}{hline `wrest'}"
        local bnames : colfullnames `b'
        forvalues i = 1/`K' {
            local nm : word `i' of `bnames'
            local isbase = ustrregexm("`nm'", "(^|#)[0-9]+b\.")
            local isomit = ustrregexm("`nm'", "(^|#)([0-9]+b?|c)?o[A-Za-z0-9]*\.")
            local show = cond("`hr'" == "", exp(`b'[1,`i']), `b'[1,`i'])
            if `isbase' {
                di as txt %12s abbrev("`nm'", 12) " {c |} " as res %10.0g `show' as txt "  (base)"
            }
            else if `isomit' {
                di as txt %12s abbrev("`nm'", 12) " {c |} " as res %10.0g `show' as txt "  (omitted)"
            }
            else {
                di as txt %12s abbrev("`nm'", 12) " {c |} " as res %10.0g `show'
            }
        }
        di as txt "{hline `w1'}{c +}{hline `wrest'}"
    }

    * theta row: estimate and SE as stcox_fr; profile-likelihood CI (plci)
    * in the CI columns
    local w1fmt = `w1' - 1
    di as txt %`w1fmt's "theta" " {c |}  " as res %9.0g e(theta) _c
    if !missing(e(se_theta)) di as res _s(2) %9.0g e(se_theta) _c
    else                     di _s(11) _c
    if `showpl' di as res _s(21) %9.0g e(theta_pl_lower) _s(3) %9.0g e(theta_pl_upper)
    else        di
    di as txt "{hline `w1'}{c BT}{hline `wrest'}"

    Footer, `hr' showpl(`showpl') level(`level')
end


program define Footer
    syntax [, noHR showpl(integer 0) level(cilevel)]
    * LR test of theta=0 (rev 50 posts e(ll_c), e(chi2_c), e(p_c)); as stcox_fr
    if "`e(ll_c)'" != "" {
        if ((e(chi2_c) > 0.005) & (e(chi2_c) < 1e4)) | (e(chi2_c) == 0) {
            local fmt "%8.2f"
        }
        else local fmt "%8.2e"
        local chi : di `fmt' e(chi2_c)
        local chi = trim("`chi'")
        di in smcl as txt "LR test of theta=0: " as txt ///
            "{help j_chibar##|_new:chibar2(01)} = " ///
            as res "`chi'" _col(56) as txt "Prob >= chibar2 = " ///
            as res %5.3f e(p_c)
    }
    di
    local type = cond("`hr'" == "", "hazard ratios", "coefficients")
    * rev 52: no coefficient notes for a model without covariates
    tempname bchk
    local Kf 0
    capture matrix `bchk' = e(b)
    if !_rc local Kf = colsof(`bchk')
    if `Kf' == 0 {
        if "`e(vce)'" == "bootstrap" {
            di as txt "Note: Bootstrap over the shared() groups, " ///
                as res e(boot_reps_used) as txt " of " as res e(boot_reps) ///
                as txt " replications used (`e(boot_start)' start)."
        }
    }
    else if "`e(vce)'" == "oim" {
        di as txt "Note: Standard errors of `type' allow for the estimation of theta" ///
            " (Louis)."
    }
    else if inlist("`e(vce)'", "robust", "cluster") {
        di as txt "Note: Robust standard errors of `type' are conditional on the" ///
            " frailties."
    }
    else if "`e(vce)'" == "bootstrap" {
        di as txt "Note: Bootstrap over the shared() groups, " ///
            as res e(boot_reps_used) as txt " of " as res e(boot_reps) ///
            as txt " replications used (`e(boot_start)' start)."
    }
    else di as txt "Note: No standard errors computed (vce(none))."
    local src "`e(se_theta_source)'"
    if "`src'" == "plse" {
        di as txt "Note: Standard error of theta from the profile likelihood."
    }
    else if substr("`src'", 1, 3) == "oim" {
        di as txt "Note: Standard error of theta from the observed information" ///
            " (Louis);"
        di as txt "      it is known to be too small, see help stcoxem."
    }
    else if substr("`src'", 1, 9) == "bootstrap" {
        di as txt "Note: Standard error of theta from the bootstrap."
    }
    if `showpl' {
        di as txt "Note: `level'% CI of theta from the profile likelihood."
    }
    if "`e(strata)'" != "" {
        di as txt "Note: Stratified by " as res "`e(strata)'" as txt " (" ///
            as res e(N_strata) as txt " strata)."
    }
    if "`e(wtype)'" != "" {
        local wv = strtrim(substr("`e(wexp)'", 2, .))
        di as txt "Note: `e(wtype)' " as res "`wv'" as txt ///
            ": weight k counts a group as k copies."
    }
    if "`e(sh_warn)'" == "sh_warn" {
        di as txt "Warning: observations within subject belong to different" ///
            " frailty groups."
    }
    if "`e(ll_c_conv)'" == "0" {
        di as txt "Note: comparison Cox model (theta=0) did not converge;" ///
            " LR test approximate."
    }
    if missing(e(H0mean)) {
        di as txt "Note: e(H0mean) is missing: the baseline hazard at x=0 is" ///
            " out of range;"
        di as txt "      centre covariates near 0 (see help stcoxem)."
    }
    if e(converged) != 1 {
        di as err "Warning: convergence not achieved (" e(iter) " EM iterations)."
    }
end


* =====================================================================
*                              MATA CODE
* =====================================================================
mata:
mata set matastrict off

// ---------------------------------------------------------------------
// Map an arbitrary numeric id vector to consecutive group indices 1..G.
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// Collinearity check (rev 47b): centre each column, scale to unit length,
// invsym() the cross-product (a correlation-type matrix); a zero diagonal
// marks a column that is constant or a linear combination of earlier
// columns (same keep-the-first rule as _rmcoll). Positions are returned in
// the caller's local collinpos.
// ---------------------------------------------------------------------
void stcoxem_collin(string scalar xvars, string scalar touse)
{
    real matrix X, XX
    real rowvector nrm, nrm0, bad
    real scalar j
    string scalar out

    X    = st_data(., tokens(xvars), touse)
    nrm0 = sqrt(colsum(X:^2))
    X    = X :- mean(X)
    nrm  = sqrt(colsum(X:^2))
    // constant column: centred length is rounding noise of the raw length
    bad  = (nrm :<= 1e-12 :* (nrm0 :+ 1))
    for (j=1; j<=cols(X); j++) {
        if (!bad[j]) X[.,j] = X[.,j] / nrm[j]
    }
    XX = cross(X, X)
    bad = bad :| (diagonal(invsym(XX))' :== 0)
    out = ""
    for (j=1; j<=cols(X); j++) {
        if (bad[j]) out = out + " " + strofreal(j)
    }
    st_local("collinpos", out)
}


// ---------------------------------------------------------------------
// rev 49: expand the fitted b (1 x p) and V (p x p) to the full e(b)
// stripe (K terms incl. base and omitted levels, zero there). pos holds
// the positions of the fitted terms; bn == "" or Vn == "" -> zeros.
// Vout == "" -> no V is written.
// ---------------------------------------------------------------------
void stcoxem_expandbv(string scalar bn, string scalar Vn, string scalar pos,
        real scalar K, string scalar bout, string scalar Vout)
{
    real rowvector ix, bF, b
    real matrix VF, V
    real scalar j, k

    ix = strtoreal(tokens(pos))
    bF = J(1, K, 0)
    if (bn != "") {
        b = st_matrix(bn)
        for (j=1; j<=cols(ix); j++) bF[ix[j]] = b[j]
    }
    st_matrix(bout, bF)
    if (Vout != "") {
        VF = J(K, K, 0)
        if (Vn != "") {
            V = st_matrix(Vn)
            for (j=1; j<=cols(ix); j++) {
                for (k=1; k<=cols(ix); k++) VF[ix[j], ix[k]] = V[j, k]
            }
        }
        st_matrix(Vout, VF)
    }
}


real vector stcoxem_group(real vector id, G)
{
    transmorphic A
    real scalar r, key, n
    real vector gidx

    A = asarray_create("real", 1)
    n = rows(id)
    gidx = J(n,1,.)
    G = 0
    for (r=1; r<=n; r++) {
        key = id[r]
        if (!asarray_contains(A,key)) {
            G++
            asarray(A, key, G)
        };
        gidx[r] = asarray(A, key)
    }
    return(gidx)
}

// ---------------------------------------------------------------------
// rev 48: are the shared() groups nested within the vce clusters?
// Sets the caller's locals nestok (1/0) and nclus (number of clusters).
// ---------------------------------------------------------------------
void stcoxem_nested(string scalar sharedvar, string scalar clusvar,
        string scalar touse)
{
    real vector g, c, cmap
    real scalar G, Gc, r, ok

    g = stcoxem_group(st_data(., sharedvar, touse), G)
    c = stcoxem_group(st_data(., clusvar, touse), Gc)
    cmap = J(G,1,.)
    ok = 1
    for (r=1; r<=rows(g); r++) {
        if (cmap[g[r]] == .) cmap[g[r]] = c[r]
        else if (cmap[g[r]] != c[r]) {
            ok = 0
            break
        }
    }
    st_local("nestok", strofreal(ok))
    st_local("nclus", strofreal(Gc))
}


// ---------------------------------------------------------------------
// Stable log(1+x) and stable log[Gamma(a+d)/Gamma(a)].
// ---------------------------------------------------------------------
real vector stcoxem_log1p(real vector x)
{
    real vector out, sm
    sm  = abs(x) :< 1e-4
    out = ln(1 :+ x)
    out = sm :* (x :- x:^2:/2 :+ x:^3:/3 :- x:^4:/4) :+ (1 :- sm) :* out
    return(out)
}

real vector stcoxem_lgammaratio(real vector d, real scalar a)
{
    // rev 47a: the asymptotic series (needed for large a, where
    // lngamma(a+d)-lngamma(a) cancels) assumes d small relative to a; with
    // recurrent events d can be large, so use it only where d^2 < a/1000.
    real vector asym, sel
    if (a > 1e4) {
        asym = d:*ln(a) :+ (d:*(d:-1)):/(2*a) :- (d:*(d:-1):*(2*d:-1)):/(12*a^2)
        sel  = (d:^2 :< a/1000)
        return( sel:*asym :+ (1:-sel):*(lngamma(a :+ d) :- lngamma(a)) )
    }
    else {
        return( lngamma(a :+ d) :- lngamma(a) )
    };
}

real scalar stcoxem_theta_step(real vector w1, real vector w2,
        real vector Wi, real scalar theta0, real scalar maxit, real scalar tol)
{
    real scalar theta, g, gp, mbar, it, step, a

    mbar = sum(Wi :* (w1 - w2)) / sum(Wi)
    theta = theta0
    for (it=1; it<=maxit; it++) {
        a = 1/theta
        if (a > 1e4) {
            g  = -1/(2*a) - 1/(12*a^2) + 1/(120*a^4) - 1 + mbar
            gp = -0.5 - 1/(6*a) + 1/(30*a^3)
        }
        else {
            g  = ln(theta) + digamma(a) - 1 + mbar
            gp = a - trigamma(a)*a^2
        };
        step = g/gp
        while (theta - step <= 1e-8 | theta - step > 1e6) step = step/2
        theta = theta - step
        if (abs(step) < tol) break;
    }
    return(theta)
}


// ---------------------------------------------------------------------
// Aitken (Delta-squared) EM acceleration (accelerate(aitken)). Factored
// into its own flat function (rev 12) to avoid a Mata compiler "illegal
// arglist" load error from deep nesting. Untyped arguments are aliases
// written back to the caller. After two consecutive raw EM steps,
// extrapolate to the implied fixed point (guarded against tiny
// denominators and jumps > 20x the raw step), then reset the history.
// jumped (rev 21) = 1 if theta/beta were overwritten this call, so the
// caller can resync Hi/Di/H0diff/offset to the post-jump parameters.
// ---------------------------------------------------------------------
void stcoxem_aitken_update(real scalar useaitken, real scalar p,
        theta, beta, theta_h1, theta_h2, beta_h1, beta_h2,
        n_theta_att, n_theta_acc, n_beta_att, n_beta_acc, jumped)
{
    real scalar d1t, d2t, denomt, theta_acc, kk2, cand
    real vector d1b, d2b, denomb, beta_acc
    real scalar theta_before
    real vector beta_before

    jumped = 0
    if (!useaitken) return
    theta_before = theta
    if (p>0) beta_before = beta

    if (theta_h2 < .) {
        n_theta_att = n_theta_att + 1
        d1t    = theta - theta_h1
        d2t    = theta_h1 - theta_h2
        denomt = d1t - d2t
        // rev 55: only a contracting sequence (|d1/d2| < 1)
        if (abs(denomt) > 1e-10 & abs(d1t) < abs(d2t)) {
            theta_acc = theta - d1t^2/denomt
            if (theta_acc > 1e-6 & theta_acc < 1e6 & abs(theta_acc-theta) < 20*abs(d1t)+1e-8) {
                theta = theta_acc
                n_theta_acc = n_theta_acc + 1
            }
        }
        if (p>0) {
            d1b = beta - beta_h1
            d2b = beta_h1 - beta_h2
            denomb = d1b - d2b
            beta_acc = beta
            for (kk2=1; kk2<=p; kk2++) {
                n_beta_att = n_beta_att + 1
                if (abs(denomb[kk2]) > 1e-10 & abs(d1b[kk2]) < abs(d2b[kk2])) {
                    cand = beta[kk2] - d1b[kk2]^2/denomb[kk2]
                    if (abs(cand-beta[kk2]) < 20*abs(d1b[kk2])+1e-8) {
                        beta_acc[kk2] = cand
                        n_beta_acc = n_beta_acc + 1
                    }
                }
            }
            beta = beta_acc
        }
        theta_h1 = .
        theta_h2 = .
        if (p>0) {
            beta_h1 = J(p,1,.)
            beta_h2 = J(p,1,.)
        }
    }
    else {
        theta_h2 = theta_h1
        theta_h1 = theta
        if (p>0) {
            beta_h2 = beta_h1
            beta_h1 = beta
        }
    }

    if (theta != theta_before) jumped = 1
    if (p>0) {
        if (any(beta :!= beta_before)) jumped = 1
    }
}


// ---------------------------------------------------------------------
// SQUAREM (Varadhan & Roland 2008) EM acceleration -- default since
// rev 25. Same 3-point raw-iterate history and reset-after-jump cadence
// as Aitken, but extrapolates the full (theta,beta) vector jointly with
// one scalar step length: with x0,x1,x2 the three iterates,
//     r = x1-x0, v = (x2-x1)-r, alpha = -sqrt(r'r / v'v),
//     x_new = x0 - 2*alpha*r + alpha^2*v.
// alpha is floored at -1 (never weaker than plain EM) and capped at -100;
// jumps larger than 20x the raw step, or with theta outside (1e-6,1e6),
// are rejected. Plain (unscaled) SQUAREM does not account for beta and
// theta being on different scales.
// ---------------------------------------------------------------------
void stcoxem_squarem_update(real scalar accelmode, real scalar p,
        theta, beta, theta_h1, theta_h2, beta_h1, beta_h2,
        n_att, n_acc, jumped)
{
    real scalar theta_before, sqnorm_r, sqnorm_v, alpha, rawnorm, jumpnorm, theta_new
    real vector beta_before, x0v, x1v, x2v, rvec, vvec, xnewv, beta_new

    jumped = 0
    if (accelmode != 2) return
    theta_before = theta
    if (p>0) beta_before = beta

    if (theta_h2 < .) {
        n_att = n_att + 1

        x0v = (p>0 ? (theta_h2 \ beta_h2) : theta_h2)
        x1v = (p>0 ? (theta_h1 \ beta_h1) : theta_h1)
        x2v = (p>0 ? (theta    \ beta   ) : theta)

        rvec = x1v - x0v
        vvec = (x2v - x1v) - rvec
        sqnorm_r = sum(rvec:^2)
        sqnorm_v = sum(vvec:^2)

        if (sqnorm_v > 1e-14) {
            alpha = -sqrt(sqnorm_r/sqnorm_v)
            if (alpha > -1)   alpha = -1
            if (alpha < -100) alpha = -100

            xnewv = x0v - 2*alpha*rvec + alpha^2*vvec
            theta_new = xnewv[1]
            rawnorm  = sqrt(sqnorm_r)
            jumpnorm = sqrt(sum((xnewv-x2v):^2))

            if (theta_new > 1e-6 & theta_new < 1e6 & jumpnorm < 20*rawnorm+1e-8) {
                theta = theta_new
                if (p>0) {
                    beta_new = xnewv[2::(p+1)]
                    beta = beta_new
                }
                n_acc = n_acc + 1
            }
        }

        theta_h1 = .
        theta_h2 = .
        if (p>0) {
            beta_h1 = J(p,1,.)
            beta_h2 = J(p,1,.)
        }
    }
    else {
        theta_h2 = theta_h1
        theta_h1 = theta
        if (p>0) {
            beta_h2 = beta_h1
            beta_h1 = beta
        }
    }

    if (theta != theta_before) jumped = 1
    if (p>0) {
        if (any(beta :!= beta_before)) jumped = 1
    }
}


// ---------------------------------------------------------------------
// Step-function lookup: given breakpoints bvals (sorted ascending) with
// cumulative values bcum, return for each query point (supplied already
// SORTED via qorder/qvals_sorted) the cumulative value at the largest
// breakpoint <= the query value (0 if none). O(n+m) via merge.
// ---------------------------------------------------------------------
real vector stcoxem_steplookup(real vector qorder, real vector qvals_sorted,
        real vector bvals, real vector bcum)
{
    real scalar n, m, k, bp
    real vector out

    n = rows(qorder)
    m = rows(bvals)
    out = J(n,1,0)
    bp = 0
    for (k=1; k<=n; k++) {
        while (bp < m) {
            if (bvals[bp+1] <= qvals_sorted[k]) bp++
            else break;
        }
        out[qorder[k]] = (bp>0 ? bcum[bp] : 0)
    }
    return(out)
}


// ---------------------------------------------------------------------
// One O(n) risk-set sweep over a SINGLE stratum's rows (entry/exit
// two-pointer algorithm; left truncation via t0; Breslow ties). Given
// eta = X*beta + offset, returns:
//   score, info    : Cox partial-likelihood score (p x1), information (p x p)
//   pl_ll          : this stratum's partial log-likelihood
//   rowH0diff      : n x1, [H0(t_i) - H0(t0_i)] per row
//   rowDbar1       : n x p, [cumDbar1(t_i) - cumDbar1(t0_i)], where
//                    cumDbar1(t) = running sum over event times of
//                    [D_tau/S0(tau)]*Sbar1(tau) (delta-method Louis/robust)
//   rowSbar1       : n x p, risk-set-average covariate at record i's own
//                    event time (0 for non-event rows) (Lin-Wei residual)
// rowDbar1/rowSbar1 only computed when needsens!=0 (and p>0).
// ---------------------------------------------------------------------
void stcoxem_stratum_pass(
        real matrix X, real vector t0, real vector t, real vector d,
        real vector wt, real vector eta, real scalar needsens,
        score, info, real scalar pl_ll, rowH0diff, rowDbar1, rowSbar1)
{
    real scalar n, p, i, j, k, m, ep, tau, devents, S0
    real vector r, order_entry, order_exit, S1, xevents, Sbar1
    real matrix S2, xx
    real vector evtimes, incr, cumincr
    real matrix incrV, cumincrV
    real scalar nev, etaevents, idx, maxev, mm, jj
    real scalar dosens, nrisk

    n = rows(t)
    p = cols(X)
    r = wt :* exp(eta)
    dosens = (needsens != 0 & p > 0)

    // rev 47a: Mata's order() breaks ties randomly (different permutation on
    // every call), which changed the summation order of S0/S1/S2 between
    // identical calls -> row number as tie-breaker makes the sweep
    // deterministic.
    order_entry = order((t0, (1::n)), (1,2))
    order_exit  = order((t,  (1::n)), (1,2))
    nrisk = 0

    score = J(p,1,0)
    info  = J(p,p,0)
    pl_ll = 0

    S0 = 0
    S1 = J(p,1,0)
    S2 = J(p,p,0)

    maxev   = n
    evtimes = J(maxev,1,.)
    incr    = J(maxev,1,.)
    if (dosens) {
        incrV    = J(maxev,p,.)
        rowSbar1 = J(n,p,0)
    };
    nev     = 0

    ep = 1
    i  = 1
    while (i <= n) {
        idx = order_exit[i]
        tau = t[idx]

        // ---- Step A: add newly at-risk rows (t0 < tau) ----
        while (ep <= n) {
            if (t0[order_entry[ep]] < tau) {
                j = order_entry[ep]
                S0 = S0 + r[j]
                S1 = S1 + r[j]*X[j,.]'
                xx = X[j,.]'*X[j,.]
                S2 = S2 + r[j]*xx
                nrisk++
                ep++
            }
            else break;
        }

        // ---- Step B: gather the block of rows exiting exactly at tau ----
        k = i
        devents   = 0
        xevents   = J(p,1,0)
        etaevents = 0
        while (k <= n) {
            if (t[order_exit[k]] == tau) {
                j = order_exit[k]
                if (d[j] == 1) {
                    devents   = devents + wt[j]
                    xevents   = xevents + wt[j]*X[j,.]'
                    etaevents = etaevents + wt[j]*eta[j]
                };
                k++
            }
            else break;
        }

        // ---- Step C: Breslow increment + score/info/ll ----
        if (devents > 0) {
            nev++
            evtimes[nev] = tau
            incr[nev]    = devents / S0
            Sbar1 = S1 / S0
            if (dosens) {
                incrV[nev,.] = incr[nev] * Sbar1'
                for (mm=i; mm<k; mm++) {
                    jj = order_exit[mm]
                    if (d[jj] == 1) rowSbar1[jj,.] = Sbar1';
                }
            };
            score = score + xevents - devents*Sbar1
            info  = info  + devents*( S2/S0 - Sbar1*Sbar1' )
            pl_ll = pl_ll + etaevents - devents*ln(S0)
        };

        // ---- Step D: remove the block that exits at tau ----
        for (m=i; m<k; m++) {
            j = order_exit[m]
            S0 = S0 - r[j]
            S1 = S1 - r[j]*X[j,.]'
            xx = X[j,.]'*X[j,.]
            S2 = S2 - r[j]*xx
            nrisk--
        }
        // rev 47a: risk set empty -> reset to exact zero (drop the rounding
        // residue of the add-then-subtract running sums)
        if (nrisk == 0) {
            S0 = 0
            S1 = J(p,1,0)
            S2 = J(p,p,0)
        };

        i = k
    }

    evtimes = evtimes[1::nev]
    incr    = incr[1::nev]
    cumincr = runningsum(incr)

    rowH0diff = stcoxem_steplookup(order_exit, t[order_exit], evtimes, cumincr) -
                stcoxem_steplookup(order_entry, t0[order_entry], evtimes, cumincr)

    if (dosens) {
        incrV    = incrV[1::nev,.]
        // runningsum() only accepts a vector: cumulate column by column
        cumincrV = J(nev,p,.)
        real scalar col
        for (col=1; col<=p; col++) {
            cumincrV[.,col] = runningsum(incrV[.,col])
        }
        rowDbar1 = J(n,p,0)
        for (col=1; col<=p; col++) {
            rowDbar1[.,col] = stcoxem_steplookup(order_exit, t[order_exit], evtimes, cumincrV[.,col]) -
                              stcoxem_steplookup(order_entry, t0[order_entry], evtimes, cumincrV[.,col])
        }
    }
    else {
        rowDbar1 = J(n,0,0)
        rowSbar1 = J(n,0,0)
    };
}


// ---------------------------------------------------------------------
// Core Klein-EM fit (rev 17): shared by the main fit, every bootstrap
// replicate and every profile-likelihood refit. Typed arguments are
// read-only inputs; untyped (aliased) arguments are outputs. Pass
// dolog/dodots = 0/0 for silent use. fixtheta=1 holds theta at theta0
// (profile-likelihood evaluation: re-maximize over beta only).
// ---------------------------------------------------------------------
void stcoxem_em_core(
        real vector t0, real vector t, real vector d, real vector wt,
        real vector gidx, real scalar G, real vector strat, real scalar S,
        real matrix X, real scalar p, real scalar N,
        real scalar maxiter, real scalar coxiter, real scalar tol,
        real scalar ltol, real scalar accelmode, real scalar fixtheta,
        real scalar dolog, real scalar dodots,
        real vector beta0, real scalar theta0,
        beta, theta, iter, converged, H0diff_all, Hi, Di, ll, offset,
        n_theta_att, n_theta_acc, n_beta_att, n_beta_acc, n_sq_att, n_sq_acc)
{
    real vector Wi, Wsum, Wcnt
    real scalar r, s, it2
    real vector betaold, u1row, w1, w2
    real matrix Xs
    real vector t0s, ts, ds, wts, etas, scoretot
    real matrix infotot
    real scalar thetaold, llold, Dtot
    transmorphic strat_rows

    Dtot = sum(wt :* d)
    Wsum = J(G,1,0)
    Wcnt = J(G,1,0)
    for (r=1; r<=N; r++) {
        Wsum[gidx[r]] = Wsum[gidx[r]] + wt[r]
        Wcnt[gidx[r]] = Wcnt[gidx[r]] + 1
    }
    Wi = Wsum :/ Wcnt

    strat_rows = asarray_create("real", 1)
    for (s=1; s<=S; s++) {
        asarray(strat_rows, s, selectindex(strat:==s))
    }

    beta  = beta0
    theta = theta0
    llold = .
    ll    = .
    converged = 0

    real scalar theta_h1, theta_h2
    real vector beta_h1, beta_h2
    theta_h1 = .
    theta_h2 = .
    if (p>0) {
        beta_h1 = J(p,1,.)
        beta_h2 = J(p,1,.)
    }

    n_theta_att = 0; n_theta_acc = 0; n_beta_att = 0; n_beta_acc = 0
    n_sq_att = 0; n_sq_acc = 0

    if (dodots) printf("  EM iterations: ");

    for (iter=1; iter<=maxiter; iter++) {

        // ---- M-step uses the offset from the PREVIOUS E-step (0 on iter 1) ----
        if (iter == 1) offset = J(N,1,0);

        // ---- M-step, part 1: coxiter Newton steps for beta only ----
        // (rev 47a: betaold = beta at the START of the EM cycle, so the
        //  convergence check sees the whole cycle's change for coxiter>1)
        if (p>0) betaold = beta
        real scalar unident
        unident = 0
        for (it2=1; it2<=coxiter; it2++) {
            real vector etafull
            if (p>0) etafull = X*beta :+ offset
            else     etafull = offset;

            scoretot = J(p,1,0)
            infotot  = J(p,p,0)
            for (s=1; s<=S; s++) {
                real vector ridx
                ridx = asarray(strat_rows, s)
                if (rows(ridx)==0) continue;
                Xs  = (p>0 ? X[ridx,.] : J(rows(ridx),0,.))
                t0s = t0[ridx]; ts = t[ridx]; ds = d[ridx]; wts = wt[ridx]
                etas = etafull[ridx]

                real vector score_s, rowH0diff_s2
                real matrix info_s, rowDbar1_s2, rowSbar1_s2
                real scalar pl_ll_s
                stcoxem_stratum_pass(Xs, t0s, ts, ds, wts, etas, 0,
                    score_s, info_s, pl_ll_s, rowH0diff_s2, rowDbar1_s2, rowSbar1_s2)

                if (p>0) {
                    scoretot = scoretot + score_s
                    infotot  = infotot  + info_s
                };
            }
            // rev 47a: invsym() zeroes a row/column it cannot pivot on
            // (unidentified coefficient); flag it instead of silently
            // freezing that coefficient (kept flat: deep nesting breaks
            // Mata's compiler)
            if (p>0) {
                real matrix Iinv
                Iinv = invsym(infotot)
                unident = any(diagonal(Iinv) :== 0)
                if (!unident) beta = beta + Iinv * scoretot
            };
            if (unident) break;
        }
        if (unident) {
            beta = J(p,1,.)
            converged = 0
            break
        };

        // ---- M-step, part 2: one more pass at the final beta for a
        //      baseline hazard (and pl_ll, info) consistent with it ----
        real vector etafinal2
        real scalar pl_total
        if (p>0) etafinal2 = X*beta :+ offset
        else     etafinal2 = offset;
        H0diff_all = J(N,1,0)
        scoretot = J(p,1,0)
        infotot  = J(p,p,0)
        pl_total = 0
        for (s=1; s<=S; s++) {
            real vector ridx4
            ridx4 = asarray(strat_rows, s)
            if (rows(ridx4)==0) continue;
            Xs  = (p>0 ? X[ridx4,.] : J(rows(ridx4),0,.))
            t0s = t0[ridx4]; ts = t[ridx4]; ds = d[ridx4]; wts = wt[ridx4]
            etas = etafinal2[ridx4]

            real vector score_s3, rowH0diff_s3
            real matrix info_s3, rowDbar1_s3, rowSbar1_s3
            real scalar pl_ll_s3
            stcoxem_stratum_pass(Xs, t0s, ts, ds, wts, etas, 0,
                score_s3, info_s3, pl_ll_s3, rowH0diff_s3, rowDbar1_s3, rowSbar1_s3)

            if (p>0) {
                scoretot = scoretot + score_s3
                infotot  = infotot  + info_s3
            };
            H0diff_all[ridx4] = rowH0diff_s3
            pl_total = pl_total + pl_ll_s3
        }

        // ---- Hi, Di: cluster-level (weighted) cumulative hazard / events ----
        real vector Hrow
        if (p>0) Hrow = H0diff_all :* exp(X*beta)
        else     Hrow = H0diff_all;
        Hi = J(G,1,0)
        Di = J(G,1,0)
        for (r=1; r<=N; r++) {
            Hi[gidx[r]] = Hi[gidx[r]] + (wt[r]/Wi[gidx[r]])*Hrow[r]
            Di[gidx[r]] = Di[gidx[r]] + (wt[r]/Wi[gidx[r]])*d[r]
        }

        // ---- E-step ----
        w1 = (1/theta :+ Di) :/ (1/theta :+ Hi)
        w2 = digamma(1/theta :+ Di) :- ln(1/theta :+ Hi)
        u1row  = w1[gidx]
        offset = ln(u1row)

        // ---- theta M-step (skipped when fixtheta) ----
        thetaold = theta
        if (!fixtheta) theta = stcoxem_theta_step(w1, w2, Wi, theta, 100, 1e-10)

        // ---- EM acceleration: 2=squarem (default), 1=aitken, 0=none ----
        real scalar jumped
        if (accelmode == 1) {
            stcoxem_aitken_update(1, p, theta, beta, theta_h1, theta_h2, beta_h1, beta_h2,
                n_theta_att, n_theta_acc, n_beta_att, n_beta_acc, jumped)
        }
        else if (accelmode == 2) {
            stcoxem_squarem_update(accelmode, p, theta, beta, theta_h1, theta_h2, beta_h1, beta_h2,
                n_sq_att, n_sq_acc, jumped)
        }
        else {
            jumped = 0
        }

        // ---- resync after an accelerated jump (rev 21): Hi, Di, H0diff_all
        //      and offset were computed from the PRE-jump beta/theta ----
        if (jumped) {
            real vector etajump, Hrowj
            real scalar pl_totalj
            if (p>0) etajump = X*beta :+ offset
            else     etajump = offset
            H0diff_all = J(N,1,0)
            pl_totalj = 0
            for (s=1; s<=S; s++) {
                real vector ridxj
                ridxj = asarray(strat_rows, s)
                if (rows(ridxj)==0) continue
                Xs  = (p>0 ? X[ridxj,.] : J(rows(ridxj),0,.))
                t0s = t0[ridxj]; ts = t[ridxj]; ds = d[ridxj]; wts = wt[ridxj]
                etas = etajump[ridxj]

                real vector score_sj, rowH0diff_sj
                real matrix info_sj, rowDbar1_sj, rowSbar1_sj
                real scalar pl_ll_sj
                stcoxem_stratum_pass(Xs, t0s, ts, ds, wts, etas, 0,
                    score_sj, info_sj, pl_ll_sj, rowH0diff_sj, rowDbar1_sj, rowSbar1_sj)
                H0diff_all[ridxj] = rowH0diff_sj
                pl_totalj = pl_totalj + pl_ll_sj
            }
            pl_total = pl_totalj

            if (p>0) Hrowj = H0diff_all :* exp(X*beta)
            else     Hrowj = H0diff_all
            Hi = J(G,1,0)
            Di = J(G,1,0)
            for (r=1; r<=N; r++) {
                Hi[gidx[r]] = Hi[gidx[r]] + (wt[r]/Wi[gidx[r]])*Hrowj[r]
                Di[gidx[r]] = Di[gidx[r]] + (wt[r]/Wi[gidx[r]])*d[r]
            }

            w1 = (1/theta :+ Di) :/ (1/theta :+ Hi)
            w2 = digamma(1/theta :+ Di) :- ln(1/theta :+ Hi)
            u1row  = w1[gidx]
            offset = ln(u1row)
        }

        // ---- marginal log-likelihood monitor. REV 31: the offset-partial-
        //      likelihood trick (offset=ln(w1)) leaves a spurious
        //      sum(Di*ln(w1)) term in pl_total's VALUE (harmless for beta,
        //      but it corrupts curvature in theta), so it is subtracted.
        //      Rev 50 adds the constant D = sum(wt*d), so ll is on stcox's
        //      scale (lrtest, estat ic, LR test of theta=0). ----
        real scalar a
        real vector lgratio
        a = 1/theta
        lgratio = stcoxem_lgammaratio(Di, a)
        // rev 40: cluster-level terms scaled by Wi (cluster replication
        // weight); Hi/Di are per-copy (weight-normalized) quantities.
        // rev 50: + Dtot = sum(wt*d) puts ll on stcox's scale. At the
        // Breslow point the expression before Dtot equals the marginal
        // log likelihood minus sum_k d_k ln d_k; stcox, shared() reports
        // a scale that equals the Breslow partial likelihood as theta->0,
        // where the marginal ll tends to PL + sum d ln d - D. Hence
        // ll_stcox = (expression) + D: a data-only constant, valid with
        // ties, strata, late entry; under pweights D is weighted (k copies).
        ll = pl_total - sum(Wi :* Di :* ln(w1)) + sum( Wi :* ( lgratio :- Di:*ln(a)
                              :- (Di :+ a) :* stcoxem_log1p(Hi:/a) ) ) + Dtot

        if (dolog) {
            printf("  EM iter %3.0f:  logL = %12.4f   theta = %8.5f\n", iter, ll, theta)
        }
        else if (dodots) {
            printf(".")
        };

        // ---- convergence (rev 47a, Stata convention): parameter change
        //      |new-old|/(|old|+1) below tol AND relative change of the
        //      log likelihood |ll-llold|/(|llold|+1) below ltol ----
        real scalar parchange, llchange
        parchange = abs(theta-thetaold)/(abs(thetaold)+1)
        if (p>0) parchange = max((parchange, max(abs(beta-betaold):/(abs(betaold):+1))));
        llchange = (llold < . ? abs(ll-llold)/(abs(llold)+1) : .)
        if (iter>1 & parchange < tol & llchange < ltol) {
            converged = 1
            break
        };
        llold = ll
    }
    if (dodots) printf("\n");
}


// ---------------------------------------------------------------------
// rev 50: comparison model for the LR test of theta=0 -- plain Cox
// (Breslow ties, same sample/strata/weights/X) by Newton-Raphson with step
// halving. stcox_fr fits it with -cox ..., nocoef norefine-; same model,
// validated against -stcox- in the rev 50 tests. Returns the log partial
// likelihood (stcox's scale); conv_out = 1 if converged.
// ---------------------------------------------------------------------
real scalar stcoxem_cox_pass(
        real vector t0, real vector t, real vector d, real vector wt,
        transmorphic strat_rows, real scalar S, real matrix X, real scalar p,
        real vector beta, score, info)
{
    real scalar s, llv, pl_s
    real vector ridx, eta, score_s, rowH0_s
    real matrix info_s, rD_s, rS_s

    if (p>0) eta = X*beta
    else     eta = J(rows(t),1,0);
    llv   = 0
    score = J(p,1,0)
    info  = J(p,p,0)
    for (s=1; s<=S; s++) {
        ridx = asarray(strat_rows, s)
        if (rows(ridx)==0) continue;
        stcoxem_stratum_pass((p>0 ? X[ridx,.] : J(rows(ridx),0,.)),
            t0[ridx], t[ridx], d[ridx], wt[ridx], eta[ridx], 0,
            score_s, info_s, pl_s, rowH0_s, rD_s, rS_s)
        llv = llv + pl_s
        if (p>0) {
            score = score + score_s
            info  = info  + info_s
        };
    }
    return(llv)
}

real scalar stcoxem_cox_ll(
        real vector t0, real vector t, real vector d, real vector wt,
        real vector strat, real scalar S, real matrix X, real scalar p,
        real scalar tol, real scalar ltol, real scalar maxit, conv_out)
{
    real scalar s, it, h, llv, llnew, parchange, llchange
    real vector beta, bnew, step, score, scnew
    real matrix info, infnew
    transmorphic strat_rows

    strat_rows = asarray_create("real", 1)
    for (s=1; s<=S; s++) {
        asarray(strat_rows, s, selectindex(strat:==s))
    }
    beta = J(p,1,0)
    llv  = stcoxem_cox_pass(t0, t, d, wt, strat_rows, S, X, p, beta, score, info)
    conv_out = 1
    if (p==0) return(llv);

    conv_out = 0
    for (it=1; it<=maxit; it++) {
        step = cholsolve(info, score)
        if (hasmissing(step)) step = invsym(info)*score;
        // step halving until the log likelihood does not decrease
        for (h=0; h<=30; h++) {
            bnew  = beta + step/(2^h)
            llnew = stcoxem_cox_pass(t0, t, d, wt, strat_rows, S, X, p, bnew, scnew, infnew)
            if (llnew < .) {
                if (llnew >= llv - 1e-12*(abs(llv)+1)) break;
            };
        }
        if (llnew >= .) return(.);
        parchange = max(abs(bnew-beta):/(abs(beta):+1))
        llchange  = abs(llnew-llv)/(abs(llv)+1)
        beta  = bnew
        llv   = llnew
        score = scnew
        info  = infnew
        if (parchange < tol & llchange < ltol) {
            conv_out = 1
            break
        };
    }
    return(llv)
}


// ---------------------------------------------------------------------
// Profile log-likelihood at a FIXED theta (rev 26): one call to
// stcoxem_em_core with fixtheta=1. beta_warm is the warm start (chained by
// callers); beta_out/conv_out are aliased outputs.
// ---------------------------------------------------------------------
real scalar stcoxem_profile_ll(
        real vector t0, real vector t, real vector d, real vector wt,
        real vector gidx, real scalar G, real vector strat, real scalar S,
        real matrix X, real scalar p, real scalar N,
        real scalar maxiter, real scalar coxiter, real scalar tol,
        real scalar ltol, real scalar accelmode, real scalar theta_fixed,
        real vector beta_warm,
        beta_out, conv_out)
{
    real scalar theta_r, iter_r, converged_r, ll_r
    real vector H0diff_r, Hi_r, Di_r, offset_r
    real scalar nt_att_r, nt_acc_r, nb_att_r, nb_acc_r, nsq_att_r, nsq_acc_r

    stcoxem_em_core(t0, t, d, wt, gidx, G, strat, S, X, p, N,
        maxiter, coxiter, tol, ltol, accelmode, 1, 0, 0,
        beta_warm, theta_fixed,
        beta_out, theta_r, iter_r, converged_r, H0diff_r, Hi_r, Di_r, ll_r, offset_r,
        nt_att_r, nt_acc_r, nb_att_r, nb_acc_r, nsq_att_r, nsq_acc_r)

    conv_out = converged_r
    return(ll_r)
}


// ---------------------------------------------------------------------
// One profile-likelihood CI bound for theta (rev 26): bracket-then-bisect
// on theta with geometric steps/geometric-mean bisection. direction<0 =
// lower bound, >0 = upper bound. Returns theta_floor or 1e6 if the search
// leaves a sane range. Each evaluation warm-starts from the previous
// converged beta.
// ---------------------------------------------------------------------
real scalar stcoxem_profile_bound(
        real vector t0, real vector t, real vector d, real vector wt,
        real vector gidx, real scalar G, real vector strat, real scalar S,
        real matrix X, real scalar p, real scalar N,
        real scalar maxiter, real scalar coxiter, real scalar tol,
        real scalar ltol, real scalar accelmode, real scalar theta_hat,
        real vector beta_hat,
        real scalar target_ll, real scalar direction, real scalar theta_floor,
        real scalar showprog, n_nonconverged)
{
    real scalar factor, theta_in, theta_out, ll_out, ll_mid, mid, conv, step
    real vector beta_cur, beta_new

    factor   = (direction < 0 ? (1/1.6) : 1.6)
    theta_in = theta_hat
    theta_out = theta_hat
    beta_cur = beta_hat

    step = 0
    while (step < 60) {
        theta_out = theta_out * factor
        if (theta_out < theta_floor | theta_out > 1e6) {
            return(direction < 0 ? theta_floor : 1e6)
        }
        ll_out = stcoxem_profile_ll(t0,t,d,wt,gidx,G,strat,S,X,p,N,
            maxiter,coxiter,tol,ltol,accelmode,theta_out,beta_cur, beta_new,conv)
        beta_cur = beta_new
        if (showprog) printf(".")
        if (!conv) n_nonconverged = n_nonconverged + 1
        if (ll_out < target_ll) break
        theta_in = theta_out
        step++
    }

    step = 0
    while (step < 60) {
        mid = sqrt(theta_in * theta_out)
        ll_mid = stcoxem_profile_ll(t0,t,d,wt,gidx,G,strat,S,X,p,N,
            maxiter,coxiter,tol,ltol,accelmode,mid,beta_cur, beta_new,conv)
        beta_cur = beta_new
        if (showprog) printf(".")
        if (!conv) n_nonconverged = n_nonconverged + 1
        if (ll_mid >= target_ll) theta_in = mid
        else                     theta_out = mid
        if (abs(theta_in - theta_out) < 1e-5*theta_hat + 1e-8) break
        step++
    }
    return(theta_in)
}


// ---------------------------------------------------------------------
// Full EM driver.
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// rev 51: ln u_i at fixed beta (through xb) and theta -- the EM E-step
// iterated to its fixed point: offset -> Breslow baseline (per stratum) ->
// Hi, Di -> w1 = (a+Di)/(a+Hi) -> offset = ln w1. Same Hi/Di/Wi
// conventions as stcoxem_em_core. xb is centred (u is invariant to it).
// Writes ln u to outvar (touse rows); sets locals lnu_conv, lnu_iter.
// ---------------------------------------------------------------------
void stcoxem_lnu(string scalar xbvar, string scalar sharedvar,
        string scalar stratvar, string scalar wtvar, string scalar touse,
        string scalar outvar)
{
    real vector t0, t, d, eta, gidx, strat, wt, Wsum, Wcnt, Wi
    real vector offset, offold, H0diff_all, Hrow, Hi, Di, w1, ridx
    real vector score_s, rowH0diff_s
    real matrix info_s, rowDbar1_s, rowSbar1_s
    real scalar N, G, S, s, r, it, maxit, conv, theta, a, pl_ll_s
    transmorphic strat_rows

    theta = st_numscalar("e(theta)")
    t0    = st_data(., "_t0", touse)
    t     = st_data(., "_t",  touse)
    d     = st_data(., "_d",  touse)
    eta   = st_data(., xbvar, touse)
    strat = st_data(., stratvar, touse)
    wt    = st_data(., wtvar, touse)
    gidx  = stcoxem_group(st_data(., sharedvar, touse), G)
    N = rows(t)
    S = max(strat)
    eta = eta :- mean(eta)

    Wsum = J(G,1,0)
    Wcnt = J(G,1,0)
    for (r=1; r<=N; r++) {
        Wsum[gidx[r]] = Wsum[gidx[r]] + wt[r]
        Wcnt[gidx[r]] = Wcnt[gidx[r]] + 1
    }
    Wi = Wsum :/ Wcnt
    strat_rows = asarray_create("real", 1)
    for (s=1; s<=S; s++) {
        asarray(strat_rows, s, selectindex(strat:==s))
    }

    a = 1/theta
    offset = J(N,1,0)
    maxit = 10000
    conv = 0
    for (it=1; it<=maxit; it++) {
        H0diff_all = J(N,1,0)
        for (s=1; s<=S; s++) {
            ridx = asarray(strat_rows, s)
            if (rows(ridx)==0) continue
            stcoxem_stratum_pass(J(rows(ridx),0,.), t0[ridx], t[ridx], d[ridx],
                wt[ridx], eta[ridx] :+ offset[ridx], 0,
                score_s, info_s, pl_ll_s, rowH0diff_s, rowDbar1_s, rowSbar1_s)
            H0diff_all[ridx] = rowH0diff_s
        }
        Hrow = H0diff_all :* exp(eta)
        Hi = J(G,1,0)
        Di = J(G,1,0)
        for (r=1; r<=N; r++) {
            Hi[gidx[r]] = Hi[gidx[r]] + (wt[r]/Wi[gidx[r]])*Hrow[r]
            Di[gidx[r]] = Di[gidx[r]] + (wt[r]/Wi[gidx[r]])*d[r]
        }
        w1 = (a :+ Di) :/ (a :+ Hi)
        offold = offset
        offset = ln(w1[gidx])
        if (max(abs(offset :- offold)) < 1e-12) {
            conv = 1
            break
        }
    }
    st_store(., outvar, touse, offset)
    st_local("lnu_conv", strofreal(conv))
    st_local("lnu_iter", strofreal(it <= maxit ? it : maxit))
}


void stcoxem_driver(
        string scalar xvars, string scalar sharedvar, string scalar stratvar,
        string scalar wtvar, string scalar touse, string scalar subjvar,
        real scalar S, real scalar maxiter, real scalar coxiter,
        real scalar tol, real scalar ltol, real scalar dolog, real scalar dodots,
        real scalar accelmode, real scalar vcemode,
        string scalar genvar, real scalar reps, real scalar bsseed,
        real scalar usebswarm, real scalar doplci, real scalar level,
        real scalar doplse, string scalar ufvar, string scalar clusvar,
        real scalar dolrtest)
{
    real vector t0, t, d, cid, strat, wt, subj
    real matrix X
    real scalar N, p, G, s, r, iter, converged
    real vector gidx, Wi, Wsum, Wcnt
    real vector beta
    real scalar theta, se_theta, ll
    real vector Hi, Di, H0diff_all
    string vector xnames
    real matrix V, Vbeta
    real vector offset
    real matrix infotot, rowDbar1_all, rowSbar1_all
    real scalar n_theta_att, n_theta_acc, n_beta_att, n_beta_acc, n_sq_att, n_sq_acc
    transmorphic strat_rows

    t0    = st_data(., "_t0", touse)
    t     = st_data(., "_t",  touse)
    d     = st_data(., "_d",  touse)
    cid   = st_data(., sharedvar, touse)
    strat = st_data(., stratvar, touse)
    wt    = st_data(., wtvar, touse)
    subj  = st_data(., subjvar, touse)
    if (strtrim(xvars) != "") {
        xnames = tokens(xvars)
        X = st_data(., xnames, touse)
        p = cols(X)
    }
    else {
        X = J(rows(t),0,.)
        p = 0
    };
    N = rows(t)

    // rev 47a: centre the covariates (estimation-sample column means).
    // Cox partial likelihood, beta, theta, ll and all variances are
    // invariant to this shift; the Breslow baseline is then the hazard at
    // x = xbar and is rescaled to x = 0 at output (exp(-xbar*beta)).
    real rowvector xbar
    if (p>0) {
        xbar = mean(X)
        X = X :- xbar
    };

    gidx = stcoxem_group(cid, G)
    real scalar GS
    real vector sgidx
    sgidx = stcoxem_group(subj, GS)

    // ---- main fit: cold-started at beta=0, theta=0.5 ----
    stcoxem_em_core(t0, t, d, wt, gidx, G, strat, S, X, p, N,
        maxiter, coxiter, tol, ltol, accelmode, 0, dolog, dodots, J(p,1,0), 0.5,
        beta, theta, iter, converged, H0diff_all, Hi, Di, ll, offset,
        n_theta_att, n_theta_acc, n_beta_att, n_beta_acc, n_sq_att, n_sq_acc)
    if (p>0) {
        if (hasmissing(beta)) {
            errprintf("coefficients not identified during the EM (covariates" +
                " collinear, or without variation within the risk sets);" +
                " remove covariates and refit\n")
            exit(459)
        }
    }

    // ---- rev 50: comparison model (plain Cox) for the LR test of theta=0 ----
    if (dolrtest) {
        // (names llcmp*: ll_c is already the plse centre point below)
        real scalar llcmp, llcmp_conv
        llcmp = stcoxem_cox_ll(t0, t, d, wt, strat, S, X, p, 1e-10, 1e-12, 100, llcmp_conv)
        st_numscalar("r(ll_c)", llcmp)
        st_numscalar("r(ll_c_conv)", llcmp_conv)
    }

    // Aitken/SQUAREM accept/reject counts are kept in e() for anyone who
    // wants to inspect them (ereturn list), but not printed on every fit.
    st_numscalar("r(aitken_theta_att)", n_theta_att)
    st_numscalar("r(aitken_theta_acc)", n_theta_acc)
    st_numscalar("r(aitken_beta_att)", n_beta_att)
    st_numscalar("r(aitken_beta_acc)", n_beta_acc)
    st_numscalar("r(squarem_att)", n_sq_att)
    st_numscalar("r(squarem_acc)", n_sq_acc)

    // ---- plci(): profile-likelihood CI for theta (rev 26). Independent
    //      of vce(); sidesteps the Louis decomposition entirely. Each
    //      bound = bracket search + bisection, each step a full EM refit. ----
    if (doplci) {
        real scalar target_ll, theta_pl_lower, theta_pl_upper, n_pl_nonconv
        real scalar showprog_pl

        showprog_pl = (dolog | dodots)
        target_ll = ll - invchi2(1, level/100)/2
        n_pl_nonconv = 0

        if (showprog_pl) printf("  Profile-likelihood CI for theta: ")
        theta_pl_lower = stcoxem_profile_bound(t0, t, d, wt, gidx, G, strat, S,
            X, p, N, maxiter, coxiter, tol, ltol, accelmode, theta, beta,
            target_ll, -1, 1e-6, showprog_pl, n_pl_nonconv)
        theta_pl_upper = stcoxem_profile_bound(t0, t, d, wt, gidx, G, strat, S,
            X, p, N, maxiter, coxiter, tol, ltol, accelmode, theta, beta,
            target_ll, 1, 1e-6, showprog_pl, n_pl_nonconv)
        if (showprog_pl) printf("\n")

        printf("  Profile-likelihood CI for theta: [%9.6f, %9.6f]\n",
            theta_pl_lower, theta_pl_upper)
        if (n_pl_nonconv > 0) {
            printf("  (%g of the profile refits above did not converge" +
                " within maxiter() -- treat these bounds as provisional;" +
                " consider raising maxiter())\n", n_pl_nonconv)
        }
        st_numscalar("r(theta_pl_lower)", theta_pl_lower)
        st_numscalar("r(theta_pl_upper)", theta_pl_upper)
        st_numscalar("r(theta_pl_nonconv)", n_pl_nonconv)
    }

    // ---- plse (rev 28-30): numerical 2nd derivative of the profile
    //      log-likelihood at the MLE, se = 1/sqrt(-ll''). Relative steps
    //      hrel1=0.01 and hrel2=0.005 (x theta); refits use a dedicated
    //      tight tolerance tol_pl; the CENTER point is refit at tol_pl too.
    //      Reported (Richardson-extrapolated) only if both curvatures are
    //      negative and agree within 25%; otherwise missing, with the raw
    //      numbers printed. ----
    if (doplse) {
        real scalar hrel1, hrel2, h1, h2, tol_pl
        real scalar ll_p1, ll_m1, ll_p2, ll_m2
        real scalar conv_p1, conv_m1, conv_p2, conv_m2, n_prof_nonconv
        real scalar D1, D2, Drich, reldiff, se_theta_profile
        real vector beta_p1, beta_m1, beta_p2, beta_m2

        hrel1 = 0.01
        hrel2 = 0.005
        h1 = theta * hrel1
        h2 = theta * hrel2
        tol_pl = min((tol, 1e-10))
        real scalar ltol_pl
        ltol_pl = min((ltol, 1e-12))
        n_prof_nonconv = 0

        real scalar ll_c, conv_c
        real vector beta_c
        ll_c = stcoxem_profile_ll(t0,t,d,wt,gidx,G,strat,S,X,p,N,
            maxiter,coxiter,tol_pl,ltol_pl,accelmode, theta, beta, beta_c, conv_c)
        if (!conv_c) n_prof_nonconv = n_prof_nonconv + 1

        ll_p1 = stcoxem_profile_ll(t0,t,d,wt,gidx,G,strat,S,X,p,N,
            maxiter,coxiter,tol_pl,ltol_pl,accelmode, theta+h1, beta_c, beta_p1, conv_p1)
        if (!conv_p1) n_prof_nonconv = n_prof_nonconv + 1
        ll_m1 = stcoxem_profile_ll(t0,t,d,wt,gidx,G,strat,S,X,p,N,
            maxiter,coxiter,tol_pl,ltol_pl,accelmode, theta-h1, beta_c, beta_m1, conv_m1)
        if (!conv_m1) n_prof_nonconv = n_prof_nonconv + 1
        ll_p2 = stcoxem_profile_ll(t0,t,d,wt,gidx,G,strat,S,X,p,N,
            maxiter,coxiter,tol_pl,ltol_pl,accelmode, theta+h2, beta_c, beta_p2, conv_p2)
        if (!conv_p2) n_prof_nonconv = n_prof_nonconv + 1
        ll_m2 = stcoxem_profile_ll(t0,t,d,wt,gidx,G,strat,S,X,p,N,
            maxiter,coxiter,tol_pl,ltol_pl,accelmode, theta-h2, beta_c, beta_m2, conv_m2)
        if (!conv_m2) n_prof_nonconv = n_prof_nonconv + 1

        D1 = (ll_p1 - 2*ll_c + ll_m1) / (h1^2)
        D2 = (ll_p2 - 2*ll_c + ll_m2) / (h2^2)

        // se_theta_profile (Richardson-extrapolated numerical 2nd derivative
        // of the profile log-likelihood) is computed here and reported via
        // e(se_theta)/e(se_theta_profile); the intermediate diagnostics
        // (step-size agreement, curvature sign, nonconvergence count) are
        // no longer printed on every fit -- e(theta_profile_nonconv) still
        // reports how many of the 5 internal refits failed to converge.
        if (D1 < 0 & D2 < 0) {
            reldiff = abs(D2 - D1) / abs(D1)
            if (reldiff <= 0.25) {
                Drich = (4*D2 - D1) / 3
                if (Drich < 0) se_theta_profile = 1 / sqrt(-Drich)
                else           se_theta_profile = .
            }
            else {
                se_theta_profile = .
            }
        }
        else {
            se_theta_profile = .
        }
        st_numscalar("r(se_theta_profile)", se_theta_profile)
        st_numscalar("r(theta_profile_nonconv)", n_prof_nonconv)
    }

    // ---- Wi and strat_rows are local to stcoxem_em_core; rebuild here ----
    Wsum = J(G,1,0)
    Wcnt = J(G,1,0)
    for (r=1; r<=N; r++) {
        Wsum[gidx[r]] = Wsum[gidx[r]] + wt[r]
        Wcnt[gidx[r]] = Wcnt[gidx[r]] + 1
    }
    Wi = Wsum :/ Wcnt

    strat_rows = asarray_create("real", 1)
    for (s=1; s<=S; s++) {
        asarray(strat_rows, s, selectindex(strat:==s))
    }

    // ---- ONE extra pass at the converged (beta,offset), only if oim/
    //      robust was requested: infotot plus the rowDbar1/rowSbar1
    //      sensitivity outputs. ----
    if (p>0 & (vcemode==1 | vcemode==2)) {
        real vector etasens
        etasens = X*beta :+ offset
        rowDbar1_all = J(N,p,0)
        rowSbar1_all = J(N,p,0)
        infotot = J(p,p,0)
        for (s=1; s<=S; s++) {
            real vector ridx5
            ridx5 = asarray(strat_rows, s)
            if (rows(ridx5)==0) continue;
            real vector t0sn, tsn, dsn, wtsn, etasn
            t0sn = t0[ridx5]; tsn = t[ridx5]; dsn = d[ridx5]; wtsn = wt[ridx5]
            etasn = etasens[ridx5]

            real vector score_s4, rowH0diff_s4
            real matrix info_s4, rowDbar1_s4, rowSbar1_s4
            real scalar pl_ll_s4
            stcoxem_stratum_pass(X[ridx5,.], t0sn, tsn, dsn, wtsn, etasn, 1,
                score_s4, info_s4, pl_ll_s4, rowH0diff_s4, rowDbar1_s4, rowSbar1_s4)

            infotot = infotot + info_s4
            rowDbar1_all[ridx5,.] = rowDbar1_s4
            rowSbar1_all[ridx5,.] = rowSbar1_s4
        }
    };

    // ---- Covariance for beta AND theta: joint Louis (1982) EM correction.
    //      l_com = sum_i log f_gamma(u_i;theta) + PL(beta; offset=log u), so
    //      the beta/theta COMPLETE-information cross term is exactly zero.
    //      theta's own missing-information piece is exact (Gamma moments of
    //      u_i|data); beta's piece and the beta-theta cross term use a
    //      first-order (delta-method) linearization of the partial-
    //      likelihood score in u around its posterior mean, via the
    //      per-cluster sensitivity dS_beta/du_i built from rowDbar1_all/
    //      H0diff_all.
    //      KNOWN, EXPLAINED (rev 26/31): theta's oim/robust SE sits ~13-20%
    //      below -stcox-/bootstrap/frailtyEM because the decomposition
    //      treats the nonparametric Breslow baseline as known (semiparametric
    //      information deficit; Louis 1982, Gill 1989, Nielsen et al. 1992).
    //      Ctt/Mtt were re-derived three times and proven exactly correct at
    //      EM convergence (via the theta M-step stationarity identity), so
    //      this is NOT a coding bug. Use plse (default), plci or bootstrap
    //      for theta's SE. Beta's oim SE matches -stcox- to ~1.6%.
    if (p>0) {
        st_matrix("stcoxem_b", beta')
        st_matrixcolstripe("stcoxem_b", (J(p,1,""), xnames'))
    }
    else {
        st_matrix("stcoxem_b", J(1,0,.))
    };

    if (vcemode==1 | vcemode==2) {
        real scalar a2, gp2, Ctt, Mtt, ai, bi, kk
        real vector missing_i, r0
        real matrix infl_clust, Cbb, Mbb, Mbt, Ijoint, Vjoint

        a2 = 1/theta
        if (a2 > 1e4) gp2 = -0.5 - 1/(6*a2) + 1/(30*a2^3)
        else          gp2 = a2 - trigamma(a2)*a2^2;
        Ctt = -(sum(Wi)/theta^2) * gp2

        missing_i = J(G,1,.)
        for (kk=1; kk<=G; kk++) {
            ai = Di[kk] + a2
            bi = Hi[kk] + a2
            missing_i[kk] = ai/bi^2 - 2/bi + trigamma(ai)
        }
        // Mtt = Var(complete-data score | data) = sum_i Wi^2*Var(score_i|data)
        // (Wi^2, fixed rev 16; only matters under pweight).
        Mtt = (1/theta^4) * sum(Wi :* missing_i)

        if (p>0) {
            // per-record sensitivity dS_beta/du_i, summed by cluster
            // (deliberately WITHOUT the offset -- it is a derivative wrt u)
            real matrix infl_rec
            r0       = wt :* exp(X*beta)
            infl_rec = r0 :* (rowDbar1_all - X:*H0diff_all)
            infl_clust = J(G,p,0)
            for (r=1; r<=N; r++) {
                infl_clust[gidx[r],.] = infl_clust[gidx[r],.] + infl_rec[r,.]
            }

            Cbb = infotot
            Mbb = J(p,p,0)
            Mbt = J(p,1,0)
            for (kk=1; kk<=G; kk++) {
                ai = Di[kk] + a2
                bi = Hi[kk] + a2
                Mbb = Mbb + (ai/bi^2)/Wi[kk] * (infl_clust[kk,.]' * infl_clust[kk,.])
                Mbt = Mbt + (1/theta^2)*(ai/bi^2 - 1/bi) * infl_clust[kk,.]'
            }

            Ijoint = J(p+1,p+1,.)
            Ijoint[1::p, 1::p] = Cbb - Mbb
            Ijoint[1::p, p+1]  = -Mbt
            Ijoint[p+1, 1::p]  = (-Mbt)'
            Ijoint[p+1, p+1]   = Ctt - Mtt

            Vjoint = invsym(Ijoint)
            Vbeta  = Vjoint[1::p, 1::p]
            if (Vjoint[p+1,p+1] > 0) se_theta = sqrt(Vjoint[p+1,p+1])
            else                     se_theta = .;

            // theta SE ignoring beta-theta coupling, and the Schur complement
            // term that couples them; kept in e() for inspection, not printed.
            real scalar se_theta_naive, schur_term
            if (Ctt - Mtt > 0) se_theta_naive = sqrt(1/(Ctt-Mtt))
            else               se_theta_naive = .;
            schur_term = Mbt' * invsym(Cbb - Mbb) * Mbt
            st_numscalar("r(se_theta_naive)", se_theta_naive)
            st_numscalar("r(schur_term)", schur_term)

            st_matrix("stcoxem_V", Vbeta)
            st_matrixcolstripe("stcoxem_V", (J(p,1,""), xnames'))
            st_matrixrowstripe("stcoxem_V", (J(p,1,""), xnames'))
        }
        else {
            if (Ctt - Mtt > 0) se_theta = sqrt(1/(Ctt-Mtt))
            else               se_theta = .;
        };
        st_numscalar("r(se_theta)", se_theta)

        // vce(robust): cluster-robust (Lin & Wei 1989 style) sandwich for
        // beta from the per-record score residual
        //   U_r = wt_r*d_r*(x_r - Sbar1(t_r))                              [event term]
        //       + wt_r*exp(x_r'b+offset_r)*[rowDbar1_r - x_r*rowH0diff_r]  [at-risk term]
        // REV 38: the at-risk term uses the hazard INCLUDING the frailty
        // offset (r0_robust); infl_rec above (no offset) is a different
        // object -- a derivative wrt u -- and stays untouched for Mbb.
        // REV 48: scores summed by vce cluster (clusvar; default = the
        // shared() groups, i.e. vce(robust)) and the sandwich multiplied by
        // Gc/(Gc-1), Stata's small-sample factor for cluster-robust VCEs.
        // The model-based (oim) V is kept as V_modelbased (stcox convention).
        real scalar Gc
        real vector cgidx
        if (vcemode == 2 & clusvar != "") cgidx = stcoxem_group(st_data(., clusvar, touse), Gc)
        else {
            cgidx = gidx
            Gc = G
        };
        if (vcemode == 2) st_numscalar("r(N_clust_vce)", Gc);
        if (vcemode == 2 & p>0) {
            real matrix Ucluster, Bread, Meat, eventterm, atriskterm, Utot
            real vector r0_robust
            st_matrix("stcoxem_Vmb", Vbeta)
            r0_robust  = wt :* exp(X*beta :+ offset)
            eventterm  = (wt :* d) :* (X - rowSbar1_all)
            atriskterm = r0_robust :* (rowDbar1_all - X:*H0diff_all)
            Utot       = eventterm + atriskterm
            Ucluster   = J(Gc,p,0)
            for (r=1; r<=N; r++) {
                Ucluster[cgidx[r],.] = Ucluster[cgidx[r],.] + Utot[r,.]
            }
            Bread = invsym(infotot)
            Meat  = cross(Ucluster, Ucluster)
            Vbeta = (Gc/(Gc-1)) * (Bread * Meat * Bread)
            st_matrix("stcoxem_V", Vbeta)
            st_matrixcolstripe("stcoxem_V", (J(p,1,""), xnames'))
            st_matrixrowstripe("stcoxem_V", (J(p,1,""), xnames'))
        };
    }
    else if (vcemode == 3) {
        // ---- vce(bootstrap) (rev 17): nonparametric CLUSTER bootstrap over
        //      shared() groups. G clusters drawn with replacement, each
        //      DRAWN SLOT getting a fresh cluster id (like Stata's
        //      "bootstrap, cluster() idcluster()"). Each replicate calls
        //      stcoxem_em_core silently, cold-started by default (rev 24;
        //      bswarmstart opts into a warm start). ----
        real scalar rep, slot, nboot_used, theta_r, iter_r, converged_r, ll_r
        real scalar theta_bar
        real vector draw, ridx_all, newgid_all, ridx_c
        real vector bs_t0, bs_t, bs_d, bs_wt, bs_strat, bs_gidx
        real matrix bs_X
        real vector beta_r, H0diff_r, Hi_r, Di_r, offset_r
        real vector theta_reps, beta_bar_row
        real matrix beta_reps, beta_centered
        real scalar nt_att_r, nt_acc_r, nb_att_r, nb_acc_r, nsq_att_r, nsq_acc_r
        real scalar showprog, prog_every
        real scalar n_maxiter, n_badtheta, n_badbeta, sum_iter
        real vector theta_all, iter_all, outcome_all
        real matrix beta_all
        transmorphic clus_rows

        if (bsseed >= 0) rseed(bsseed)

        clus_rows = asarray_create("real", 1)
        for (r=1; r<=G; r++) {
            asarray(clus_rows, r, selectindex(gidx:==r))
        }

        theta_reps = J(reps,1,.)
        beta_reps  = J(reps, (p>0 ? p : 1), .)
        nboot_used = 0
        n_maxiter  = 0
        n_badtheta = 0
        n_badbeta  = 0
        sum_iter   = 0
        // raw per-replicate results (rev 23); outcome codes:
        // 1=usable, 2=hit maxiter, 3=invalid theta, 4=invalid beta
        theta_all   = J(reps,1,.)
        iter_all    = J(reps,1,.)
        outcome_all = J(reps,1,.)
        beta_all    = J(reps, (p>0 ? p : 1), .)

        showprog = (dolog | dodots)
        prog_every = max((1, floor(reps/10)))
        if (showprog) printf("  Bootstrap replications (%g): ", reps)

        for (rep=1; rep<=reps; rep++) {
            draw = ceil(runiform(G,1):*G)

            ridx_all   = J(0,1,.)
            newgid_all = J(0,1,.)
            for (slot=1; slot<=G; slot++) {
                ridx_c     = asarray(clus_rows, draw[slot])
                ridx_all   = ridx_all \ ridx_c
                newgid_all = newgid_all \ J(rows(ridx_c),1,slot)
            }

            bs_t0    = t0[ridx_all]
            bs_t     = t[ridx_all]
            bs_d     = d[ridx_all]
            bs_wt    = wt[ridx_all]
            bs_strat = strat[ridx_all]
            bs_gidx  = newgid_all
            bs_X     = (p>0 ? X[ridx_all,.] : J(rows(ridx_all),0,.))

            stcoxem_em_core(bs_t0, bs_t, bs_d, bs_wt, bs_gidx, G, bs_strat, S,
                bs_X, p, rows(ridx_all), maxiter, coxiter, tol, ltol, accelmode, 0, 0, 0,
                (usebswarm ? beta : J(p,1,0)), (usebswarm ? theta : 0.5),
                beta_r, theta_r, iter_r, converged_r, H0diff_r, Hi_r, Di_r, ll_r,
                offset_r, nt_att_r, nt_acc_r, nb_att_r, nb_acc_r, nsq_att_r, nsq_acc_r)

            sum_iter = sum_iter + iter_r
            theta_all[rep] = theta_r
            iter_all[rep]  = iter_r
            if (p>0 & !hasmissing(beta_r)) beta_all[rep,.] = beta_r'
            if (theta_r>=. | theta_r<=0) {
                n_badtheta++
                outcome_all[rep] = 3
            }
            else if (p>0 & hasmissing(beta_r)) {
                n_badbeta++
                outcome_all[rep] = 4
            }
            else if (converged_r != 1) {
                n_maxiter++
                outcome_all[rep] = 2
            }
            else {
                nboot_used++
                theta_reps[nboot_used] = theta_r
                if (p>0) beta_reps[nboot_used,.] = beta_r'
                outcome_all[rep] = 1
            };

            if (showprog) {
                // displayflush() so dots appear one at a time
                printf(".")
                if (mod(rep,prog_every)==0 & rep!=reps) printf(" %3.0f%%", 100*rep/reps)
                displayflush()
            }
        }
        if (showprog) printf(" 100%%\n")

        real scalar avg_iter
        avg_iter = sum_iter / reps
        printf("  Average EM iterations per bootstrap replicate = %6.1f" +
               " (maxiter = %g)\n", avg_iter, maxiter)

        if (nboot_used < 2) {
            printf("  WARNING: vce(bootstrap) had only %g of %g replicates" +
                   " converge to a usable estimate; results below are" +
                   " unreliable -- consider more reps or checking why" +
                   " replicates are failing to converge.\n", nboot_used, reps)
            se_theta = .
            if (p>0) Vbeta = J(p,p,.)
        }
        else {
            theta_reps = theta_reps[1::nboot_used]
            theta_bar  = mean(theta_reps)
            se_theta   = sqrt( sum((theta_reps :- theta_bar):^2) / (nboot_used-1) )
            if (p>0) {
                beta_reps     = beta_reps[1::nboot_used,.]
                beta_bar_row  = mean(beta_reps)
                beta_centered = beta_reps :- beta_bar_row
                Vbeta = (beta_centered' * beta_centered) / (nboot_used-1)
            }
        }
        if (nboot_used < reps) {
            printf("  Bootstrap outcomes: %g usable, %g hit maxiter without" +
                   " converging, %g had an invalid theta, %g had an invalid" +
                   " beta component.\n", nboot_used, n_maxiter, n_badtheta, n_badbeta)
        }

        st_numscalar("r(se_theta)", se_theta)
        st_numscalar("r(boot_reps)", reps)
        st_numscalar("r(boot_reps_used)", nboot_used)
        st_numscalar("r(boot_reps_maxiter)", n_maxiter)
        st_numscalar("r(boot_reps_badtheta)", n_badtheta)
        st_numscalar("r(boot_reps_badbeta)", n_badbeta)
        st_numscalar("r(boot_avg_iter)", avg_iter)
        st_matrix("boot_theta_all", theta_all)
        st_matrix("boot_iter_all", iter_all)
        st_matrix("boot_outcome_all", outcome_all)
        if (p>0) {
            st_matrix("boot_beta_all", beta_all)
            st_matrixcolstripe("boot_beta_all", (J(p,1,""), xnames'))
        }
        if (p>0) {
            st_matrix("stcoxem_V", Vbeta)
            st_matrixcolstripe("stcoxem_V", (J(p,1,""), xnames'))
            st_matrixrowstripe("stcoxem_V", (J(p,1,""), xnames'))
        }
    };

    // ---- mean cumulative baseline hazard: SUM each subject's record-level
    //      increments, THEN average across subjects (not records) ----
    real scalar H0mean
    real vector H0_subj
    // rev 47a: baseline back from x = xbar to x = 0
    if (p>0) H0diff_all = H0diff_all :* exp(-(xbar*beta))
    H0_subj = J(GS,1,0)
    for (r=1; r<=N; r++) {
        H0_subj[sgidx[r]] = H0_subj[sgidx[r]] + H0diff_all[r]
    }
    H0mean = mean(H0_subj)
    st_store(., genvar, touse, H0diff_all)
    if (ufvar != "") st_store(., ufvar, touse, exp(offset))

    st_numscalar("r(theta)", theta)
    st_numscalar("r(ll)", ll)
    st_numscalar("r(N_orig)", N)
    st_numscalar("r(N_subj)", GS)
    st_numscalar("r(N_clust)", G)
    st_numscalar("r(H0mean)", H0mean)
    st_numscalar("r(iter)", iter<=maxiter ? iter : maxiter)
    st_numscalar("r(converged)", converged)

    // ---- memory footprint estimate (lower bound, not peak RSS) ----
    real scalar mem_bytes
    mem_bytes = 8 * (rows(X)*cols(X) + 6*N)
    if ((vcemode==1 | vcemode==2) & p>0) mem_bytes = mem_bytes + 8*3*N*p;
    if (vcemode==3 & p>0) mem_bytes = mem_bytes + 8*reps*p;
    st_numscalar("r(mem_est_mb)", mem_bytes / (1024*1024))
}

end
