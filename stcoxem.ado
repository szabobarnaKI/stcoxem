*! stcoxem.ado — Semiparametric Cox model with shared gamma frailty,
*! fit by the Klein (1992) EM algorithm. Stata 19.0 / Mata.
*!
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
*!        [ strata(varname) vce(none|oim|robust|bootstrap) coxiter(#)
*!          generate(newvar) maxiter(#) tolerance(#) nolog dots
*!          accelerate(squarem|aitken|none) noaitken level(#) nohr
*!          reps(#) bsseed(#) bswarmstart plci noplse ufrailty(newvar) ]
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

program define stcoxem, eclass
    version 19.0

    syntax [varlist(min=0 numeric fv ts)] [if] [in] , ///
        SHared(varname) ///
        [ STRata(varname) VCE(string) COXiter(integer 1) GENerate(name) ///
          MAXiter(integer 200) TOLerance(real 1e-6) noLOg DOTS ///
          ACCELerate(string) noAItken LEVEL(cilevel) noHR ///
          REPS(integer 50) BSSEED(integer -1) BSWarmstart PLci noPLse ///
          UFrailty(name) ]

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

    if "`vce'" == "" local vce "none"
    if !inlist("`vce'", "none", "oim", "robust", "bootstrap") {
        di as err "vce() must be none, oim, robust, or bootstrap"
        exit 198
    }
    local vcemode = cond("`vce'"=="none", 0, cond("`vce'"=="oim", 1, cond("`vce'"=="robust", 2, 3)))
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
        di as txt "(weights: none)"
    }

    * ---- sample ----------------------------------------------------------
    marksample touse
    quietly replace `touse' = 0 if _st != 1
    * rev 43: strok -- markout otherwise treats ANY string variable in its
    * varlist as unconditionally missing (undocumented in -help markout-'s
    * one-line synopsis, but see the strok entry), which silently zeroed out
    * `touse' whenever strata() (or, in principle, shared()) was a string.
    markout `touse' `shared' `strata' _t0 _t _d `wt', strok
    quietly count if `touse'
    if r(N) == 0 {
        di as err "no observations"
        exit 2000
    }

    * ---- factor-variable / time-series expansion ------------------------
    * xvars_full (fvexpand) and xvars_full_revar (fvrevar, DEFAULT mode) are
    * built from the same expression, one token per term in the same order;
    * base-level terms are dropped from BOTH lists at the same position.
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
        local xvars_expand
        local xvars_revar
        forvalues i = 1/`nfull' {
            local v  : word `i' of `xvars_full'
            local rv : word `i' of `xvars_full_revar'
            if !regexm("`v'", "^[0-9]+b\.") {
                local xvars_expand `xvars_expand' `v'
                local xvars_revar  `xvars_revar'  `rv'
            }
        }
        local p         : word count `xvars_revar'
        local p_expand  : word count `xvars_expand'
        if `p_expand' != `p' {
            di as err "internal error: xvars_revar/xvars_expand length" ///
                " mismatch after base-level filtering -- please report" ///
                " this along with your varlist."
            exit 198
        }
    }
    else {
        local xvars_revar
        local xvars_expand
        local p 0
    }

    * ---- numeric group ids for shared() and strata() -------------------
    capture confirm numeric variable `shared'
    if _rc {
        tempvar sharedid
        quietly egen long `sharedid' = group(`shared') if `touse'
    }
    else local sharedid `shared'

    if "`strata'" != "" {
        * rev 42: ALWAYS remap to consecutive codes 1..S. The Mata code loops
        * s=1..S and selects strat==s, so a numeric strata variable coded
        * 0/1, negative, fractional or with gaps would have silently lost
        * (some of) its strata from every risk set.
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
    local dohr       = ("`hr'" == "")

    * ---- plci(): profile-likelihood CI for theta ----
    local doplci = ("`plci'" != "")
    if `doplci' local chi2crit = invchi2(1, `level'/100)
    else        local chi2crit = .

    * ---- plse: default ON; a declared noPLse option populates local `plse'
    *      with the string "noplse" when given (rev 37) ----
    local doplse = cond("`plse'" == "noplse", 0, 1)

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

    * ================= hand off to Mata =================================
    mata: stcoxem_driver( ///
        "`xvars_revar'", "`sharedid'", "`stratid'", "`wt'", "`touse'", "`subjid'", ///
        `S', `maxiter', `coxiter', `tolerance', `dolog', `dodots', `accelmode', ///
        `vcemode', "`generate'", `reps', `bsseed', `usebswarm', `doplci', `chi2crit', ///
        `doplse', "`ufrailty'")

    if `dogen' {
        label variable `generate' "record-level Breslow baseline cumulative hazard, stcoxem"
    }
    if `doufrailty' {
        label variable `ufrailty' "record-level posterior mean frailty E[u|data] (shared() cluster), stcoxem"
    }

    * ================= ereturn + display =================================
    if `p' > 0 & `vcemode' >= 1 {
        ereturn post b V, esample(`touse') obs(`=r(N_orig)')
    }
    else if `p' > 0 {
        ereturn post b, esample(`touse') obs(`=r(N_orig)')
    }
    else {
        ereturn post, esample(`touse') obs(`=r(N_orig)')
    }
    if `p' > 0 {
        tempname btmp
        matrix `btmp' = e(b)
        matrix colnames `btmp' = `xvars_expand'
        if `vcemode' >= 1 {
            tempname Vtmp
            matrix `Vtmp' = e(V)
            matrix colnames `Vtmp' = `xvars_expand'
            matrix rownames `Vtmp' = `xvars_expand'
            ereturn repost b = `btmp' V = `Vtmp', rename
        }
        else {
            ereturn repost b = `btmp', rename
        }
        ereturn local depvar "_t"
    }
    ereturn local  cmd       "stcoxem"
    ereturn local  title     "Semiparametric Cox shared gamma-frailty model (EM)"
    ereturn local  sharedvar "`shared'"
    ereturn local  vce       "`vce'"
    ereturn local  accelerate "`accelerate'"
    if "`strata'" != "" ereturn local stratavar "`strata'"
    if "`wtype'" == "pweight" ereturn local wexp "`wvar'"
    ereturn scalar N_clust   = r(N_clust)
    ereturn scalar N_orig    = r(N_orig)
    ereturn scalar N_subj    = r(N_subj)
    ereturn scalar n_strata  = `S'
    ereturn scalar theta     = r(theta)
    * e(se_theta): plse's profile SE when it ran and converged reliably,
    * else the analytic Louis/EM value (always kept in e(se_theta_oim)).
    if `vcemode' >= 1 {
        ereturn scalar se_theta_oim = r(se_theta)
        local use_plse_se = 0
        if `doplse' {
            local se_theta_pl = r(se_theta_profile)
            if `se_theta_pl' < . local use_plse_se = 1
        }
        if `use_plse_se' {
            ereturn scalar se_theta = `se_theta_pl'
            ereturn local  se_theta_source "plse"
        }
        else {
            ereturn scalar se_theta = r(se_theta)
            ereturn local  se_theta_source = cond(`doplse', ///
                "oim (plse did not converge reliably for this fit)", "oim (noplse)")
        }
    }
    ereturn scalar ll        = r(ll)
    ereturn scalar H0mean    = r(H0mean)
    ereturn scalar iter      = r(iter)
    ereturn scalar converged = r(converged)
    ereturn scalar mem_est   = r(mem_est_mb)
    if `accelmode' == 1 {
        ereturn scalar aitken_theta_att = r(aitken_theta_att)
        ereturn scalar aitken_theta_acc = r(aitken_theta_acc)
        if `p' > 0 {
            ereturn scalar aitken_beta_att = r(aitken_beta_att)
            ereturn scalar aitken_beta_acc = r(aitken_beta_acc)
        }
    }
    else if `accelmode' == 2 {
        ereturn scalar squarem_att = r(squarem_att)
        ereturn scalar squarem_acc = r(squarem_acc)
    }
    if (`vcemode'==1 | `vcemode'==2) & `p' > 0 {
        ereturn scalar se_theta_naive = r(se_theta_naive)
        ereturn scalar schur_term     = r(schur_term)
    }
    if `doplci' {
        ereturn scalar theta_pl_lower  = r(theta_pl_lower)
        ereturn scalar theta_pl_upper  = r(theta_pl_upper)
        ereturn scalar theta_pl_nonconv = r(theta_pl_nonconv)
    }
    if `doplse' {
        ereturn scalar se_theta_profile = r(se_theta_profile)
        ereturn scalar theta_profile_nonconv = r(theta_profile_nonconv)
    }
    if `vcemode' == 3 {
        ereturn scalar boot_reps          = r(boot_reps)
        ereturn scalar boot_reps_used     = r(boot_reps_used)
        ereturn scalar boot_reps_maxiter  = r(boot_reps_maxiter)
        ereturn scalar boot_reps_badtheta = r(boot_reps_badtheta)
        ereturn scalar boot_reps_badbeta  = r(boot_reps_badbeta)
        ereturn scalar boot_avg_iter      = r(boot_avg_iter)
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
            ereturn matrix boot_beta_all = `boot_beta_all'
        }
    }

    di ""
    di as txt "Semiparametric Cox frailty model (gamma frailty), Klein EM algorithm"
    di as txt "{hline 70}"
    di as txt "Number of records (obs) = " as res %9.0f e(N_orig)
    di as txt "Number of frailty groups = " as res %9.0f e(N_clust) ///
        as txt "  (" as res "`shared'" as txt ")"
    if `nosubjid' {
        di as txt "(no id() found in -stset-; treating each record as its own" ///
            " subject for H0mean -- set your data up with" ///
            " -stset ..., id(subjectvar)- if subjects have multiple records)"
    }
    else di as txt "Number of subjects      = " as res %9.0f e(N_subj)
    if "`strata'" != "" di as txt "Number of strata = " as res %9.0f e(n_strata) ///
        as txt "  (" as res "`strata'" as txt ")"
    if "`wtype'" == "pweight" di as txt "Weights: pweight = " as res "`wvar'" as txt " (from -stset-)"
    di as txt "Cox partial log-likelihood + frailty terms = " as res %9.3f e(ll)
    di as txt "EM iterations            = " as res %9.0f e(iter) ///
       as txt "   converged = " as res cond(e(converged)==1, "yes", "no") ///
       as txt "   accelerate = " as res "`accelerate'"
    * Est. memory footprint (lower bound) is no longer printed here -- still
    * available via e(mem_est).

    * ================= coefficient + theta table ==========================
    local critval = invnormal(1 - (100-`level')/200)
    local hrlabel = cond(`dohr', "Haz. ratio", "Coef.")
    local vcelabel = cond(`vcemode'==0, "none (no SEs computed)", ///
                      cond(`vcemode'==1, "oim (joint beta/theta Louis/EM correction)", ///
                      cond(`vcemode'==2, "robust (cluster-robust beta; joint Louis theta)", ///
                                         "bootstrap (nonparametric cluster bootstrap over shared() groups, `reps' reps)")))

    tempname bmat
    matrix `bmat' = e(b)
    local ncoef = colsof(`bmat')
    local cnames : colnames `bmat'
    if `vcemode' >= 1 & `p' > 0 {
        tempname Vmat
        matrix `Vmat' = e(V)
    }

    di as txt "vce: " as res "`vcelabel'"
    if `vcemode' == 3 di as txt "Bootstrap replications used = " as res %4.0f e(boot_reps_used) as txt " of " as res %4.0f e(boot_reps)
    if `vcemode' == 3 & !`usebswarm' di as txt "(replicates cold-started at beta=0, theta=0.5; use bswarmstart to warm-start from the main fit instead, faster but can bias se_theta upward)"
    if `vcemode' == 3 & `usebswarm' di as txt "(bswarmstart: replicates warm-started from the main fit; can converge a resample with a true near-zero theta to a spurious high-theta fixed point; compare against the cold-started default if results look implausible)"
    di as txt "{hline 16}{c TT}{hline 60}"
    di as txt %15s abbrev("_t",15) " {c |} " %10s "`hrlabel'" "  " %9s "Std. err." ///
        "  " %6s "z" "  " %6s "P>|z|" "  " %20s "[`level'% Conf. Interval]"
    di as txt "{hline 16}{c +}{hline 60}"

    forvalues i = 1/`ncoef' {
        local nm : word `i' of `cnames'
        local bval = `bmat'[1,`i']
        local shown = cond(`dohr', exp(`bval'), `bval')
        if `vcemode' >= 1 & `p' > 0 {
            local se = sqrt(`Vmat'[`i',`i'])
            local z  = `bval'/`se'
            local pval = 2*normal(-abs(`z'))
            local lo = `bval' - `critval'*`se'
            local hi = `bval' + `critval'*`se'
            if `dohr' {
                local se_shown = exp(`bval')*`se'
                local lo_shown = exp(`lo')
                local hi_shown = exp(`hi')
            }
            else {
                local se_shown = `se'
                local lo_shown = `lo'
                local hi_shown = `hi'
            }
            di as txt %15s abbrev("`nm'",15) " {c |} " as res %10.4g `shown' "  " %9.4g `se_shown' ///
                "  " %6.2f `z' "  " %6.3f `pval' "  " %9.4g `lo_shown' "  " %9.4g `hi_shown'
        }
        else {
            di as txt %15s abbrev("`nm'",15) " {c |} " as res %10.4g `shown' ///
                as txt "  " %9s "." "  " %6s "." "  " %6s "." "  " %9s "." "  " %9s "."
        }
    }

    di as txt "{hline 16}{c +}{hline 60}"
    if `vcemode' >= 1 {
        local se_th = e(se_theta)
        local zth   = e(theta)/`se_th'
        local pth   = 2*normal(-abs(`zth'))
        local loth  = e(theta) - `critval'*`se_th'
        local hith  = e(theta) + `critval'*`se_th'
        di as txt %15s "theta" " {c |} " as res %10.4g e(theta) "  " %9.4g `se_th' ///
            "  " %6.2f `zth' "  " %6.3f `pth' "  " %9.4g `loth' "  " %9.4g `hith'
    }
    else {
        di as txt %15s "theta" " {c |} " as res %10.4g e(theta) ///
            as txt "  " %9s "." "  " %6s "." "  " %6s "." "  " %9s "." "  " %9s "."
    }
    di as txt "{hline 16}{c BT}{hline 60}"
    if `dohr' di as txt "(theta shown on its natural scale, not exponentiated)"

end


* =====================================================================
*                              MATA CODE
* =====================================================================
mata:
mata set matastrict off

// ---------------------------------------------------------------------
// Map an arbitrary numeric id vector to consecutive group indices 1..G.
// ---------------------------------------------------------------------
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
    if (a > 1e4) {
        return( d:*ln(a) :+ (d:*(d:-1)):/(2*a) :- (d:*(d:-1):*(2*d:-1)):/(12*a^2) )
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
        if (abs(denomt) > 1e-10) {
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
                if (abs(denomb[kk2]) > 1e-10) {
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
    real scalar dosens

    n = rows(t)
    p = cols(X)
    r = wt :* exp(eta)
    dosens = (needsens != 0 & p > 0)

    order_entry = order(t0, 1)
    order_exit  = order(t, 1)

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
        }

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
        real scalar accelmode, real scalar fixtheta,
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
    real scalar thetaold, llold
    transmorphic strat_rows

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
            if (p>0) {
                betaold = beta
                beta = beta + invsym(infotot) * scoretot
            };
        }

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
        //      ll is then correct up to an additive tie-pattern constant:
        //      fine for within-model differences (trace, plci, plse), NOT
        //      for cross-model LR/AIC comparisons. ----
        real scalar a
        real vector lgratio
        a = 1/theta
        lgratio = stcoxem_lgammaratio(Di, a)
        // rev 40: cluster-level terms scaled by Wi (cluster replication
        // weight); Hi/Di are per-copy (weight-normalized) quantities.
        ll = pl_total - sum(Wi :* Di :* ln(w1)) + sum( Wi :* ( lgratio :- Di:*ln(a)
                              :- (Di :+ a) :* stcoxem_log1p(Hi:/a) ) )

        if (dolog) {
            printf("  EM iter %3.0f:  logL = %12.4f   theta = %8.5f\n", iter, ll, theta)
        }
        else if (dodots) {
            printf(".")
        };

        // ---- convergence: parameter change ----
        real scalar parchange
        parchange = abs(theta-thetaold)/(abs(thetaold)+1e-8)
        if (p>0) parchange = max((parchange, max(abs(beta-betaold):/(abs(betaold):+1e-8))));
        if (iter>1 & parchange < tol) {
            converged = 1
            break
        };
        llold = ll
    }
    if (dodots) printf("\n");
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
        real scalar accelmode, real scalar theta_fixed, real vector beta_warm,
        beta_out, conv_out)
{
    real scalar theta_r, iter_r, converged_r, ll_r
    real vector H0diff_r, Hi_r, Di_r, offset_r
    real scalar nt_att_r, nt_acc_r, nb_att_r, nb_acc_r, nsq_att_r, nsq_acc_r

    stcoxem_em_core(t0, t, d, wt, gidx, G, strat, S, X, p, N,
        maxiter, coxiter, tol, accelmode, 1, 0, 0,
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
        real scalar accelmode, real scalar theta_hat, real vector beta_hat,
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
            maxiter,coxiter,tol,accelmode,theta_out,beta_cur, beta_new,conv)
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
            maxiter,coxiter,tol,accelmode,mid,beta_cur, beta_new,conv)
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
void stcoxem_driver(
        string scalar xvars, string scalar sharedvar, string scalar stratvar,
        string scalar wtvar, string scalar touse, string scalar subjvar,
        real scalar S, real scalar maxiter, real scalar coxiter,
        real scalar tol, real scalar dolog, real scalar dodots,
        real scalar accelmode, real scalar vcemode,
        string scalar genvar, real scalar reps, real scalar bsseed,
        real scalar usebswarm, real scalar doplci, real scalar chi2crit,
        real scalar doplse, string scalar ufvar)
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

    gidx = stcoxem_group(cid, G)
    real scalar GS
    real vector sgidx
    sgidx = stcoxem_group(subj, GS)

    // ---- main fit: cold-started at beta=0, theta=0.5 ----
    stcoxem_em_core(t0, t, d, wt, gidx, G, strat, S, X, p, N,
        maxiter, coxiter, tol, accelmode, 0, dolog, dodots, J(p,1,0), 0.5,
        beta, theta, iter, converged, H0diff_all, Hi, Di, ll, offset,
        n_theta_att, n_theta_acc, n_beta_att, n_beta_acc, n_sq_att, n_sq_acc)

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
        target_ll = ll - chi2crit/2
        n_pl_nonconv = 0

        if (showprog_pl) printf("  Profile-likelihood CI for theta: ")
        theta_pl_lower = stcoxem_profile_bound(t0, t, d, wt, gidx, G, strat, S,
            X, p, N, maxiter, coxiter, tol, accelmode, theta, beta,
            target_ll, -1, 1e-6, showprog_pl, n_pl_nonconv)
        theta_pl_upper = stcoxem_profile_bound(t0, t, d, wt, gidx, G, strat, S,
            X, p, N, maxiter, coxiter, tol, accelmode, theta, beta,
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
        n_prof_nonconv = 0

        real scalar ll_c, conv_c
        real vector beta_c
        ll_c = stcoxem_profile_ll(t0,t,d,wt,gidx,G,strat,S,X,p,N,
            maxiter,coxiter,tol_pl,accelmode, theta, beta, beta_c, conv_c)
        if (!conv_c) n_prof_nonconv = n_prof_nonconv + 1

        ll_p1 = stcoxem_profile_ll(t0,t,d,wt,gidx,G,strat,S,X,p,N,
            maxiter,coxiter,tol_pl,accelmode, theta+h1, beta_c, beta_p1, conv_p1)
        if (!conv_p1) n_prof_nonconv = n_prof_nonconv + 1
        ll_m1 = stcoxem_profile_ll(t0,t,d,wt,gidx,G,strat,S,X,p,N,
            maxiter,coxiter,tol_pl,accelmode, theta-h1, beta_c, beta_m1, conv_m1)
        if (!conv_m1) n_prof_nonconv = n_prof_nonconv + 1
        ll_p2 = stcoxem_profile_ll(t0,t,d,wt,gidx,G,strat,S,X,p,N,
            maxiter,coxiter,tol_pl,accelmode, theta+h2, beta_c, beta_p2, conv_p2)
        if (!conv_p2) n_prof_nonconv = n_prof_nonconv + 1
        ll_m2 = stcoxem_profile_ll(t0,t,d,wt,gidx,G,strat,S,X,p,N,
            maxiter,coxiter,tol_pl,accelmode, theta-h2, beta_c, beta_m2, conv_m2)
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
        st_matrix("b", beta')
        st_matrixcolstripe("b", (J(p,1,""), xnames'))
    }
    else {
        st_matrix("b", J(1,0,.))
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

            st_matrix("V", Vbeta)
            st_matrixcolstripe("V", (J(p,1,""), xnames'))
            st_matrixrowstripe("V", (J(p,1,""), xnames'))
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
        if (vcemode == 2 & p>0) {
            real matrix Ucluster, Bread, Meat, eventterm, atriskterm, Utot
            real vector r0_robust
            r0_robust  = wt :* exp(X*beta :+ offset)
            eventterm  = (wt :* d) :* (X - rowSbar1_all)
            atriskterm = r0_robust :* (rowDbar1_all - X:*H0diff_all)
            Utot       = eventterm + atriskterm
            Ucluster   = J(G,p,0)
            for (r=1; r<=N; r++) {
                Ucluster[gidx[r],.] = Ucluster[gidx[r],.] + Utot[r,.]
            }
            Bread = invsym(infotot)
            Meat  = cross(Ucluster, Ucluster)
            Vbeta = Bread * Meat * Bread
            st_matrix("V", Vbeta)
            st_matrixcolstripe("V", (J(p,1,""), xnames'))
            st_matrixrowstripe("V", (J(p,1,""), xnames'))
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
                bs_X, p, rows(ridx_all), maxiter, coxiter, tol, accelmode, 0, 0, 0,
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
            st_matrix("V", Vbeta)
            st_matrixcolstripe("V", (J(p,1,""), xnames'))
            st_matrixrowstripe("V", (J(p,1,""), xnames'))
        }
    };

    // ---- mean cumulative baseline hazard: SUM each subject's record-level
    //      increments, THEN average across subjects (not records) ----
    real scalar H0mean
    real vector H0_subj
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
