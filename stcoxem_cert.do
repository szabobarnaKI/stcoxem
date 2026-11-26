* stcoxem_cert.do -- certification script for stcoxem 1.0.0 (rev 55, 2026-11-26)
* rev 55: Aitken probe removed from c0; c11 asserts that accelerate(aitken),
*   squarem and none converge to the same theta on simulated and registry data.
* rev 53: simulated data without delayed entry by default (stcox refuses
*   shared() with delayed entry: 498); true frailty renamed utrue; c0 convergence
*   probe (informational); c6 and c10 use delayed entry with forceshared.
*
* Runs the package as a user would get it: the files are found on the ado
* path (autoload), not -run-. Simulated data (fixed seed), so the script is
* self-contained; an optional block reruns the key comparison on the user's
* registry test set if $REALDATA points to it.
*
* Before running, put in $TDIR (from claude_communication ado/stcoxem/):
*   stcoxem.ado  stcoxem_p.ado  stcoxem.sthlp  stcoxem_postestimation.sthlp
*
*   c1 e() layout, version, help files present
*   c2 estimates vs stcox, shared(): b, theta, ll, ll_c, chi2_c, se_theta, V
*   c3 predict vs stcox, shared() and vs ufrailty()/generate()
*   c4 vce(robust)/vce(cluster) == stcox, offset(ln u) vce(cluster)
*   c5 options: strata, none, bootstrap, plci, noplse, nolrtest, coxiter,
*      accelerate, replay/display options, no covariates, collinearity
*   c6 integer pweights == expanded data (k copies, new group ids)
*   c7 postestimation: estat, lrtest vs stcox's, margins, test, stcurve (PNG)
*   c8 error codes
*   c9 (optional) registry test set vs stcox, shared()
*   c10 delayed entry: 498 without forceshared; with it vs stcox forceshared
*   c0 (informational) EM iterations by tolerance and acceleration
*   c11 squarem, aitken, none all converge to the tight solution

version 19
clear all
discard
set more off

global TDIR     "C:\Barna\PhD/ado development/stcoxem/most recent"
global REALDATA "$TDIR/test_set.dta"
global TIGHT    "tolerance(1e-11) maxiter(5000)"
global FAILED   ""
global PNGS     ""
local logf "$TDIR/stcoxem_cert.log"

capture log close tlog
log using "`logf'", text replace name(tlog)
di as txt "Stata " c(stata_version) "  " c(os) "  " c(current_date) " " c(current_time)

* ---- package from $TDIR, first on the ado path ----
adopath ++ "$TDIR"
capture noisily which stcoxem
capture noisily which stcoxem_p

* ======================= helpers ========================================
* simulated data: 250 groups of 1-5 subjects, one record each, gamma
* frailty (theta 0.6), Weibull baseline, censoring, times rounded to 0.05
* (ties); "simdata late": 30% delayed entry
capture program drop simdata
program define simdata
    args late
    clear
    set seed 20261126
    quietly set obs 250
    gen long gid = _n
    gen double utrue = rgamma(1/0.6, 0.6)
    gen byte n = 1 + floor(runiform()*5)
    quietly expand n
    sort gid
    gen long id = _n
    gen double x1 = rnormal()
    gen byte x2 = runiform() < 0.4
    gen byte x3 = 1 + floor(3*runiform())
    gen double eta = 0.5*x1 - 0.4*x2 + 0.3*(x3==2) + 0.6*(x3==3)
    gen double t = (-ln(runiform())/(0.1*utrue*exp(eta)))^(1/1.2)
    gen double c = 12*runiform()
    gen double y = ceil(min(t, c)*20)/20
    gen byte d = t <= c
    gen double t0 = cond(runiform() < 0.3, floor(0.5*y*20)/20, 0)
    if "`late'" == "" quietly replace t0 = 0
    gen long cl = ceil(gid/10)
    gen byte w = 1 + mod(gid, 3)
    drop t c eta n
    quietly stset y, failure(d) id(id) time0(t0)
end

* check "label" condition
capture program drop check
program define check
    gettoken lbl 0 : 0
    di as txt "  `lbl'"
    assert `0'
end

* cmps a b tol "label": reldif of two scalars/expressions
capture program drop cmps
program define cmps
    args a b tol lbl
    tempname x y
    scalar `x' = `a'
    scalar `y' = `b'
    di as txt "  " %-40s "`lbl'" as res %14.0g `x' %14.0g `y' ///
        as txt "  reldif " %9.2e reldif(`x', `y')
    assert reldif(`x', `y') < `tol'
end

* cmpv a b tol "label": variables, same missing pattern, max reldif < tol
capture program drop cmpv
program define cmpv
    args a b tol lbl
    tempvar rd
    quietly count if missing(`a') != missing(`b')
    local nm = r(N)
    quietly count if !missing(`a')
    local nn = r(N)
    quietly gen double `rd' = reldif(`a', `b') if !missing(`a') & !missing(`b')
    quietly summarize `rd', meanonly
    local mx = cond(r(N) == 0, 0, r(max))
    di as txt "  " %-40s "`lbl'" " n=" %5.0f `nn' "  max reldif=" %9.2e `mx' ///
        "  missing-pattern diffs=" `nm'
    assert `nm' == 0 & `nn' > 0 & `mx' < `tol'
end

* cmpm A B tol "label": mreldif of two matrices (same stripe asserted)
capture program drop cmpm
program define cmpm
    args A B tol lbl
    local ca : colfullnames `A'
    local cb : colfullnames `B'
    di as txt "  " %-40s "`lbl'" "  mreldif " %9.2e mreldif(`A', `B')
    assert `"`ca'"' == `"`cb'"'
    assert mreldif(`A', `B') < `tol'
end

* rcis "label" #rc command
capture program drop rcis
program define rcis
    gettoken lbl 0 : 0
    gettoken want 0 : 0
    capture noisily `0'
    local got = _rc
    di as txt "  `lbl': rc = `got' (expected `want')"
    assert `got' == `want'
end

capture program drop runtest
program define runtest
    di as res _n "===== `0' ====="
    capture noisily `0'
    if _rc {
        di as err "`0' FAILED (rc = " _rc ")"
        global FAILED $FAILED `0'
    }
    else di as res "`0' passed"
end

* ======================= tests ==========================================
* C1: version, e() layout, help files
capture program drop c1
program define c1
    findfile stcoxem.ado
    local f `"`r(fn)'"'
    mata: l_ = cat(st_local("f")); st_local("v1", l_[1])
    di as txt "  first line: `v1'"
    check "a: version 1.0.0 loaded" strpos(`"`v1'"', "version 1.0.0") > 0
    findfile stcoxem_p.ado
    findfile stcoxem.sthlp
    findfile stcoxem_postestimation.sthlp
    simdata
    stcoxem x1 x2 i.x3, shared(gid) nolog
    check "b: e(cmd) stcox_fr"      "`e(cmd)'" == "stcox_fr"
    check "c: e(cmd2) stcox"        "`e(cmd2)'" == "stcox"
    check "d: e(estimator) em"      "`e(estimator)'" == "em"
    check "e: e(predict) stcox_p"   "`e(predict)'" == "stcox_p"
    check "f: e(estat_cmd)"         "`e(estat_cmd)'" == "stcox_estat"
    check "g: e(shared) gid"        "`e(shared)'" == "gid"
    check "h: e(vce) oim (default)" "`e(vce)'" == "oim"
    check "i: e(se_theta_source)"   "`e(se_theta_source)'" == "plse"
    check "j: e(N_g) 250"           e(N_g) == 250
    check "k: converged"            e(converged) == 1
    check "l: cmdline"              substr("`e(cmdline)'", 1, 8) == "stcoxem "
    check "m: datasignature"        "`e(datasignature)'" != ""
    local cn : colfullnames e(b)
    di as txt "  stripe: `cn'"
    check "n: stripe has base level" strpos("`cn'", "1b.x3") > 0
    check "o: e(H0mean) nonmissing"  !missing(e(H0mean))
end

* C2: estimates vs stcox, shared() frailty(gamma)
capture program drop c2
program define c2
    simdata
    quietly stcox x1 x2 i.x3, shared(gid) breslow
    tempname bs Vs
    matrix `bs' = e(b)
    matrix `Vs' = e(V)
    local th = e(theta)
    local seth = e(se_theta)
    local ll = e(ll)
    local llc = e(ll_c)
    local chic = e(chi2_c)
    local Ns = e(N_sub)
    local Nf = e(N_fail)
    local rk = e(risk)
    stcoxem x1 x2 i.x3, shared(gid) nolog
    tempname be Ve
    matrix `be' = e(b)
    matrix `Ve' = e(V)
    cmpm `be' `bs' 1e-5 "a: b vs stcox"
    cmps e(theta) `th' 1e-4 "b: theta vs stcox"
    cmps e(ll) `ll' 1e-6 "c: ll vs stcox"
    cmps e(ll_c) `llc' 1e-8 "d: ll_c vs stcox"
    cmps e(chi2_c) `chic' 1e-3 "e: chi2_c vs stcox"
    cmps e(se_theta) `seth' 2e-3 "f: se_theta (plse) vs stcox"
    di as txt "  g: (info) se_theta_oim " e(se_theta_oim) "  plse " e(se_theta)
    cmps e(N_sub) `Ns' 1e-12 "h: N_sub"
    cmps e(N_fail) `Nf' 1e-12 "i: N_fail"
    cmps e(risk) `rk' 1e-12 "j: risk"
    * oim SEs of b within 10% of stcox's
    local k = colsof(`be')
    forvalues j = 1/`k' {
        if `Vs'[`j',`j'] > 0 {
            cmps sqrt(`Ve'[`j',`j']) sqrt(`Vs'[`j',`j']) 0.10 "k: SE[`j'] vs stcox (10%)"
        }
    }
end

* C3: predict -- Path A vs stcox, shared() and vs stcoxem's own variables
capture program drop c3
program define c3
    simdata
    quietly stcox x1 x2 i.x3, shared(gid) breslow
    foreach s in xb hr effects basechazard mgale csnell deviance {
        quietly predict double st_`s', `s'
    }
    quietly predict double st_mgp, mgale partial
    quietly predict double st_sc*, scores
    stcoxem x1 x2 i.x3, shared(gid) nolog $TIGHT ufrailty(u_em) generate(dh_em)
    foreach s in xb hr effects basechazard mgale csnell deviance {
        quietly predict double em_`s', `s'
    }
    quietly predict double em_mgp, mgale partial
    quietly predict double em_sc*, scores
    cmpv em_xb st_xb 1e-5 "a: xb"
    cmpv em_hr st_hr 1e-5 "b: hr"
    cmpv em_effects st_effects 1e-3 "c: effects vs stcox"
    cmpv em_basechazard st_basechazard 1e-4 "d: basechazard vs stcox"
    cmpv em_mgale st_mgale 1e-4 "e: mgale vs stcox"
    cmpv em_mgp st_mgp 1e-4 "f: mgale partial vs stcox"
    cmpv em_csnell st_csnell 1e-4 "g: csnell vs stcox"
    cmpv em_deviance st_deviance 1e-4 "h: deviance vs stcox"
    local k = colsof(e(b))
    forvalues j = 1/`k' {
        capture confirm variable em_sc`j'
        if !_rc {
            quietly count if !missing(em_sc`j')
            if r(N) > 0 cmpv em_sc`j' st_sc`j' 1e-4 "i: scores[`j'] vs stcox"
        }
    }
    gen double lnu_em = ln(u_em)
    cmpv em_effects lnu_em 1e-7 "j: effects == ln ufrailty()"
    * basechazard at _t == running sum of generate() over the distinct times
    * (single-record subjects with delayed entry: H0(_t) = H0(_t0) + dh)
    gen double h0t0 = em_basechazard - dh_em
    quietly count if _t0 == 0 & abs(h0t0) > 1e-8*(1+em_basechazard)
    check "k: basechazard == generate() when _t0==0" r(N) == 0
    check "l: e() intact after predict" "`e(cmd)'" == "stcox_fr" & "`e(estimator)'" == "em"
end

* C4: vce(robust), vce(cluster) == stcox with offset ln u and vce(cluster)
capture program drop c4
program define c4
    simdata
    stcoxem x1 x2 i.x3, shared(gid) nolog $TIGHT vce(robust) ufrailty(u) nolrtest
    tempname bR VR
    matrix `bR' = e(b)
    matrix `VR' = e(V)
    check "a: e(vcetype) Robust"   "`e(vcetype)'" == "Robust"
    check "b: e(clustvar) gid"     "`e(clustvar)'" == "gid"
    check "c: e(N_clust) 250"      e(N_clust) == 250
    capture confirm matrix e(V_modelbased)
    check "d: e(V_modelbased) posted" _rc == 0
    gen double lnu = ln(u)
    quietly stcox x1 x2 i.x3, offset(lnu) breslow vce(cluster gid)
    tempname bS VS
    matrix `bS' = e(b)
    matrix `VS' = e(V)
    cmpm `bR' `bS' 1e-7 "e: b vs offset stcox"
    cmpm `VR' `VS' 1e-6 "f: robust V vs stcox vce(cluster gid)"
    drop u lnu
    stcoxem x1 x2 i.x3, shared(gid) nolog $TIGHT vce(cluster cl) ufrailty(u) nolrtest
    matrix `VR' = e(V)
    check "g: e(N_clust) 25" e(N_clust) == 25
    gen double lnu = ln(u)
    quietly stcox x1 x2 i.x3, offset(lnu) breslow vce(cluster cl)
    matrix `VS' = e(V)
    cmpm `VR' `VS' 1e-6 "h: cluster V vs stcox vce(cluster cl)"
end

* C5: options and display
capture program drop c5
program define c5
    simdata
    rcis "a: strata(x2)" 0 stcoxem x1 i.x3, shared(gid) strata(x2) nolog
    check "   N_strata 2" e(N_strata) == 2
    rcis "b: vce(none)" 0 stcoxem x1 x2 i.x3, shared(gid) vce(none) nolog
    capture confirm matrix e(V)
    check "   no e(V)" _rc != 0
    check "   no se_theta" missing(e(se_theta))
    local thn = e(theta)
    rcis "c: vce(bootstrap)" 0 stcoxem x1 x2 i.x3, shared(gid) vce(bootstrap) reps(5) bsseed(1) nolog
    check "   vcetype Bootstrap" "`e(vcetype)'" == "Bootstrap"
    check "   N_reps <= 5" e(N_reps) <= 5 & e(N_reps) >= 2
    check "   se_theta source bootstrap" "`e(se_theta_source)'" == "bootstrap"
    cmps e(theta) `thn' 1e-12 "   theta same as vce(none)"
    rcis "d: plci" 0 stcoxem x1 x2 i.x3, shared(gid) plci nolog
    check "   lower < theta < upper" e(theta_pl_lower) < e(theta) & e(theta) < e(theta_pl_upper)
    rcis "e: noplse" 0 stcoxem x1 x2 i.x3, shared(gid) noplse nolog
    check "   se_theta == se_theta_oim" e(se_theta) == e(se_theta_oim)
    check "   source oim" substr("`e(se_theta_source)'", 1, 3) == "oim"
    rcis "f: nolrtest" 0 stcoxem x1 x2 i.x3, shared(gid) nolrtest nolog
    check "   no ll_c" missing(e(ll_c))
    stcoxem x1 x2 i.x3, shared(gid) $TIGHT nolog
    tempname b1
    matrix `b1' = e(b)
    local t1 = e(theta)
    foreach o in "coxiter(3)" "accelerate(aitken)" "accelerate(none)" {
        quietly stcoxem x1 x2 i.x3, shared(gid) $TIGHT nolog `o'
        tempname b2
        matrix `b2' = e(b)
        cmpm `b2' `b1' 1e-7 "g: `o': b unchanged"
        cmps e(theta) `t1' 1e-7 "   `o': theta unchanged"
    }
    rcis "h: dots" 0 stcoxem x1 x2 i.x3, shared(gid) dots
    rcis "i: replay" 0 stcoxem
    rcis "j: replay nohr level(90) cformat" 0 stcoxem, nohr level(90) cformat(%9.3f)
    rcis "k: replay noheader nocoef" 0 stcoxem, noheader nocoef
    rcis "l: replay by stcox" 0 stcox
    rcis "m: no covariates" 0 stcoxem, shared(gid) nolog
    capture confirm matrix e(b)
    check "   no e(b)" _rc != 0
    check "   theta > 0" e(theta) > 0
    gen double x4 = 2*x1 + 3
    rcis "n: collinear term omitted" 0 stcoxem x1 x4 x2, shared(gid) nolog
    local cn : colfullnames e(b)
    check "   o.x4 in stripe" strpos("`cn'", "o.x4") > 0
    gen double x1big = x1 + 1e6
    rcis "o: covariate far from 0" 0 stcoxem x1big x2, shared(gid) nolog
    di as txt "   (informational) e(H0mean) = " e(H0mean)
end

* C6: integer pweights == expanded data with new group ids
capture program drop c6
program define c6
    simdata late
    quietly stset y [pweight=w], failure(d) id(id) time0(t0)
    stcoxem x1 x2 i.x3, shared(gid) nolog $TIGHT forceshared
    check "a: e(predict) stcoxem_p" "`e(predict)'" == "stcoxem_p"
    check "b: e(wtype) pweight" "`e(wtype)'" == "pweight"
    tempname bw
    matrix `bw' = e(b)
    local thw = e(theta)
    local llw = e(ll)
    quietly predict double eff_w, effects
    quietly predict double bch_w, basechazard
    quietly predict double mg_w, mgale partial
    * expanded
    gen long rid = _n
    quietly expand w
    bysort rid: gen int copy = _n
    gen long ngid = gid*10 + copy
    gen long nid = id*10 + copy
    quietly stset y, failure(d) id(nid) time0(t0)
    stcoxem x1 x2 i.x3, shared(ngid) nolog $TIGHT forceshared
    tempname bx
    matrix `bx' = e(b)
    cmpm `bw' `bx' 1e-7 "c: b weighted == expanded"
    cmps `thw' e(theta) 1e-7 "d: theta weighted == expanded"
    cmps `llw' e(ll) 1e-9 "e: ll weighted == expanded"
    quietly predict double eff_x, effects
    quietly predict double bch_x, basechazard
    quietly predict double mg_x, mgale partial
    keep if copy == 1
    cmpv eff_w eff_x 1e-6 "f: effects weighted == expanded"
    cmpv bch_w bch_x 1e-6 "g: basechazard weighted == expanded"
    cmpv mg_w mg_x 1e-6 "h: mgale partial weighted == expanded"
end

* C7: postestimation commands; lrtest and estat ic vs stcox's
capture program drop c7
program define c7
    simdata
    quietly stcox x1 x2 i.x3, shared(gid) breslow
    estimates store sF
    quietly stcox x1 i.x3, shared(gid) breslow
    estimates store sR
    quietly lrtest sF sR
    local lrs = r(chi2)
    stcoxem x1 x2 i.x3, shared(gid) nolog
    estimates store eF
    quietly estat ic
    tempname S
    matrix `S' = r(S)
    local aic = `S'[1,5]
    quietly stcoxem x1 i.x3, shared(gid) nolog
    estimates store eR
    quietly lrtest eF eR
    cmps r(chi2) `lrs' 1e-3 "a: lrtest chi2 vs stcox"
    estimates restore sF
    quietly estat ic
    matrix `S' = r(S)
    cmps `aic' `S'[1,5] 1e-6 "b: AIC vs stcox"
    estimates restore eF
    rcis "c: estat phtest" 0 estat phtest
    rcis "d: estat vce" 0 estat vce
    rcis "e: estat summarize" 0 estat summarize
    rcis "f: test" 0 test x1 x2
    rcis "g: lincom" 0 lincom x1 + x2, hr
    rcis "h: nlcom" 0 nlcom exp(_b[x1])
    rcis "i: margins" 0 margins x3
    rcis "j: estimates table" 0 estimates table eF eR, stats(theta se_theta ll)
    rcis "k: estat gofplot" 0 estat gofplot
    capture noisily graph export "$TDIR/stcoxem_cert_gofplot.png", replace width(1200)
    if !_rc global PNGS $PNGS stcoxem_cert_gofplot.png
    rcis "l: stcurve" 0 stcurve, survival at(x2=0) at(x2=1)
    capture noisily graph export "$TDIR/stcoxem_cert_stcurve.png", replace width(1200)
    if !_rc global PNGS $PNGS stcoxem_cert_stcurve.png
    graph drop _all
end

* C8: error codes
capture program drop c8
program define c8
    simdata
    rcis "a: shared() missing" 198 stcoxem x1
    rcis "b: vce(foo)" 198 stcoxem x1, shared(gid) vce(foo)
    rcis "c: vce(cluster) no variable" 198 stcoxem x1, shared(gid) vce(cluster)
    rcis "d: groups not nested in clusters" 498 stcoxem x1, shared(gid) vce(cluster x3)
    rcis "e: nolog with dots" 198 stcoxem x1, shared(gid) nolog dots
    rcis "f: reps(1)" 198 stcoxem x1, shared(gid) vce(bootstrap) reps(1)
    rcis "g: unknown option" 198 stcoxem x1, shared(gid) foo
    quietly stcoxem x1, shared(gid) nolog
    quietly replace x1 = x1 + 1
    rcis "h: predict after data change" 459 predict double z, mgale
    quietly regress x1 x2
    rcis "i: replay after regress" 301 stcoxem
    quietly stset y [iweight=w], failure(d) id(id) time0(t0)
    rcis "j: iweight from stset" 101 stcoxem x1, shared(gid)
    quietly stset, clear
    rcis "k: not stset" 119 stcoxem x1, shared(gid)
end

* C9 (optional): registry test set vs stcox, shared()
capture program drop c9
program define c9
    capture confirm file "$REALDATA"
    if _rc {
        di as txt "  (skipped: $REALDATA not found)"
        exit
    }
    use "$REALDATA", clear
    keep ID followup_years outcome drug GENDER DM age_entry
    quietly keep if ID<=200000
    quietly stset followup_years, enter(time 0) exit(time .) failure(outcome==1) id(ID)
    quietly stcox i.drug i.GENDER i.DM age_entry, shared(ID) breslow
    tempname bs
    matrix `bs' = e(b)
    local th = e(theta)
    local ll = e(ll)
    stcoxem i.drug i.GENDER i.DM age_entry, shared(ID) nolog
    tempname be
    matrix `be' = e(b)
    cmpm `be' `bs' 1e-5 "a: b vs stcox"
    cmps e(theta) `th' 1e-4 "b: theta vs stcox"
    cmps e(ll) `ll' 1e-6 "c: ll vs stcox"
end

* C10: delayed entry -- as stcox: 498 without forceshared
capture program drop c10
program define c10
    simdata late
    rcis "a: delayed entry, no forceshared" 498 stcoxem x1 x2 i.x3, shared(gid) nolog
    quietly stcox x1 x2 i.x3, shared(gid) breslow forceshared
    tempname bs
    matrix `bs' = e(b)
    local th = e(theta)
    local ll = e(ll)
    quietly predict double st_eff, effects
    quietly predict double st_mg, mgale
    stcoxem x1 x2 i.x3, shared(gid) nolog forceshared $TIGHT
    tempname be
    matrix `be' = e(b)
    cmpm `be' `bs' 1e-5 "b: b vs stcox forceshared"
    cmps e(theta) `th' 1e-4 "c: theta vs stcox forceshared"
    cmps e(ll) `ll' 1e-6 "d: ll vs stcox forceshared"
    quietly predict double em_eff, effects
    quietly predict double em_mg, mgale
    cmpv em_eff st_eff 1e-3 "e: effects vs stcox forceshared"
    cmpv em_mg st_mg 1e-4 "f: mgale vs stcox forceshared"
end

* C0 (informational, never fails): EM iterations to convergence
capture program drop c0
program define c0
    foreach dat in sim real {
        if "`dat'" == "sim" {
            simdata
            local spec x1 x2 i.x3, shared(gid)
        }
        else {
            capture confirm file "$REALDATA"
            if _rc continue
            use "$REALDATA", clear
            keep ID followup_years outcome drug GENDER DM age_entry
            quietly keep if ID<=200000
            quietly stset followup_years, enter(time 0) exit(time .) failure(outcome==1) id(ID)
            local spec i.drug i.GENDER i.DM age_entry, shared(ID)
        }
        quietly stcoxem `spec' vce(none) nolrtest $TIGHT
        local thT = e(theta)
        di as txt _n "  `dat': theta at tolerance(1e-11) = " %18.0g `thT' "  iter = " e(iter)
        di as txt "  " %-9s "accel" %10s "tol" %8s "iter" %6s "conv" %10s "sq att" %8s "sq acc" ///
            %14s "theta reldif"
        foreach acc in squarem aitken none {
            foreach tol in 1e-6 1e-7 1e-8 {
                quietly stcoxem `spec' vce(none) nolrtest maxiter(5000) ///
                    tolerance(`tol') accelerate(`acc')
                di as txt "  " %-9s "`acc'" %10s "`tol'" %8.0f e(iter) %6.0f e(converged) ///
                    %10.0f e(squarem_att) %8.0f e(squarem_acc) %14.2e reldif(e(theta), `thT')
            }
        }
        quietly stcoxem `spec' nolrtest
        di as txt "  default fit: iter " e(iter) " converged " e(converged) ///
            "  plse nonconverged refits " e(theta_profile_nonconv)
    }
end

* C11: every acceleration converges to the same solution (maxiter 2000)
capture program drop c11
program define c11
    foreach dat in sim simlate real {
        if "`dat'" == "sim" {
            simdata
            local spec x1 x2 i.x3, shared(gid)
        }
        else if "`dat'" == "simlate" {
            simdata late
            local spec x1 x2 i.x3, shared(gid) forceshared
        }
        else {
            capture confirm file "$REALDATA"
            if _rc continue
            use "$REALDATA", clear
            keep ID followup_years outcome drug GENDER DM age_entry
            quietly keep if ID<=200000
            quietly stset followup_years, enter(time 0) exit(time .) failure(outcome==1) id(ID)
            local spec i.drug i.GENDER i.DM age_entry, shared(ID)
        }
        quietly stcoxem `spec' vce(none) nolrtest $TIGHT
        local thT = e(theta)
        tempname bT
        matrix `bT' = e(b)
        foreach acc in squarem aitken none {
            quietly stcoxem `spec' vce(none) nolrtest accelerate(`acc') maxiter(2000)
            tempname bA
            matrix `bA' = e(b)
            di as txt "  `dat' `acc': iter " e(iter) " converged " e(converged)
            check "    converged" e(converged) == 1
            cmps e(theta) `thT' 1e-5 "    `dat' `acc': theta vs tight"
            cmpm `bA' `bT' 1e-5 "    `dat' `acc': b vs tight"
        }
    }
end

* ======================= run ============================================
local tl = cond("$TONLY" == "", "0 1 2 3 4 5 6 7 8 9 10 11", "$TONLY")
foreach k of local tl {
    runtest c`k'
}

local rc 0
if "$FAILED" != "" {
    di as err _n "FAILED: $FAILED"
    local rc 9
}
else di as res _n "STCOXEM 1.0.0 CERTIFICATION PASSED"
log close tlog
capture noisily ghlog using "`logf'", repo(szabobarnaKI/claude_communication) path(logs/stcoxem/)
foreach f of global PNGS {
    capture noisily ghlog using "$TDIR/`f'", repo(szabobarnaKI/claude_communication) path(logs/stcoxem/)
}
adopath - "$TDIR"
exit `rc'
