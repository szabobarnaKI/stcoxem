*! version 1.0.0  26nov2026  Barna Szabo-Soderberg (stcoxem rev 52)
*! stcoxem_p.ado -- predict after stcoxem with weights (rev 51, 2026-11-26)
*!
*! Set as e(predict) by stcoxem for pweighted fits only (unweighted fits use
*! stcox_p directly). stcox's own refit after shared() would re-solve a
*! weighted penalized fit, which is a different estimator from stcoxem's
*! "weight k = k copies of the group" (probe 2/3, 2026-11-25). Hence:
*!  - effects, atfrailty: stcoxem's own frailties, ln u rebuilt at e(b) and
*!    e(theta) by -stcoxem _lnu-;
*!  - refit statistics (mgale, csnell, deviance, basechazard, basehc,
*!    basesurv, scores/esr, schoenfeld, scaledsch, dfbeta, ldisplace, lmax):
*!    stcox_p itself, run on a plain-Cox view of the results (e(cmd)="cox",
*!    no e(shared), e(offset) = ln u), so its refit is
*!    -cox ... [pweight], matfrom(e(b)) iter(0) offset(ln u)-: the Breslow
*!    baseline and residuals at stcoxem's b, theta and frailties;
*!  - xb, stdp, hr: stcox_p on the results as posted.
*! Anything else (invalid combinations) goes to stcox_p for its own error.

program stcoxem_p, sort
    version 19.0
    if "`e(cmd2)'" != "stcox" | "`e(estimator)'" != "em" {
        error 301
    }
    local zero : copy local 0

    syntax [anything] [if] [in] [, XB Index STDP HR CSNell DEViance ///
        CCSNell CMGale MGale BASEChazard basehc BASESurv EFFects ESR    ///
        SCOres SCHoenfeld SCAledsch DFBeta LDisplace LMax PARTial       ///
        ATFRAILty1 ATFRAILty(string) noOFFset ]

    local type "`xb' `index' `stdp' `hr' `csnell' `deviance' `ccsnell'"
    local type "`type' `cmgale' `mgale' `basechazard' `basehc' `basesurv'"
    local type "`type' `effects' `scores' `esr' `schoenfeld' `dfbeta'"
    local type "`type' `ldisplace' `scaledsch' `lmax'"
    local type : list retokenize type
    if `: word count `type'' > 1 {
        stcox_p `zero'
        exit
    }
    local atf = ("`atfrailty1'" != "" | `"`atfrailty'"' != "")

    * ---- effects: stcoxem's own log frailties ----
    if "`type'" == "effects" & !`atf' {
        if "`partial'`offset'" != "" {
            stcox_p `zero'
            exit
        }
        Effects `anything' `if' `in'
        exit
    }

    * ---- atfrailty / atfrailty(): base statistic times the frailty, as
    *      stcox_p (GenSinglePred), with stcoxem's own frailties ----
    if `atf' {
        if ("`atfrailty1'" != "" & `"`atfrailty'"' != "") | ///
           !inlist("`type'", "basesurv", "basehc", "basechazard") | ///
           "`partial'`offset'" != "" {
            stcox_p `zero'
            exit
        }
        local afopt `atfrailty1'
        if `"`atfrailty'"' != "" local afopt atfrailty(`atfrailty')
        AtFrail `anything' `if' `in', type(`type') `afopt'
        exit
    }

    * ---- statistics without a refit: stcox_p on the posted results ----
    local norefit xb index stdp hr
    if "`type'" == "" | `: list type in norefit' {
        stcox_p `zero'
        exit
    }

    * ---- refit statistics: stcox_p on a plain-Cox view with offset ln u
    *      (held results come back when this program ends, also on error)
    tempvar lnu
    quietly stcoxem _lnu `lnu'
    tempname hold
    _estimates hold `hold', copy restore
    MaskCox `lnu'
    stcox_p `zero'
end


program MaskCox, eclass
    args lnu
    ereturn local cmd     "cox"
    ereturn local shared  ""
    ereturn local offset  "`lnu'"
    ereturn local predict "stcox_p"
end


program Effects
    syntax newvarname [if] [in]
    marksample touse, novarlist
    tempvar lnu
    quietly stcoxem _lnu `lnu'
    if e(theta) < 1e-12 {
        gen `typlist' `varlist' = . if `touse'
    }
    else {
        gen `typlist' `varlist' = `lnu' if `touse'
    }
    label var `varlist' "Shared log-frailty"
end


program AtFrail
    syntax newvarname [if] [in], type(string) [ ATFRAILty1 ATFRAILty(string) ]
    marksample touse, novarlist
    tempvar frailty tvar
    if "`atfrailty1'" != "" {
        tempvar lnu
        quietly stcoxem _lnu `lnu'
        quietly gen double `frailty' = exp(`lnu')
    }
    else {
        * validation as stcox_p
        capture confirm number `atfrailty'
        if _rc {
            capture confirm numeric variable `atfrailty'
            if _rc {
                di as err "{p}a positive number or a "
                di as err "numeric variable expected in "
                di as err "option {bf:atfrailty()}{p_end}"
                exit 198
            }
            quietly summarize `atfrailty', meanonly
            if r(min) <= 0 {
                di as err "{p}only positive values allowed " ///
                    "in variable {bf:`atfrailty'}{p_end}"
                exit 198
            }
        }
        else if `atfrailty' <= 0 {
            di as err "{p}a positive number or a "
            di as err "numeric variable expected in "
            di as err "option {bf:atfrailty()}{p_end}"
            exit 198
        }
        quietly gen double `frailty' = `atfrailty'
    }
    quietly predict double `tvar', `type'
    if "`type'" == "basesurv" {
        gen `typlist' `varlist' = `tvar'^`frailty' if `touse'
    }
    else if "`type'" == "basehc" {
        gen `typlist' `varlist' = 1-(1-`tvar')^`frailty' if `touse'
    }
    else {
        gen `typlist' `varlist' = `tvar'*`frailty' if `touse'
    }
    quietly levelsof `frailty'
    if r(r) == 1 local fr_note " at specified frailty"
    else         local fr_note " at specified frailties"
    if "`type'" == "basechazard" {
        label var `varlist' "Baseline cumulative hazard`fr_note'"
    }
    else if "`type'" == "basehc" {
        label var `varlist' "Baseline hazard contribution`fr_note'"
    }
    else {
        label var `varlist' "Baseline survivor function`fr_note'"
    }
end
