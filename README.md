# stcoxem

`stcoxem` fits the Cox proportional hazards model with a shared gamma frailty,
the model of `stcox, shared()`, by the EM algorithm of Klein (1992),
implemented in Mata for large data. Stata 19 or later.

* Estimates (coefficients, frailty variance theta, frailties, baseline hazard,
  log likelihood) agree with `stcox, shared()`.
* Also allows `strata()`, `vce(robust)`, `vce(cluster)`, `vce(bootstrap)`,
  `vce(none)` (fastest, point estimates only), and pweights (set with `stset`).
* Standard error of theta from the profile likelihood (default), and an
  optional profile-likelihood confidence interval (`plci`).
* Results are stored as `stcox, shared()` stores them, so `predict`, `estat`,
  `stcurve`, `margins`, `test`, `lrtest`, and the other postestimation
  commands of `stcox` work after `stcoxem`.

## Installation

```stata
net install stcoxem, from("https://raw.githubusercontent.com/szabobarnaKI/stcoxem/main/") replace
```

or

```stata
net from "https://szabobarnaki.github.io/stcoxem/"
```

## Example

```stata
webuse catheter
stset time, failure(infect)
stcoxem age female, shared(patient)
predict nu, effects
estat phtest
```

See `help stcoxem` and `help stcoxem postestimation`.

## Files

| File | Content |
|---|---|
| `stcoxem.ado` | the command (with its Mata code) |
| `stcoxem_p.ado` | `predict` after fits with pweights |
| `stcoxem.sthlp` | help |
| `stcoxem_postestimation.sthlp` | help for postestimation |

## Reference

Klein, J. P. 1992. Semiparametric estimation of random effects using the Cox
model based on the EM algorithm. *Biometrics* 48: 795-806.

## Author

Barna Szabo-Soderberg, Karolinska Institutet
