# stcoxem
stcoxem fits a Cox proportional hazards model with gamma frailty using the Expectation-Maximization method in Stata.

stcoxem fits a Cox proportional hazards model with gamma-distributed frailty using the Expectation-Maximization method.
The E-M method provides significant speed benefit over the Stata built-in stcox command, but preserves accuracy.

To install the command, run the following command in Stata:
##
    net from "https://szabobarnaki.github.io/stcoxem/"
