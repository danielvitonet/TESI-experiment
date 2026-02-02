## Description
> This code comes without technical support of any kind. The code is free to use, provided that the paper is cited properly.

The repository contains code files to estimate two different models: MF-BAVART model and flexBART model.
Each model has its own folder with corresponding code files for estimation.
Code files for MF-BAVART reside within the folder `MF-BAVART`.
Code files for flexBART reside within the folder `flexBART`.


### Code files for MF-BAVART
These files create a function mfbavart(...) to estimate the MF-BAVART model.
- `mfbavart_func.R` contains the main function
- `aux_func.R` collects several auxiliary functions
- `example.R` contains an example code for using the function.

In addition to the baseline model, the code also includes an option to introduce stochastic volatility (SV) in the error terms. Several parts of the original code used for the paper in directory "replication" have been replaced to improve computational efficiency.

Some codes and helper functions are taken or adapted from the "mfbvar" package. Thanks to Vincent Dorie (mtn. of "dbarts") and Sebastian Ankargren (mtn. of "mfbvar") for technical support regarding their excellent packages.

### Inputs for mfbavart(...):
- data            a list that contains ts-objects of different frequencies in its 
                  M (number of endogenous variables) slots, such that high-frequency (monthly)
                  series are ordered first and followed by low-frequency series (quarterly)
- itr             intertermporal restriction ("lvl" or "grw") of length 
                  corresponding to number of low frequency series
- p               numeric, lag-length of the VAR (minimum of 5 if itr=="grw", and 3 if itr=="lvl")
- fhorz           numeric, forecast horizon in months (3 per quarter)
- cons            TRUE/FALSE, whether a constant should be included
- exact           TRUE/FALSE, whether BART fit is stored in output values or 
                  filtered data based on approximation
- sv              TRUE/FALSE, whether structural errors feature stochastic volatility
- var.thrsh       numeric, threshold for resampling coefficients by draw (for sampler stability)
- max.count.var   numeric, maximum number of resampling steps
- cgm.level       numeric \in (0,1), \alpha in the paper (probability of terminal node)
- cgm.exp         numeric > 0, \beta in the paper (probability of terminal node)
- sd.mu           numeric, \gamma in the paper
- num.trees       numeric, number of trees for BART, S in the paper
- prior.sig       numeric of length 2, [1] nu_j, [2] v in the paper, 
- nburn           numeric, number of burnins
- nsave           numeric, number of draws for posterior/predictive inference
- thinfac         numeric, thinning factor
- quiet           TRUE/FALSE, whether progress bar is indicated during sampling

Function returns:
- Y              [nsave,T,M] array that contains the latent states (see also option "exact")
- fcst           [nsave,fhorz,M] array that contains forecasts
- Yq             [nsave,T+fhorz,M] array that contains aggregated filtered and forecasted series



### Code files for flexBART
These files create a function flexBART(...) to estimate flexBART model plus .
The folder `estimation_functions` contains the estimation functions for different models:
- `flexBART.R` contains the main function for flexBART
- `qr.R` contains the main function for bayesian quantile regression
The rest of the files are:
- `aux_func.R` collects several auxiliary functions
- `fcst_script.R` contains an example code for using the function.

### Inputs for flexBART(...)
- Yraw            numeric matrix with time observations in rows and variables in columns (T × M).
                  standardized in-place using column means/sds, with na.rm=TRUE allowed, and then used to form the lagged design matrix and response matrix. 
                  This implies:
                  - Rows = time points.
                  - Columns = multiple series.
                  - NAs are allowed (used with na.rm for means/sds).
                  Evidence: standardization and creation of X/Y rely on matrix structure and columnwise means/sds.

- nburn           scalar integer (number of burn-in MCMC draws).
                  burn-in length and BART control parameter n.burn.

- nsave           scalar integer (number of saved draws).
                  defines total samples ntot, thinning count nthin, and BART control n.samples.

- thinfac         (default 1) numeric scalar; thinning multiplier.
                  Used in nthin <- round(thinfac * nsave) to select how many draws are stored.

- prior           string flag; the code checks for "Minn" to switch to a Minnesota prior branch.
                  Used as "Minn" triggers Minnesota prior logic for VAR coefficients. Otherwise HS prior machinery is used (default).

- prior.sig       numeric vector of length 2 (shape and scale parameters) for the BART residual prior (chisq).
                  Used as chisq(prior.sig[[1]], prior.sig[[2]]) for each equation’s BART residual prior (unless overridden by sv). In sv == "SV" or "heteroBART" it is internally overridden to c(10000^50, 0.5).

- model           (default "mixBART") string flag; "mixBART" or "BART" enables BART-driven mean dynamics; other values imply linear VAR-only behavior.
                  Used as branch for whether tree predictions are used when sampling and forecasting.

- sv              string flag. Recognized values include "SV", "heteroBART", and "homo" (as seen in script usage).
                  Used as controls stochastic volatility modeling and associated priors; "SV" builds stochastic volatility priors; "heteroBART" uses BART to model volatility. This also changes how prior.sig is set. 

- fc.approx       string flag; valid values used in code are "exact" and "approx".
                  Used as switches between exact sampling of the forecast mean (with BART mean + VAR) vs. an approximate form using (A_approx + PHI_draw).

- restr.var       NULL or a column name matching colnames(Y); the code uses which(colnames(Y)==restr.var).
                  Used as if non-NULL, forecasts are conditional on specified variable values across quantile grid range.conditional, which adds an extra R dimension to outputs.

- fhorz           scalar integer; number of forecast steps.
                  Used as forecast horizon dimension in arrays and loop bounds.

- quiet           logical scalar. 
                  Used as if FALSE, shows progress bar and diagnostic plotting during sampling. 

- pr.mean         numeric matrix (M × M) used for the prior mean of the linear VAR coefficients.
                  Used as inserted into A_prior for the first M rows/columns, so it must be compatible with the dimension M = ncol(Y).

The function returns a list with two elements:
list("fcst" = fcst_store, "Hfcst" = Hfcst_store)

- fcst            numeric array of forecast draws.
                  Dimensions:
                  - If restr.var is NULL: c(nthin, fhorz, M).
                  - If restr.var is not NULL: c(nthin, fhorz, M, R) where R is the number of conditional quantile grid points. 
                  Dimnames:
                  - mcmc1..mcmcN for draw index,
                  - fhorz1..fhorzH for horizon,
                  - variable names from colnames(Y).
                  Scale: stored in the original data scale (the standardization is undone via Ysd and Ymu).
                  Unconditional branch rescales explicitly when storing.
                  Conditional branch stores Y.tp1 * Ysd + Ymu before storing.

- Hfcst           numeric array of forecast volatility draws (per series).
                  Dimensions: same as fcst (c(nthin, fhorz, M) or c(nthin, fhorz, M, R)).
                  Dimnames: identical to fcst.
                  Scale: rescaled by Ysd in unconditional branch; conditional branch stores exp(HT) * Ysd.
