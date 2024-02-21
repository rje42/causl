causl 0.8.1.9000
-------------------------------------------------------------------------------

NEW FEATURES

 * `rfrugalParam()` now allows for 'plasmode' simulation.  That is, one can add
 a `data.frame` of variables `dat` that acts as a collection of pre-generated
 covariates.
 
 * New vignette 'Plasmode' gives instructions on how to use this feature.
 
 

causl 0.8.0.9000
-------------------------------------------------------------------------------

NEW FEATURES

 * Added a `NEWS.md` file!

 * Categorical and ordinal variables are now implemented for simulation with
 `rfrugalParam()` (provided that the now default `method="inversion"` is used).  The
 parameterizations for categorical variables use the logit-contrasts of each
 level with the baseline-level; for ordinal variables it uses the logit of
 successive cumulative probabilities.  In other words, it uses 
 $\log P(X \leq i)/P(X > i)$ for
 each $i = 1,\ldots,l-1$ where there are $l$ levels.
 
 * Added a vectorized `cumsum_mat` function in `Rcpp`, that acts on the rows of 
 a matrix.
 

CHANGES

 * `rfrugalParam()` now has `method="inversion"` as the default.  The much
 slower `method="rejection"` is kept only for backwards compatibility.

 * Some functions (`rescaleVar`, `rescaleCop`, `sim_CopVal`) have had their old 
 names deprecated in favour of snake case versions (`rescale_var`, `rescale_cop`,
 `sim_copula`).  In addition the dataframe `copulaVals` has been replaced by 
 `copula_vals`.
 
 * `link_hack` has been renamed to `link_apply`.
