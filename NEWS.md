causl 0.10.4
-------------------------------------------------------------------------------

BUG FIXES

 * Plasmode method fixed for `rfrugal()`.
 
 * Corrected bug in `check_pars()`, which meant that names of supplied parameters
 were not actually being checked.



causl 0.10.0
-------------------------------------------------------------------------------

NEW FEATURES

 * Sandwich standard error estimation is now performed much more quickly.
 
 * Added `is_discrete()` function for families.

 * `fit_causl()` now initializes correctly with respect to link functions.
 
 * `dfgmCopula()` now has `log` argument.
 
 * `par2` now renamed `df` for t-distributions and copulas.
 
 * `process_family()` and `link_setup()` replaced by `process_family_link()`.


BUG FIXES

 * `fit_causl()` uses inputted link functions.



causl 0.9.10
-------------------------------------------------------------------------------

BUG FIXES

 * Fixed bug in parameter initialization for fitting Gamma distributions. 



causl 0.9.9
-------------------------------------------------------------------------------

BUG FIXES

 * Fixed bugs in simulation of ordinals and documented parameterizations for
 categorical and ordinal variables.

UPDATE

 * Improved plasmode simulation by giving more sophisticated way to obtain 
 conditional empirical quantiles.



causl 0.9.8
-------------------------------------------------------------------------------

BUG FIXES

 * Fixed bugs in testing files



causl 0.9.7
-------------------------------------------------------------------------------

UPDATE

 * Merged changes from `devel` branch.


BUG FIXES

 * Fixed bug for using plasmode data with `rfrugal()`
 
 * Fixed bugs in `pair_copula_setup()`
 
 * Fixed bug in `glm_sim()`



causl 0.9.6
-------------------------------------------------------------------------------

BUG FIX

 * Corrected quantile computation for Bernoulli random variables.



causl 0.9.5
-------------------------------------------------------------------------------

BUG FIX

* Fixed bug in `link_setup` when no links provided



causl 0.9.4
-------------------------------------------------------------------------------

NEW FEATURES

* Allow for fitting of t-distributions



causl 0.9.2
-------------------------------------------------------------------------------

NEW FEATURES

 * We now allow for specification of only certain link functions.
 
 * Added `link()` to obtain link functions from `causl_family` and `causl_cop`
 objects.
 


causl 0.9.1
-------------------------------------------------------------------------------

NEW FEATURES

 * Introduced `causl_model()` and `rfrugal()` functions, to create a causal 
 model---an object that can be used to simulate from or, in future, fit a 
 frugally parameterized model.

 * Added vignette showing how to add a custom family to simulation methods.
 

CHANGES

 * `rfrugalParam` calls new function `rfrugal` to simulate data.

 * Redundant argument `seed` removed from `rfrugalParam`.  Argument `careful` 
 moved to `control`.
 
 * Other `control` arguments removed: `trace`, `warn`, `max_wt`.
 
 * `control` argument removed from `sim_inversion()` and `sim_rejection()`.

 * `"multi_copula"` changed to `"inversion_mv"`.



causl 0.9.0
-------------------------------------------------------------------------------

NEW FEATURES

 * Introducing ability to simulate from a multivariate copula using the 
 inversion method.  Choose `method = "multi_copula"`.



causl 0.8.9
-------------------------------------------------------------------------------

NEW FEATURES

 * Made it possible to use a custom distribution to simulate via an object of 
 class `causl_family`.


CHANGES

 * Began migration of `univarDens()` to `glm_dens()`.

 * Changed name `dGaussDiscCop2()` to `dGaussDiscCop()`.
 
 * Made changes to `fitCausal`, including renaming as `fit_causl` and removing
 some arguments from function and splitting off identification of discrete
 variables.
 
 * Modified `par2` in simulations for t-distribution to `df`.


BUG FIXES

 * Fixed bug in discrete `causal_family` types.
 
 * Fixed bug in `pair_copula_setup()`.
 
 * Added `link` to output from `binomial_causl_fam`.



causl 0.8.8
-------------------------------------------------------------------------------

NEW FEATURES

 * Added log-link for binomial family


CHANGES
 
 * Suppressed messages in tests
 
 * Edited to remove quantiles from simulated data under inversion method



causl 0.8.7
-------------------------------------------------------------------------------

CHANGES

 * Added `process_prespecified`



causl 0.8.6.9000
-------------------------------------------------------------------------------

CHANGES

 * Minor fix to `process_family`

 * New test for plasmode simulation with one pregenerated variable.



causl 0.8.5.9000
-------------------------------------------------------------------------------

CHANGES

 * Fixed bug in plasmode dataset generation that caused failure if only one 
 variable was used.



causl 0.8.4.9000
-------------------------------------------------------------------------------

CHANGES

 * Further reorganization, mostly to export additional functions.

 

causl 0.8.3.9000
-------------------------------------------------------------------------------

CHANGES

 * `process_inputs` and its derivatives further edited for compatibility with 
 `survivl` package.



causl 0.8.2.9000
-------------------------------------------------------------------------------

CHANGES

 * `process_inputs()` function now substantially reorganized.  Uses various 
 subsidiary functions to carry out its work.


BUG FIXES

 * Fixed minor bug in `is_categorical()`



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
