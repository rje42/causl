library(causl)
# basic example
## the parameters are customizable - for example, the model families, the regression params. Here we just assume they are fixed.

gen_causl_example1 <-function(n){
    # formulae corresponding to covariates, treatments, outcomes and the dependence
    forms <- list(Z ~ 1, X ~ Z, Y ~ X, ~ 1)
    # vector of model families (3=gamma/exponential, 1=normal/Gaussian)
    fam <- c(3, 1, 1, 1)
    # list of parameters, including 'beta' (regression params) and 'phi' dispersion
    pars <- list(Z = list(beta=log(1), phi=1),   # note log-link
                X = list(beta=c(0,0.5), phi=1),
                Y = list(beta=c(-0.5,0.5), phi=1),
                cop = list(beta=1))

    ## now create a `causl_model` object
    cm <- causl_model(formulas=forms, family=fam, pars=pars)

    # now simulate 100 observations
    df <- rfrugal(n=n, causl_model=cm)
    return(df)
}


# an example of simulating data from causl with more complicated designs
gen_causl_example2 <- function(n=10000, nI=3, nX=1, nO=1, nS=1, ate=2, beta_cov=0, strength_instr=3, strength_conf=1, strength_outcome=0.2){
    O_terms <- paste0("O", 1:nO)
    O_sum <- paste(O_terms, collapse = " + ")

    X_terms <- paste0("X",1:nX)
    X_sum <- paste(X_terms, collapse = " + ")

    # Construct the formula string with the new sin() term
    formula_str <- paste0(
    "Y ~ A + ", X_sum, "+", O_sum
    )
    po_formula <- as.formula(formula_str)
    forms <- list(list(), A ~ 1, po_formula, ~ 1)

    fam <- list(rep(1, nI + nX + nO + nS), 5, 1, 1)

    pars <- list()

    # Specify the formula and parameters for each covariate type
    ## Instrumental variables (I)
    if (nI > 0) {
    for (i in seq_len(nI)) {
        forms[[1]] <- c(forms[[1]], as.formula(paste0("I", i, " ~ 1")))
        pars[paste0("I", i)] <- list(list(beta = beta_cov, phi = 1))
    }
    }

    ## Confounders (X)
    if (nX > 0) {
    for (i in seq_len(nX)) {
        forms[[1]] <- c(forms[[1]], as.formula(paste0("X", i, " ~ 1")))
        pars[paste0("X", i)] <- list(list(beta = beta_cov, phi = 1))
    }
    }

    ## Outcome variables (O)
    if (nO > 0) {
    for (i in seq_len(nO)) {
        forms[[1]] <- c(forms[[1]], as.formula(paste0("O", i, " ~ 1")))
        pars[paste0("O", i)] <- list(list(beta = beta_cov, phi = 1))
    }
    }

    ## Spurious variables (S)
    if (nS > 0) {
    for (i in seq_len(nS)) {
        forms[[1]] <- c(forms[[1]], as.formula(paste0("S", i, " ~ 1")))
        pars[paste0("S", i)] <- list(list(beta = beta_cov, phi = 1))
    }
    }

    # Specify the formula for A given covariates
    ## Add I to the propensity score formula
    if (nI > 0) {
    for (i in seq_len(nI)) {
        forms[[2]] <- update.formula(forms[[2]], paste0("A ~ . + I", i))
    }
    }

    ## Add X to the propensity score formula
    if (nX > 0) {
    for (i in seq_len(nX)) {
        forms[[2]] <- update.formula(forms[[2]], paste0("A ~ . + X", i))
    }
    }

    # Parameters for copula
    parY <- list()
    parY_names <- c()

    if (nX > 0) {
    parY <- c(parY, rep(list(list(beta = 0)), nX))
    parY_names <- c(parY_names, paste0("X", seq_len(nX)))
    }
    if (nO > 0) {
    parY <- c(parY, rep(list(list(beta = 0)), nO))
    parY_names <- c(parY_names, paste0("O", seq_len(nO)))
    }
    if (nI > 0) {
    parY <- c(parY, rep(list(list(beta = 0)), nI))
    parY_names <- c(parY_names, paste0("I", seq_len(nI)))
    }
    if (nS > 0) {
    parY <- c(parY, rep(list(list(beta = 0)), nS))
    parY_names <- c(parY_names, paste0("S", seq_len(nS)))
    }

    names(parY) <- parY_names
    pars$cop <- list(Y = parY)


    # Set parameters for A
    pars$A$beta <- c(0, rep(strength_instr, nI), rep(strength_conf, nX))


    # Set parameters for Y
    pars$Y$beta <- c(0, ate,  rep(strength_conf, nX),  rep(strength_outcome, nO))
    pars$Y$phi <- 1


    cm <- causl_model(formulas=forms, family=fam, pars=pars)

    ### Generate the data
    # Generate data
    df <- rfrugalParam(n = n, causl_model=cm)
    p <- nX + nI + nO + nS
    # Flatten the A column
    df$A <- as.vector(df$A)
    # Propensity score
    if (nI + nX == 1) {
    df$propen <- plogis(c(rep(strength_instr, nI), rep(strength_conf, nX)) * df[, 1])
    } else {
    df$propen <- plogis(rowSums(c(rep(strength_instr, nI), rep(strength_conf, nX)) * df[, c(1:(nI + nX))]))
    }

    df$mu0 = rowSums(df[,c((nI+1):(nI+nX+nO))])
    df$mu1 = df$mu0 + ate
    colnames(df) <- c(paste("X", 1:p, sep = ""), 'A', 'y', 'propen','mu0','mu1')

    # # Remove nested attributes -- if you run across prolems
    # attributes(df) <- NULL
    return(df)
}