##' List of link functions for copuulas
##'
##' @param family either the name of the family or its integer representative
##'
##' @details
##' This function returns the default link function for each possible copula.
##' These are given in the table:
##'
##' |value|family     |link                  |  link name |
##' |-----:|----------|----------------------|---|
##' |1    |gaussian  |logit((1+rho)/2)      | logit2 |
##' |2    |t         |logit((1+rho)/2)      |logit2 |
##' |3    |Clayton   |log(1+theta)      | log1p |
##' |4    |Gumbel    |log(theta - 1)      | log1m |
##' |5    |Frank     |logit(1 + theta)      | log1p |
##' |6    |Joe         |identity      | identity |
##' |11    |FGM         |logit((1+rho)/2)      | logit2 |
##'
##'
links_cop <- function (family) {
  switch (family,
          gaussian="logit2",  # 1
          t="logit2",  # 1
          Clayton="log1p",   # 3
          Gumbel="log1m",   # 4
          Frank="log1p",   # 5
          Joe="identity", # 6
          NULL, NULL, NULL, NULL,
          FGM="logit2"  # 11
  )
}

##' Simulate from vine copula
##'
##' @param dat data frame to be filled in
##' @param family family to simulate from
##' @param par matrix of parameters
##' @param par2 extra parameters for t-copula
##' @param model_matrix design matrix for covariates
##' @param link link functions for parameters (currently unused)
##'
##' @return A data frame of the same dimensions as \code{dat}.
##'
##' @export
sim_vinecop <- function (dat, family, par, par2=NULL, model_matrix, link) {

  k <- ncol(dat)
  if (length(family) != choose(k, 2)) stop("invalid family vector")

  ## get link functions
  if (!missing(link)) {
    lnks <- lapply(link, make.link)
  }
  else {
    lnks <- lapply(family, links_cop)
  }

  ## check on dimensions of model_matrix
  if (!missing(model_matrix)) {
    if (nrow(dat) != nrow(model_matrix)) stop("dat and model_matrix must have same number of rows")
    if (nrow(par) != ncol(model_matrix)) stop("rows of par and columns of model_matrix must be the same")
  }
  else {
    model_matrix <- matrix(1, nrow(dat), 1)
    if (nrow(par) != 1) stop("if model_matrix not supplied, par must have only one row")
  }

  # set up vine model
  if (missing(par2) && any(family == 2)) stop("Must specify degrees of freedom for t-copula")
  else if (any(family == 2)) {
    par2_mat <- matrix(0, k, k)
    par2_mat[lower.tri(par2_mat)] <- par2
  }
  else par2_mat <- matrix(0, k, k)

  vine <- fam_mat <- par_mat <- matrix(0, k, k)
  vine[lower.tri(vine, diag = TRUE)] <- rev(unlist(lapply(seq_len(k), seq_len)))
  err_code <- VineCopula::RVineMatrixCheck(vine)
  if (err_code != 1) stop(paste("VineCopula error code", err_code))
  fam_mat[lower.tri(fam_mat)] <- family

  cor_val <- rep(0, choose(k,2))

  ## now, if there is a covariate, get levels of all the covariates
  if (nrow(unique(model_matrix)) > 1) {
    facts <- as.matrix(apply(model_matrix, 2, function(x) as.numeric(factor(x)))) - 1
    levs <- apply(facts, 2, max)+1
    b <- cumprod(c(1, levs[-length(levs)]))
    ovl_fact <- facts %*% b
    fact <- factor(ovl_fact)
  }
  else {
    fact <- factor(rep(1,nrow(dat)))
  }

  ## now go though each group and get a suitable sample
  for (gp in levels(fact)) {
    gp_mem <- (fact == gp)
    n_gp <- sum(gp_mem)
    if (n_gp == 0) next

    eta <- c(model_matrix[match(TRUE, gp_mem), ,drop=FALSE] %*% par)

    ## transform variables accordingly
    cor_val[family == 1 | family == 2] <- 2*expit(eta[family == 1 | family == 2]) - 1  # logit2
    cor_val[family == 3] <- exp(eta[family == 3]) - 1  # log1p
    cor_val[family == 4 | family == 6] <- exp(eta[family == 4 | family == 6]) + 1  # log1m
    cor_val[family == 5] <- eta[family == 5]  # identity

    # par2[lower.tri(par2)] <- c(0,3,3)
    ## set parameters
    par_mat[lower.tri(par_mat)] <- rev(cor_val)

    RVM <- VineCopula::RVineMatrix(Matrix = vine, family = fam_mat,
                                   par = par_mat, par2 = par2_mat,
                                   names = names(dat))
    dat[gp_mem,] <- VineCopula::RVineSim(n_gp, RVM)
  }
  dat
}
