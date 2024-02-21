##' Process formulas, families and parameters
##'
##' @inheritParams rfrugalParam
##' @param kwd keyword for copula
##' @param ordering logical: should an ordering of variables be computed?
## @param ... dots from rfrugalParam
##'
process_inputs <- function (formulas, pars, family, link, dat, kwd, ordering=FALSE) {

  # ## check list of formula
  # if (ordering && !is.numeric(family[[1]]) && ("formula" %in% class(formulas[[1]]))) {
  #   process_inputs2(formulas, pars, family, link, kwd)
  # }


  ## MODIFY CODE TO TAKE IN PRE-SIMULATED VARIABLES

  ## check we have four groups of formulas
  if (length(formulas) != 4) stop("Formulas must have length 4")
  if (missing(pars)) stop("Must supply parameter values")

  ## ensure all formulas are given as lists
  if (any(sapply(formulas, class) == "formula")) {
    wh <- which(sapply(formulas, class) == "formula")
    for (i in wh) formulas[[i]] <- list(formulas[[i]])
  }

  dims <- lengths(formulas[1:3])

  if (missing(family)) {
    ## assume everything is Gaussian
    family <- lapply(lengths(formulas), rep.int, x=1)
  }
  else if (is.list(family)) {
    if (any(sapply(family, class) == "causl_family")) {
      wh <- which(sapply(family, class) == "causl_family")
      for (i in wh) family[[i]] <- list(family[[i]])
    }
    if (!all(lengths(family[1:3]) == lengths(formulas[1:3]))) stop("Mismatch in family and formulae specifications")
  }
  else if (length(family) == 4) {
    if (sum(lengths(formulas[1:3])) > 3) stop("Mismatch in family and formulae specification")
    family <- as.list(family)
  }
  else stop("family should be a list, or vector of length 4")


  # add_args <- list(...)
  # if (length(add_args) > 0) {
  #   nms <- names(add_args)
  # }

  ## check families are valid
  # if (!all(unlist(family[1:3]) %in% familyVals$val)) stop("Invalid family specification")
  if (!all(unlist(family[[4]]) %in% copula_vals$val)) stop("Invalid copula specification")

  ## check that supplied regression parameters are sufficient
  nm_pars <- names(pars)
  wh_cop <- which(nm_pars == kwd)
  if (is.na(wh_cop)) stop("No parameters specified for copula")
  nms <- lapply(pars[-wh_cop], names)
  bpres <- sapply(nms, function(x) "beta" %in% x)
  if (!all(bpres)) {
    plur <- sum(!bpres) > 1
    stop(paste(ifelse(plur, "Variables", "Variable"), paste(names(pars)[!bpres], collapse=", "), ifelse(plur, "lack", "lacks"), "a beta parameter vector"))
  }
  ### change to use the names supplied by each family

  # if (all(unlist(family) == 0)) {
  #   message("Perhaps better to simulate this using the MLLPs package")
  # }

  dZ <- dims[1]
  dX <- dims[2]
  dY <- dims[3]

  famZ <- family[[1]]
  famX <- family[[2]]
  famY <- family[[3]]
  famCop <- family[[4]]

  ## check variable names
  LHS_Z <- lhs(formulas[[1]])
  LHS_X <- lhs(formulas[[2]])
  LHS_Y <- lhs(formulas[[3]])
  rep_pars <- match(c(LHS_Z,LHS_X,LHS_Y), nm_pars, nomatch = 0L)
  if (any(rep_pars == 0)) {
    wh_nrep <- c(LHS_Z,LHS_X,LHS_Y)[rep_pars==0]
    stop(paste0("Variable ", paste(wh_nrep, collapse=", "), "not represented in the parameters list"))
  }

  formsZ <- lapply(formulas[[1]], terms)
  formsX <- lapply(formulas[[2]], terms)
  formsY <- lapply(formulas[[3]], terms)
  if (is.list(formulas[[4]]) && length(formulas[[4]]) > 0) {
    ## put in code to check that copulae formulas and parameters have matching lengths
  }

  ## useful variable summaries
  output <- c(LHS_Z, LHS_Y)
  vars <- c(LHS_Z, LHS_X, LHS_Y)
  LHSs <- list(LHS_Z, LHS_X, LHS_Y)

  ## produce dummy data.frame to check number of coefficients
  dummy_dat <- as.data.frame(rep(list(NA), sum(dims)))
  names(dummy_dat) <- vars
  if (!is.null(dat)) {
    dummy_dat <- cbind(dat[1,], dummy_dat)
  }

  for (i in 1:3) for (j in seq_along(LHSs[[i]])) {
    if (!is_categorical(family[[i]][j])) dummy_dat[[LHSs[[i]][j]]] <- 0
    else dummy_dat[[LHSs[[i]][j]]] <- factor(x=1L, levels=seq_len(pars[[LHSs[[i]][j]]]$nlevel))
  }

  for (i in seq_along(formulas[[1]])) {
    mod_mat <- model.matrix(formulas[[1]][[i]], data=dummy_dat)
    # npar <- length(attr(formsZ[[i]], "term.labels")) + attr(formsZ[[i]], "intercept")
    npar <- ncol(mod_mat)
    if (is_categorical(family[[1]][i])) {
      npar <- npar*(pars[[LHS_Z[i]]]$nlevel - 1)
      if (!is.matrix(pars[[LHS_Z[i]]]$beta)) pars[[LHS_Z[i]]]$beta <- matrix(pars[[LHS_Z[i]]]$beta, nrow=ncol(mod_mat))
    }
    if (isTRUE(is.na(npar)) || length(npar) != 1) stop(paste0("Categorical variable '", LHS_Z[i], "' requires an 'nlevel' parameter"))
    if (length(pars[[LHS_Z[i]]]$beta) != npar) stop(paste0("dimension of model matrix for ", LHS_Z[i], " does not match number of coefficients provided"))
  }
  for (i in seq_along(formulas[[2]])) {
    mod_mat <- model.matrix(formulas[[2]][[i]], data=dummy_dat)
    # npar <- length(attr(formsX[[i]], "term.labels")) + attr(formsX[[i]], "intercept")
    npar <- ncol(mod_mat)
    if (is_categorical(family[[2]][i])) {
      npar <- npar*(pars[[LHS_X[i]]]$nlevel - 1)
      if (!is.matrix(pars[[LHS_X[i]]]$beta)) pars[[LHS_X[i]]]$beta <- matrix(pars[[LHS_X[i]]]$beta, nrow=ncol(mod_mat))
    }
    if (isTRUE(is.na(npar)) || length(npar) != 1) stop(paste0("Categorical variable '", LHS_X[i], "' requires an 'nlevel' parameter"))
    if (length(pars[[LHS_X[i]]]$beta) != npar) stop(paste0("dimension of model matrix for '", LHS_X[i], "' does not match number of coefficients provided"))
  }
  for (i in seq_along(formulas[[3]])) {
    mod_mat <- model.matrix(formulas[[3]][[i]], data=dummy_dat)
    # npar <- length(attr(formsY[[i]], "term.labels")) + attr(formsY[[i]], "intercept")
    npar <- ncol(mod_mat)
    if (is_categorical(family[[3]][i])) {
      npar <- npar*(pars[[LHS_Y[i]]]$nlevel - 1)
      if (!is.matrix(pars[[LHS_Y[i]]]$beta)) pars[[LHS_Y[i]]]$beta <- matrix(pars[[LHS_Y[i]]]$beta, nrow=ncol(mod_mat))
    }
    if (isTRUE(is.na(npar)) || length(npar) != 1) stop(paste0("Categorical variable '", LHS_Y[i], "' requires an 'nlevel' parameter"))
    if (length(pars[[LHS_Y[i]]]$beta) != npar) stop(paste0("dimension of model matrix for ", LHS_Y[i], " does not match number of coefficients provided"))
  }

  ## process family variable inputs
  family <- process_family(family)

  ## check that families have necessary parameters
  for (j in 1:3) {
    if (is.numeric(family[[j]])) {
      for (i in seq_along(family[[j]])) {
        if (family[[j]][i] <= 4 && is.null(pars[[LHSs[[j]][i]]]$phi))
          stop(paste(LHSs[[j]][i], "needs a phi parameter"))
        if (family[[j]][i] == 2 && is.null(pars[[LHSs[[j]][i]]]$par2))
          stop(paste(LHSs[[j]][i], "needs a par2 parameter for degrees of freedom"))
        if (family[[j]][i] %in% 10:11 && is.null(pars[[LHSs[[j]][i]]]$nlevel))
          stop(paste(LHSs[[j]][i], "needs a parameter for number of categories"))
      }
    }
    else if (length(family[[j]]) > 0 && is(family[[j]][1], "causl_family")) {
      for (i in seq_along(family[[j]])) {
        if (!all(family[[j]][i]$pars %in% names(pars[[LHSs[[j]][i]]]))) {
          miss <- family[[j]][i]$pars[!family[[j]][i]$pars %in% names(pars[[LHSs[[j]][i]]])]
          stop(paste0("Parameters ", paste(miss, collapse=", "), " missing for variable ", LHSs[[j]][i]))
        }
      }
    }
  }


  ## check which variables are already in the formula for another, and set
  ## corresponding copula correlation to zero
  # formsZ <- merge_formulas(formulas[[1]])$old_forms
  # setmatch(LHS_Z, rhs_vars(formsZ))
  # formsY <- merge_formulas(formulas[[3]])$old_forms
  if (ordering) {
    ord_mat <- matrix(0, nrow=sum(dims), ncol=sum(dims))

    ## get adjacency matrix of dependencies
    if (dZ > 0) ord_mat[seq_len(dZ), ] <- 1*t(sapply(formsZ, function(x) vars %in% attr(x, "term.labels")))
    if (dX > 0) ord_mat[dZ + seq_len(dX), ] <- 1*t(sapply(formsX, function(x) vars %in% attr(x, "term.labels")))
    if (dY > 0) ord_mat[dZ + dX + seq_len(dY), ] <- 1*t(sapply(formsY, function(x) vars %in% attr(x, "term.labels")))

    if (!is.null(formulas[[4]][LHS_Y])) {
      ## look for possible loops in copula formulae
      cop_frms <- formulas[[4]][LHS_Y]
      if (is.list(cop_frms) && length(cop_frms) > 0 && is.list(cop_frms[[1]])) {
        pred <- lapply(cop_frms, function (x) sapply(x, lhs))
        ord_mat[dZ + dX + seq_len(dY), ] <- ord_mat[dZ + dX + seq_len(dY), ] +
          1*t(sapply(pred, function(x) vars %in% x))
        ord_mat[] <- pmin(ord_mat, 1)
      }
    }

    ord_mat[dZ + dX + seq_len(dY), ] <- ord_mat[dZ + dX + seq_len(dY), ] +
      1*t(sapply(formsY, function(x) vars %in% attr(x, "term.labels")))

    order <- topOrd(ord_mat)
    if (any(is.na(order))) stop("Formulae contain cyclic dependencies")
  }
  else {
    if (any(unlist(lapply(family, is_categorical)))) stop("Categorical variables must be simulated using 'method=\"inversion\"'")

    ## more restrictive for the rejection sampling method
    if (any(unlist(lapply(rhs_vars(formulas[[3]]),
                          function(x) any(x %in% LHS_Z))))) stop("Covariates cannot be direct predictors for outcomes under rejection sampling")
    if (any(unlist(lapply(rhs_vars(formulas[[1]]),
                          function(x) any(x %in% LHS_Y))))) stop("Outcomes cannot be direct predictors for covariates under rejection sampling")
    if (any(unlist(lapply(rhs_vars(formulas[[3]]),  # could relax this
                          function(x) any(x %in% LHS_Y))))) stop("Outcomes cannot directly predict one another under rejection sampling")

    ## check if covariates predict each other
    if (dZ > 1) {
      idx <- lapply(rhs_vars(formulas[[1]]), function(x) na.omit(match(x,LHS_Z)))
      # idx <- lapply(formsZ,
      #               function(x) match(intersect(LHS_Z, attr(x, "term.labels")[attr(x, "order") == 1]), LHS_Z))
      ignr <- do.call(cbind, mapply(function(x,y) {
        if (length(x) > 0) rbind(x,y)
        else matrix(NA, 2, 0)
      }, idx, seq_along(LHS_Z)))
      if (length(ignr) > 0) {
        wh_swt <- ignr[1,] > ignr[2,]
        ignr[,wh_swt] <- ignr[2:1, wh_swt]
        in_form <- setmatch(apply(ignr, 2, c, simplify = FALSE),
                            combn(length(LHS_Z), 2, simplify = FALSE))
      }
      else in_form <- numeric(0)
    }
    else in_form <- numeric(0)

    order <- NULL
  }

  #
  # if (any(sapply(formsZ,
  #        function(x) any(attr(x, "term.labels")[attr(x, "order") == 1] %in% LHS_Y)))) stop("Outcomes cannot be direct predictors for covariates")
  # if (any(sapply(formsY,
  #                function(x) any(attr(x, "term.labels")[attr(x, "order") == 1] %in% LHS_Z)))) stop("Covariates cannot be direct predictors for outcomes")
  # if (any(sapply(formsY,  # could weaken this?
  #                function(x) any(attr(x, "term.labels")[attr(x, "order") == 1] %in% LHS_Y)))) stop("Outcomes cannot directly predict one another")


  if (ordering) {

    ## obtain empirical quantiles from plasmode variables that will appear in copula
    if (!is.null(dat)) {
      n <- nrow(dat)

      wh_q <- setdiff(unlist(rhs_vars(formulas[[2]])),
                      c(unlist(rhs_vars(formulas[[3]])), vars))
      quantiles <- dat[wh_q]

      for (i in seq_along(wh_q)) {
        var_nm <- wh_q[i]
        quan <- (rank(dat[[var_nm]])-1/2)/nrow(dat)
        if (any(duplicated(quan))) {
          cts <- table(quan)
          vals <- as.numeric(names(cts))
          for (j in which(cts > 1)) {
            wh_j <- which(quan == vals[j])
            quan[wh_j] <- quan[wh_j] + (runif(cts[j])-1/2)/n
          }
        }
        quantiles[wh_q] <- quan
      }
    }
    else wh_q <- character(0)

    nC <- length(wh_q)

    ## put some validation code in here for inversion method

    ## get formulas in right format
    if (!is.list(formulas[[4]])) {
      formulas[[4]] <- rep(list(formulas[[4]]), dY)
    }
    else if (length(formulas[[4]]) != dY) {
      formulas[[4]] <- rep(formulas[[4]][], dY)
    }
    if (!all(sapply(formulas[[4]], is.list)) || any(lengths(formulas[[4]]) < dZ+seq_len(dY)-1)) {
      for (i in seq_len(dY)) {
        if (is.list(formulas[[4]][[i]])) {
          formulas[[4]][[i]] <- rep(formulas[[4]][[i]], nC+dZ+i-1)
        }
        else {
          formulas[[4]][[i]] <- rep(list(formulas[[4]][[i]]), nC+dZ+i-1)
        }
        lhs(formulas[[4]][[i]]) <- c(wh_q, LHS_Z[rank(order[seq_len(dZ)])], LHS_Y[rank(order[dZ+dX+seq_len(i-1)])])
        ## update to allow some Zs to come after Ys
      }
    }
    names(formulas[[4]]) <- LHS_Y[rank(order[dZ+dX+seq_len(dY)])]

    ## get families in right format
    if (!is.list(family[[4]]) || length(family[[4]]) != dY) {
      if (is.list(family[[4]])) {
        family[[4]] <- rep(family[[4]], dY)
      }
      if (!is.list(family[[4]])) {
          if (length(family[[4]]) == dY) {
            family[[4]] <- as.list(family[[4]])
          }
          else  {
            family[[4]] <- rep(as.list(family[[4]]), dY)
          }
      }
    }
    if (any(lengths(family[[4]]) != lengths(formulas[[4]]))) {
      family[[4]] <- mapply(function(x, y) rep(x, y), family[[4]], nC+dZ+seq_len(dY)-1)
    }

    ## get parameters in right format
    if (!setequal(names(pars[[kwd]]), LHS_Y)) {
      if ("beta" %in% names(pars[[kwd]])) {
        pars[[kwd]] <- rep(list(list(list(beta=pars[[kwd]]$beta))), length(LHS_Y))
        names(pars[[kwd]]) <- LHS_Y
        if (nC+dZ > 1 || dY > 1) {
          for (i in seq_len(dY)) {
            vnm <- vars[order[dZ+dX+i]]
            pars[[kwd]][[vnm]] <- rep(pars[[kwd]][[vnm]], nC+dZ+i-1)
            names(pars[[kwd]][[vnm]]) <- c(wh_q, LHS_Z, LHS_Y[rank(order[dZ+dX+seq_len(i-1)])])
          }
        }
      }
      else {
        nrep <- setdiff(LHS_Y, names(pars[[kwd]]))
        if (length(nrep) > 0) stop(paste0("Variable ", paste(nrep, collapse=", "),
                                          "not represented in the copula parameters list"))
        rep <- setdiff(names(pars[[kwd]]), LHS_Y)
        if (length(rep) > 0) stop(paste0("Variable ", paste(rep, collapse=", "),
                                         "represented in copula parameters list but not a response variable"))
        stop("Shouldn't get here")
      }
    }

    # ## TRY TO CODE UP CHECK THAT EVERY COMBINATION HAS A beta PARAMETER
    # if (any(!sapply(lapply(pars[[kwd]], names), "beta" %in% x))) {
    #   stop("")
    # }

    ## get family in right format
    if (!is.matrix(famCop)) {
      if (dY == 1) {
        family[[4]] <- list(unlist(family[[4]]))
        names(family[[4]]) <- LHS_Y
      }
      else if (is.matrix(family[[4]]) && nrow(family[[4]] == dY)) {
        family[[4]] <- apply(family[[4]], 1, c, simplify = FALSE)
      }
      else if (is.list(family[[4]]) && length(family[[4]]) == dY) {
        if (any(lengths(family[[4]]) != dZ + seq_len(dY) - 1)) stop("Incorrect format of family parameters for copula with inversion method")
      }
    }
  }
  else {
    ## ensure that copula parameters are in correct matrix format
    if (!is.matrix(pars[[kwd]]$beta)) {
      pars[[kwd]]$beta <- matrix(pars[[kwd]]$beta, ncol=choose(dZ+dY,2))
      form_cop <- terms(formulas[[4]][[1]])
      if (nrow(pars[[kwd]]$beta) != length(attr(form_cop, "term.labels"))+attr(form_cop, "intercept")) stop("Wrong number of regression parameters for copula")
    }
    if (length(in_form) > 0) {
      if (any(pars[[kwd]]$beta[,in_form] != 0)) {
        warning("Parameters in copula for objects already related by regression; setting copula parameters to zero")
        pars[[kwd]]$beta[,in_form] <- 0
      }
    }

    ## code to check if Y or Z is included in copula formula
    if (any(output %in% rhs_vars(formulas[[4]])[[1]])) stop("copula cannot depend upon Z or Y variables")
  }

  if (!exists("quantiles")) quantiles <- NULL

  ## set up link functions
  link <- link_setup(link, family[1:3], vars=list(LHS_Z,LHS_X,LHS_Y))

    return(list(formulas=formulas, pars=pars, family=family, link=link,
              LHSs=list(LHS_Z=LHS_Z, LHS_X=LHS_X, LHS_Y=LHS_Y),
              quantiles=quantiles, kwd=kwd, dim=dims, vars=vars, output=output,
              order=order))
}

##' Process input for family variables
##'
##' @param family list of family parameters
##'
process_family <- function (family) {
  lens <- lengths(family)

  for (i in seq_along(lens)) {
    if (lens[i] > 0) {
      if (is.character(family[[i]][1])) {
        fams <- unlist(family[[i]])
        family[[i]] <- lapply(fams, function(x) get_family(x)())
      }
      else if (is.numeric(family[[i]][1])) {
        next
      }
      else if (all(sapply(family[[i]], function(x) is(x, "causl_family")))) {
        next
      }
    }
  }

  return(family)
}
