##' Process formulas, families and parameters
##'
##' @inheritParams rfrugalParam
##' @param kwd keyword for copula
## @param ordering logical: should an ordering of variables be computed?
## @param ... dots from rfrugalParam
##'
process_inputs <- function (formulas, pars, family, link, dat, kwd, method="inversion") {

  ## process univariate formulas and obtain dimensions of model
  formulas <- process_formulas(formulas)
  ## check ordering is sound
  ord <- variable_order(formulas)
  dims <- lengths(formulas)

  ## obtain variable names
  LHS_Z <- lhs(formulas[[1]])
  LHS_X <- lhs(formulas[[2]])
  LHS_Y <- lhs(formulas[[3]])

  ## process family variable inputs
  family <- process_family(family=family, dims=dims)

  ### change to use the names supplied by each family

  ## useful variable summaries
  # output <- c(LHS_Z, LHS_Y)
  vars <- c(LHS_Z, LHS_X, LHS_Y)
  LHSs <- list(LHS_Z, LHS_X, LHS_Y)

  ## check univariate parameter values are appropriate
  if (missing(pars)) stop("Must supply parameter values")
  dummy_dat <- gen_dummy_dat(family=family, pars=pars, dat=dat, LHSs=LHSs, dims=dims)
  pars <- check_pars(formulas=formulas, family=family, pars=pars, dummy_dat=dummy_dat,
                     LHSs=LHSs, kwd=kwd, dims=dims)$pars

  ## different approaches based on method selected
  if (method == "rejection") pars <- check_rej(formulas=formulas, family=family, pars=pars, dims=dims, kwd=kwd)$pars
  else if (method == "inversion") {

    ## obtain empirical quantiles from any plasmode variables that will appear in copula
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

    # nC <- length(wh_q)

    ## put some validation code in here for inversion method

    tmp <- pair_copula_setup(formulas=formulas[[4]], family=family[[4]],
                             pars=pars[[kwd]], LHSs=LHSs, quans=wh_q, ord=ord)
    formulas[[4]] <- tmp$formulas
    family[[4]] <- tmp$family
    pars[[kwd]] <- tmp$pars
  }
  else stop("'method' should be \"inversion\" or \"rejection\"")

  if (!exists("quantiles")) quantiles <- NULL

  ## set up link functions
  link <- link_setup(link, family[1:3], vars=LHSs)

  return(list(formulas=formulas, pars=pars, family=family, link=link,
              LHSs=list(LHS_Z=LHS_Z, LHS_X=LHS_X, LHS_Y=LHS_Y),
              quantiles=quantiles, kwd=kwd, dim=dims[1:3], vars=vars, #output=output,
              order=ord))
}

process_formulas <- function (formulas) {
  ## check we have four groups of formulas
  if (length(formulas) != 4) stop("Formulas must have length 4")

  ## ensure all formulas are given as lists
  if (any(sapply(formulas, class) == "formula")) {
    wh <- which(sapply(formulas, class) == "formula")
    for (i in wh) formulas[[i]] <- list(formulas[[i]])
  }

  return(formulas)
}

##' Process input for family variables
##'
##' @param family list of family parameters
##'
##'
process_family <- function (family, dims) {

  if (missing(family)) {
    ## assume everything is Gaussian
    family <- lapply(dims, rep.int, x=1)
    return(family)
  }
  else if (is.list(family)) {
    lens <- lengths(family)

    if (any(sapply(family, class) == "causl_family")) {
      wh <- which(sapply(family, class) == "causl_family")
      for (i in wh) family[[i]] <- list(family[[i]])
    }
    if (!all(lens[1:3] == dims[1:3])) stop("Mismatch in family and formulae specifications")
  }
  else if (length(family) == 4) {
    if (sum(dims[1:3]) > 3) stop("Mismatch in family and formulae specification")
    family <- as.list(family)
    lens <- c(1,1,1,1)
  }
  else stop("family should be a list, or vector of length 4")

  ## deal with family names and causl_fam functions
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

  ## check copula families are valid
  # if (!all(unlist(family[1:3]) %in% familyVals$val)) stop("Invalid family specification")
  if (!all(unlist(family[[4]]) %in% copula_vals$val)) stop("Invalid copula specification")

  return(family)
}

variable_order <- function (formulas, dims) {

  if (missing(dims)) dims <- lengths(formulas)

  ## get terms objects
  formsZ <- lapply(formulas[[1]], terms)
  formsX <- lapply(formulas[[2]], terms)
  formsY <- lapply(formulas[[3]], terms)

  ord_mat <- matrix(0, nrow=sum(dims[1:3]), ncol=sum(dims[1:3]))

  LHS_Z <- lhs(formulas[[1]])
  LHS_X <- lhs(formulas[[2]])
  LHS_Y <- lhs(formulas[[3]])
  vars <- c(LHS_Z, LHS_X, LHS_Y)

  ## get adjacency matrix of dependencies
  if (dims[1] > 0) ord_mat[seq_len(dims[1]), ] <- 1*t(sapply(formsZ, function(x) vars %in% attr(x, "term.labels")))
  if (dims[2] > 0) ord_mat[dims[1] + seq_len(dims[2]), ] <- 1*t(sapply(formsX, function(x) vars %in% attr(x, "term.labels")))
  if (dims[3] > 0) ord_mat[sum(dims[1:2]) + seq_len(dims[3]), ] <- 1*t(sapply(formsY, function(x) vars %in% attr(x, "term.labels")))


  if (!is.null(formulas[[4]][LHS_Y])) {
    ## look for possible loops in copula formulae
    cop_frms <- formulas[[4]][LHS_Y]
    if (is.list(cop_frms) && length(cop_frms) > 0 && is.list(cop_frms[[1]])) {
      pred <- lapply(cop_frms, function (x) sapply(x, lhs))
      ord_mat[sum(dims[1:2]) + seq_len(dims[3]), ] <- ord_mat[sum(dims[1:2]) + seq_len(dims[3]), ] +
        1*t(sapply(pred, function(x) vars %in% x))
      ord_mat[] <- pmin(ord_mat, 1)
    }
  }

  ord_mat[sum(dims[1:2]) + seq_len(dims[3]), ] <- ord_mat[sum(dims[1:2]) + seq_len(dims[3]), ] +
    1*t(sapply(formsY, function(x) vars %in% attr(x, "term.labels")))

  ## get order for simulation, if one exists
  ord <- topOrd(ord_mat)
  if (any(is.na(ord))) stop("Formulae contain cyclic dependencies")

  return(ord)
}

##' Checks for rejection sampling
##'
check_rej <- function(formulas, family, pars, dims, kwd) {
  if (any(unlist(lapply(family, is_categorical)))) stop("Categorical variables must be simulated using 'method=\"inversion\"'")

  LHS_Z <- lhs(formulas[[1]])
  LHS_X <- lhs(formulas[[2]])
  LHS_Y <- lhs(formulas[[3]])
  if (missing(dims)) dims <- lengths(formulas)

  ## more restrictive for the rejection sampling method
  if (any(unlist(lapply(rhs_vars(formulas[[3]]),
                        function(x) any(x %in% LHS_Z))))) stop("Covariates cannot be direct predictors for outcomes under rejection sampling")
  if (any(unlist(lapply(rhs_vars(formulas[[1]]),
                        function(x) any(x %in% LHS_Y))))) stop("Outcomes cannot be direct predictors for covariates under rejection sampling")
  if (any(unlist(lapply(rhs_vars(formulas[[3]]),  # could relax this
                        function(x) any(x %in% LHS_Y))))) stop("Outcomes cannot directly predict one another under rejection sampling")

  ## check if covariates predict each other
  if (dims[1] > 1) {
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

  ## ensure that copula parameters are in correct matrix format
  if (!is.matrix(pars[[kwd]]$beta)) {
    pars[[kwd]]$beta <- matrix(pars[[kwd]]$beta, ncol=choose(dims[1]+dims[3],2))
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
  if (any(c(LHS_Z, LHS_Y) %in% unlist(rhs_vars(formulas[[4]])))) stop("copula cannot depend upon Z or Y variables")

  return(list(pars=pars))
}

gen_dummy_dat <- function (family, pars, dat, LHSs, dims) {

  if (missing(dims)) dims <- lengths(family)

  ## produce dummy data.frame to check number of coefficients
  dummy_dat <- as.data.frame(rep(list(NA), sum(dims[1:3])))
  names(dummy_dat) <- unlist(LHSs)
  if (!is.null(dat)) {
    dummy_dat <- cbind(dat[1,], dummy_dat)
  }

  ## generate dummy data for each variable to simulate
  for (i in 1:3) for (j in seq_len(dims[[i]])) {
    if (!is_categorical(family[[i]][j])) dummy_dat[[LHSs[[i]][j]]] <- 0
    else dummy_dat[[LHSs[[i]][j]]] <- factor(x=1L, levels=seq_len(pars[[LHSs[[i]][j]]]$nlevel))
  }

  return(dummy_dat)
}

## Function checks existence of beta vectors and then assesses appropriate length
check_pars <- function (formulas, family, pars, dummy_dat, LHSs, kwd, dims) {

  if (missing(dims)) dims <- lengths(formulas)

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

  ## check variable names
  vars <- unlist(LHSs)
  rep_pars <- match(vars, nm_pars, nomatch = 0L)
  if (any(rep_pars == 0)) {
    wh_nrep <- vars[rep_pars==0]
    stop(paste0("Variable ", paste(wh_nrep, collapse=", "), "not represented in the parameters list"))
  }


  for (i in seq_along(formulas[[1]])) {
    mod_mat <- model.matrix(formulas[[1]][[i]], data=dummy_dat)
    # npar <- length(attr(formsZ[[i]], "term.labels")) + attr(formsZ[[i]], "intercept")
    npar <- ncol(mod_mat)
    if (is_categorical(family[[1]][i])) {
      npar <- npar*(pars[[LHSs[[1]][i]]]$nlevel - 1)
      if (!is.matrix(pars[[LHSs[[1]][i]]]$beta)) pars[[LHSs[[1]][i]]]$beta <- matrix(pars[[LHSs[[1]][i]]]$beta, nrow=ncol(mod_mat))
    }
    if (isTRUE(is.na(npar)) || length(npar) != 1) stop(paste0("Categorical variable '", LHSs[[1]][i], "' requires an 'nlevel' parameter"))
    if (length(pars[[LHSs[[1]][i]]]$beta) != npar) stop(paste0("dimension of model matrix for ", LHSs[[1]][i], " does not match number of coefficients provided"))
  }
  for (i in seq_along(formulas[[2]])) {
    mod_mat <- model.matrix(formulas[[2]][[i]], data=dummy_dat)
    # npar <- length(attr(formsX[[i]], "term.labels")) + attr(formsX[[i]], "intercept")
    npar <- ncol(mod_mat)
    if (is_categorical(family[[2]][i])) {
      npar <- npar*(pars[[LHSs[[2]][i]]]$nlevel - 1)
      if (!is.matrix(pars[[LHSs[[2]][i]]]$beta)) pars[[LHSs[[2]][i]]]$beta <- matrix(pars[[LHSs[[2]][i]]]$beta, nrow=ncol(mod_mat))
    }
    if (isTRUE(is.na(npar)) || length(npar) != 1) stop(paste0("Categorical variable '", LHSs[[2]][i], "' requires an 'nlevel' parameter"))
    if (length(pars[[LHSs[[2]][i]]]$beta) != npar) stop(paste0("dimension of model matrix for '", LHSs[[2]][i], "' does not match number of coefficients provided"))
  }
  for (i in seq_along(formulas[[3]])) {
    mod_mat <- model.matrix(formulas[[3]][[i]], data=dummy_dat)
    # npar <- length(attr(formsY[[i]], "term.labels")) + attr(formsY[[i]], "intercept")
    npar <- ncol(mod_mat)
    if (is_categorical(family[[3]][i])) {
      npar <- npar*(pars[[LHSs[[3]][i]]]$nlevel - 1)
      if (!is.matrix(pars[[LHSs[[3]][i]]]$beta)) pars[[LHSs[[3]][i]]]$beta <- matrix(pars[[LHSs[[3]][i]]]$beta, nrow=ncol(mod_mat))
    }
    if (isTRUE(is.na(npar)) || length(npar) != 1) stop(paste0("Categorical variable '", LHSs[[3]][i], "' requires an 'nlevel' parameter"))
    if (length(pars[[LHSs[[3]][i]]]$beta) != npar) stop(paste0("dimension of model matrix for ", LHSs[[3]][i], " does not match number of coefficients provided"))
  }


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

  ## check copula regression parameters are the correct length
  for (i in seq_along(formulas[[4]])) {
    pars_cop <- pars[[kwd]]

    if (is(formulas[[4]][[i]], "formula")) {
      mod_mat <- model.matrix(formulas[[4]][[i]], data=dummy_dat)
      npar <- ncol(mod_mat)
      # if (length(pars_cop[[LHSs[[3]][i]]]$beta) != npar) stop(paste0("dimension of model matrix for ", LHSs[[3]][i], " does not match number of coefficients provided"))
    }
    else if (is.list(formulas[[4]][[i]])) {
      for (j in seq_along(formulas[[4]][[i]])) {
        mod_mat <- model.matrix(formulas[[4]][[i]][[j]], data=dummy_dat)
        npar <- ncol(mod_mat)
      }
    }

    # npar <- length(attr(formsY[[i]], "term.labels")) + attr(formsY[[i]], "intercept")
  }


  return(list(pars=pars))
}

##' Sets up copula quantities only
pair_copula_setup <- function (formulas, family, pars, LHSs, quans, ord) {

  dims <- lengths(LHSs)
  nQ <- length(quans)

  if (is.list(formulas)) {
    ## put in code to check that copula formulas and parameters have matching lengths

  }

  ## get formulas in right format
  if (!is.list(formulas)) {
    formulas <- rep(list(formulas), dims[3])
  }
  else if (length(formulas) != dims[3]) {
    formulas <- rep(formulas[], dims[3])
  }
  if (!all(sapply(formulas, is.list)) || any(lengths(formulas) < dims[1]+seq_len(dims[3])-1)) {
    for (i in seq_len(dims[3])) {
      if (is.list(formulas[[i]])) {
        formulas[[i]] <- rep(formulas[[i]], nC+dims[1]+i-1)
      }
      else {
        formulas[[i]] <- rep(list(formulas[[i]]), nQ+dims[1]+i-1)
      }
      lhs(formulas[[i]]) <- c(quans, LHSs[[1]][rank(ord[seq_len(dims[1])])],
                              LHSs[[3]][rank(ord[dims[1]+dims[2]+seq_len(i-1)])])
      ## update to allow some Zs to come after Ys
    }
  }
  names(formulas) <- LHSs[[3]][rank(ord[dims[1]+dims[2]+seq_len(dims[3])])]


  ## get parameters in right format
  if (!setequal(names(pars), LHSs[[3]])) {
    if ("beta" %in% names(pars)) {
      pars <- rep(list(list(list(beta=pars$beta))), length(LHSs[[3]]))
      names(pars) <- LHSs[[3]]
      ordY <- rank(ord[dims[1]+dims[2]+seq_len(dims[3])])
      if (nQ+dims[1] > 1 || dims[3] > 1) {
        for (i in seq_len(dims[3])) {
          vnm <- LHSs[[3]][[ordY[i]]]
          pars[[vnm]] <- rep(pars[[vnm]], nQ+dims[1]+i-1)
          names(pars[[vnm]]) <- c(quans, LHSs[[1]], LHSs[[3]][ord[seq_len(i-1)]])
        }
      }
    }
    else {
      nrep <- setdiff(LHSs[[3]], names(pars))
      if (length(nrep) > 0) stop(paste0("Variable ", paste(nrep, collapse=", "),
                                        "not represented in the copula parameters list"))
      rep <- setdiff(names(pars), LHSs[[3]])
      if (length(rep) > 0) stop(paste0("Variable ", paste(rep, collapse=", "),
                                       "represented in copula parameters list but not a response variable"))
      stop("Shouldn't get here")
    }
  }

  # ## TRY TO CODE UP CHECK THAT EVERY COMBINATION HAS A beta PARAMETER
  # if (any(!sapply(lapply(pars, names), "beta" %in% x))) {
  #   stop("")
  # }

  ## get copula families in right format
  if (!is.list(family) || length(family) != dims[3]) {
    if (is.list(family)) {
      family <- rep(family, dims[3])
    }
    if (!is.list(family)) {
      if (length(family) == dims[3]) {
        family <- as.list(family)
      }
      else  {
        family <- rep(as.list(family), dims[3])
      }
    }
  }
  if (any(lengths(family) != lengths(formulas))) {
    family <- mapply(function(x, y) rep(x, y), family, nQ+dims[1]+seq_len(dims[3])-1, SIMPLIFY = FALSE)
  }
  #
  # ## get copula families in right format
  # if (!is.matrix(famCop)) {
  #   if (dY == 1) {
  #     family <- list(unlist(family))
  #     names(family) <- LHSs[[3]]
  #   }
  #   else if (is.matrix(family) && nrow(family == dY)) {
  #     family <- apply(family, 1, c, simplify = FALSE)
  #   }
  #   else if (is.list(family) && length(family) == dY) {
  #     if (any(lengths(family) != dims[1] + seq_len(dY) - 1)) stop("Incorrect format of family parameters for copula with inversion method")
  #   }
  # }

  return(list(formulas=formulas, family=family, pars=pars))
}
