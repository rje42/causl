##' Process formulas, families and parameters
##'
##' @inheritParams rfrugalParam
##' @param kwd keyword for copula
## @param ordering logical: should an ordering of variables be computed?
## @param ... dots from rfrugalParam
##'
##' @details Function that processes and checks the validity of the main arguments
##' used for simulating data.
##'
process_inputs <- function (formulas, pars, family, link, dat, kwd, method="inversion") {

  ## process univariate formulas and obtain dimensions of model
  formulas <- process_formulas(formulas)
  ## check ordering is sound

  ord <- var_order(formulas, method=method)
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
  if (any(duplicated(vars))) stop("Repeated variable names not allowed")
  LHSs <- list(LHS_Z, LHS_X, LHS_Y)

  ## check univariate parameter values are appropriate
  if (missing(pars)) stop("Must supply parameter values")
  dummy_dat <- gen_dummy_dat(family=family, pars=pars, dat=dat, LHSs=LHSs, dims=dims)
  pars <- check_pars(formulas=formulas, family=family, pars=pars, dummy_dat=dummy_dat,
                     LHSs=LHSs, kwd=kwd, dims=dims)$pars

  ## different approaches based on method selected
  if (method == "rejection") pars <- check_rej(formulas=formulas, family=family, pars=pars, dims=dims, kwd=kwd)$pars
  else {

    ## obtain empirical quantiles from any plasmode variables that will appear in copula
    if (!is.null(dat)) {
      n <- nrow(dat)

      wh_q <- setdiff(unlist(rhs_vars(formulas[[2]])),
                      c(unlist(rhs_vars(formulas[[3]])), vars))
      quantiles <- process_prespecified(dat, prespec = wh_q)
    }
    else wh_q <- character(0)

    # nC <- length(wh_q)

    ## put some validation code in here for inversion method

    if (method == "inversion") {
      tmp <- pair_copula_setup(formulas=formulas[[4]], family=family[[4]],
                               pars=pars[[kwd]], LHSs=LHSs, quans=wh_q, ord=ord)
      formulas[[4]] <- tmp$formulas
      family[[4]] <- tmp$family
      pars[[kwd]] <- tmp$pars
    }
    else if (method == "multi_copula") {
      if (length(formulas[[4]]) != 1) stop("Only one formula for multivariate copula")
      if (length(family[[4]]) != 1) stop("Only one family for multivariate copula")
    }
    else stop("'method' should be \"inversion\", \"multi_copula\" or \"rejection\"")
  }

  if (!exists("quantiles")) quantiles <- NULL

  ## set up link functions
  link <- link_setup(link, family[1:3], vars=LHSs)

  return(list(formulas=formulas, pars=pars, family=family, link=link,
              LHSs=list(LHS_Z=LHS_Z, LHS_X=LHS_X, LHS_Y=LHS_Y),
              quantiles=quantiles, kwd=kwd, dim=dims[1:3], vars=vars, #output=output,
              order=ord))
}

##' @inheritParams process_inputs
##' @param len number of formulas
##' @describeIn process_inputs Process input for family variables
process_formulas <- function (formulas, len=4) {
  ## check we have four groups of formulas
  if (length(formulas) != len) stop(paste0("Formulas must have length ", len))

  ## ensure all formulas are given as lists
  if (any(sapply(formulas, class) == "formula")) {
    wh <- which(sapply(formulas, class) == "formula")
    for (i in wh) formulas[[i]] <- list(formulas[[i]])
  }

  return(formulas)
}

##' @inheritParams gen_dummy_dat
##' @describeIn process_inputs Process input for family variables
##' @param func_return function to use to process character arguments
##'
##' @details
##' For `causl` we use the `get_family()` function to process character based
##' arguments, but we allow for other functions to be used in packages that
##' build on this one.
##'
##'
##' @export
process_family <- function (family, dims, func_return=get_family) {

  if (missing(family)) {
    ## assume everything is Gaussian
    family <- lapply(dims, rep.int, x=1)
    return(family)
  }

  nU <- length(family) - 1

  if (is.list(family)) {
    lens <- lengths(family)

    if (any(sapply(family, class) == "causl_family")) {
      wh <- which(sapply(family, class) == "causl_family")
      for (i in wh) {
        family[[i]] <- list(family[[i]])
        lens[i] <- 1
      }
    }
    if (!all(lens[seq_len(nU)] == dims[seq_len(nU)])) stop("Mismatch in family and formulae specifications")
  }
  else if (length(family) == nU+1) {
    if (sum(dims[seq_len(nU)]) > nU) stop("Mismatch in family and formulae specification")
    family <- as.list(family)
    lens <- rep(1, nU+1)
  }
  else stop(paste0("'family' should be a list, or vector of length ", nU+1))

  ## deal with family names and causl_fam functions
  for (i in seq_len(nU)) {
    if (lens[i] > 0) {
      if (is.character(family[[i]])) {
        fams <- unlist(family[[i]])
        family[[i]] <- lapply(fams, function(x) func_return(x)())
        # family[[i]] <- lapply(fams, function(x) get_family(x)())
      }
      else if (is.numeric(family[[i]])) {
        next
      }
      else if (is.list(family[[i]])) {
        if (all(rapply(family[[i]], is.character))) {
          family[[i]] <- unlist(family[[i]])
          if (length(family[[i]]) != lens[i]) stop("Incorrect number of families specified")
          family[[i]] <- lapply(family[[i]], function(x) func_return(x)())
        }
        else if (all(rapply(family[[i]], is.numeric))) {
          family[[i]] <- unlist(family[[i]])
          if (length(family[[i]]) != lens[i]) stop("Incorrect number of families specified")
        }
        else if (all(sapply(family[[i]], function(x) is(x, "causl_family")))) {
          next
        }
        else stop("Not a valid family specification")
      }
      else stop("Not a valid family specification")
    }
  }

  ## check copula families are valid
  # if (!all(unlist(family[1:3]) %in% family_vals$val)) stop("Invalid family specification")
  if (!all(unlist(family[[nU+1]]) %in% copula_vals$val)) stop("Invalid copula specification")

  return(family)
}

##' Obtain variable ordering from formulas
##'
##' @inheritParams gen_dummy_dat
##' @inheritParams process_inputs
##' @param inc_cop logical indicating whether to include copula in the ordering
##'
var_order <- function (formulas, dims, inc_cop=TRUE, method) {

  if (missing(dims)) dims <- lengths(formulas)
  nU <- length(dims) - 1

  ## get terms objects
  forms <- lapply(formulas[seq_len(nU)], function(x) lapply(x, terms))

  ord_mat <- matrix(0, nrow=sum(dims[seq_len(nU)]), ncol=sum(dims[seq_len(nU)]))

  LHSs <- lapply(formulas[seq_len(nU)], lhs)
  vars <- unlist(LHSs)

  ## get adjacency matrix of dependencies
  for (i in seq_along(forms)) if (dims[i] > 0) {
    cvs <- sum(dims[seq_len(i-1)]) + seq_len(dims[i])  # variables for this iteration
    ord_mat[cvs, ] <- 1*t(sapply(forms[[i]], function(x) vars %in% attr(x, "term.labels")))
  }

  if (inc_cop) {
    ## now include copula parameters
    LHS_Y <- LHSs[[3]]

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

    ## add this information in
    ord_mat[sum(dims[1:2]) + seq_len(dims[3]), ] <- ord_mat[sum(dims[1:2]) + seq_len(dims[3]), ] +
      1*t(sapply(forms[[3]], function(x) vars %in% attr(x, "term.labels")))
  }

  ## get order for simulation, if one exists

  if (FALSE && method == "multi_copula") {
  # break ties with Xs going first.
  # if there are uncoditional treatments put them first
  treats <- ord_mat[(dims[1] + 1) : (dims[1] + dims[2]), , drop = FALSE]
  nuisance <- ord_mat[1:dims[1], , drop = FALSE]
  marginal_treats <- apply(treats == 0, 1, all)
  if(any(marginal_treats)){
    which_treats <- which(marginal_treats) + dims[1]
    marginal_nuisance <- apply(nuisance == 0, 1, all)
    if(any(marginal_nuisance)){
      which_nuisance <- which(marginal_nuisance)
      ord_mat[which_nuisance, which_treats] <- 0.1
    }
  }
}
  ord <- topOrd(ord_mat)

  if (any(is.na(ord))) stop("Formulae contain cyclic dependencies")

  return(ord)
}

##' @describeIn check_pars Checks for rejection sampling
##' @inheritParams process_inputs
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

##' Generate a dummy dataset
##'
##' Create a dummy dataset for the purpose of checking coefficient numbers
##'
##' @inheritParams process_inputs
##' @inheritParams rfrugalParam
##' @param LHSs left-hand sides from `formulas`
##' @param dims number of variables in each class
##'
##' @export
gen_dummy_dat <- function (family, pars, dat, LHSs, dims) {

  if (missing(dims)) dims <- lengths(family)
  nU <- length(LHSs)

  ## produce dummy data.frame to check number of coefficients
  nv <- sum(dims[seq_len(nU)])
  if (!is.null(attr(LHSs, "T"))) {
    T <- attr(LHSs, "T")[1]
    strt <- attr(LHSs, "T")[2]
    if (is.na(strt)) strt <- 0
    nv <- nv + (T-1)*sum(dims[2:4])
    nms <- LHSs[[1]]
    nms_t <- unlist(LHSs[2:4])
    nms <- c(nms, paste0(nms_t, "_", rep(seq_len(T)+strt-1, each=length(nms_t))))
  }
  else {
    nms <- unlist(LHSs)
    nv <- length(nms)
  }
  out <- as.data.frame(rep(list(NA), nv))
  names(out) <- nms

  if (!is.null(dat)) {
    out <- cbind(dat[1,,drop=FALSE], out)
  }

  ## generate dummy data for each variable to simulate
  for (i in seq_along(LHSs)) for (j in seq_len(dims[[i]])) {
    if (i > 1 && exists("strt")) {
      vnms <- paste0(LHSs[[i]][j], "_", seq_len(T)-1-strt)
      if (!is_categorical(family[[i]][j])) out[vnms] <- 0
      else out[vnms] <- factor(x=1L, levels=seq_len(pars[[LHSs[[i]][j]]]$nlevel))
    }
    else {
      if (!is_categorical(family[[i]][j])) out[[LHSs[[i]][j]]] <- 0
      else out[[LHSs[[i]][j]]] <- factor(x=1L, levels=seq_len(pars[[LHSs[[i]][j]]]$nlevel))
    }
  }

  return(out)
}

##' Check parameters for univariate families
##'
##' Checks existence of beta vectors and then assesses appropriate length
##'
##' @inheritParams process_inputs
##' @inheritParams gen_dummy_dat
##' @param dummy_dat a dummy dataset, as generated by `gen_dummy_dat()`
##'
##' @export
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
    stop(paste(ifelse(plur, "Variables", "Variable"),
               paste(names(pars)[!bpres], collapse=", "),
               ifelse(plur, "lack", "lacks"), "a beta parameter vector"))
  }

  ## check variable names
  vars <- unlist(LHSs)
  rep_pars <- match(vars, nm_pars, nomatch = 0L)
  if (any(rep_pars == 0)) {
    wh_nrep <- vars[rep_pars==0]
    stop(paste0("Variable ", paste(wh_nrep, collapse=", "), " not represented in the parameters list"))
  }


  for (j in seq_along(LHSs)) for (i in seq_along(formulas[[j]])) {
    mod_mat <- model.matrix(formulas[[j]][[i]], data=dummy_dat)
    # npar <- length(attr(formsZ[[i]], "term.labels")) + attr(formsZ[[i]], "intercept")
    npar <- ncol(mod_mat)
    if (is_categorical(family[[j]][i])) {
      npar <- npar*(pars[[LHSs[[j]][i]]]$nlevel - 1)
      if (!is.matrix(pars[[LHSs[[j]][i]]]$beta)) pars[[LHSs[[j]][i]]]$beta <- matrix(pars[[LHSs[[j]][i]]]$beta, nrow=ncol(mod_mat))
    }
    if (isTRUE(is.na(npar)) || length(npar) != 1) stop(paste0("Categorical variable '", LHSs[[j]][i], "' requires an 'nlevel' parameter"))
    if (length(pars[[LHSs[[j]][i]]]$beta) != npar) stop(paste0("dimension of model matrix for ", LHSs[[j]][i], " does not match number of coefficients provided"))
  }
  # for (i in seq_along(formulas[[2]])) {
  #   mod_mat <- model.matrix(formulas[[2]][[i]], data=dummy_dat)
  #   # npar <- length(attr(formsX[[i]], "term.labels")) + attr(formsX[[i]], "intercept")
  #   npar <- ncol(mod_mat)
  #   if (is_categorical(family[[2]][i])) {
  #     npar <- npar*(pars[[LHSs[[2]][i]]]$nlevel - 1)
  #     if (!is.matrix(pars[[LHSs[[2]][i]]]$beta)) pars[[LHSs[[2]][i]]]$beta <- matrix(pars[[LHSs[[2]][i]]]$beta, nrow=ncol(mod_mat))
  #   }
  #   if (isTRUE(is.na(npar)) || length(npar) != 1) stop(paste0("Categorical variable '", LHSs[[2]][i], "' requires an 'nlevel' parameter"))
  #   if (length(pars[[LHSs[[2]][i]]]$beta) != npar) stop(paste0("dimension of model matrix for '", LHSs[[2]][i], "' does not match number of coefficients provided"))
  # }
  # for (i in seq_along(formulas[[3]])) {
  #   mod_mat <- model.matrix(formulas[[3]][[i]], data=dummy_dat)
  #   # npar <- length(attr(formsY[[i]], "term.labels")) + attr(formsY[[i]], "intercept")
  #   npar <- ncol(mod_mat)
  #   if (is_categorical(family[[3]][i])) {
  #     npar <- npar*(pars[[LHSs[[3]][i]]]$nlevel - 1)
  #     if (!is.matrix(pars[[LHSs[[3]][i]]]$beta)) pars[[LHSs[[3]][i]]]$beta <- matrix(pars[[LHSs[[3]][i]]]$beta, nrow=ncol(mod_mat))
  #   }
  #   if (isTRUE(is.na(npar)) || length(npar) != 1) stop(paste0("Categorical variable '", LHSs[[3]][i], "' requires an 'nlevel' parameter"))
  #   if (length(pars[[LHSs[[3]][i]]]$beta) != npar) stop(paste0("dimension of model matrix for ", LHSs[[3]][i], " does not match number of coefficients provided"))
  # }


  ## check that families have necessary parameters
  for (j in seq_along(LHSs)) {
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
  cop_pos <- length(LHSs) + 1   # should be 4 for causl
  for (i in seq_along(formulas[[cop_pos]])) {
    pars_cop <- pars[[kwd]]

    if (is(formulas[[cop_pos]][[i]], "formula")) {
      mod_mat <- model.matrix(formulas[[cop_pos]][[i]], data=dummy_dat)
      npar <- ncol(mod_mat)
      # if (length(pars_cop[[LHSs[[3]][i]]]$beta) != npar) stop(paste0("dimension of model matrix for ", LHSs[[3]][i], " does not match number of coefficients provided"))
    }
    else if (is.list(formulas[[cop_pos]][[i]])) {
      for (j in seq_along(formulas[[cop_pos]][[i]])) {
        mod_mat <- model.matrix(formulas[[cop_pos]][[i]][[j]], data=dummy_dat)
        npar <- ncol(mod_mat)
      }
    }

    # npar <- length(attr(formsY[[i]], "term.labels")) + attr(formsY[[i]], "intercept")
  }


  return(list(pars=pars))
}

##' Sets up copula quantities only
##'
## @inheritParams process_inputs
##' @param formulas list of formulas for copula only
##' @param family list of families for copula only
##' @param pars list of copula parameters
##' @param LHSs left-hand sides for all variables
##' @param quans character vector of already existing variables to include
##' @param ord topological ordering
##'
##' @export
pair_copula_setup <- function (formulas, family, pars, LHSs, quans, ord) {

  dims <- lengths(LHSs)
  dZ <- dims[length(dims)-2]
  dX <- dims[length(dims)-1]
  dY <- last(dims)
  LHS_Z <- LHSs[[length(LHSs)-2]]
  LHS_Y <- LHSs[[length(LHSs)]]

  nQ <- length(quans)

  if (is.list(formulas)) {
    ## put in code to check that copula formulas and parameters have matching lengths

  }

  ## get formulas in right format
  if (!is.list(formulas)) {
    formulas <- rep(list(formulas), dY)
  }
  else if (length(formulas) != dY) {
    formulas <- rep(formulas[], dY)
  }
  if (!all(sapply(formulas, is.list)) || any(lengths(formulas) < dZ+seq_len(dY)-1)) {
    for (i in seq_len(dY)) {
      if (is.list(formulas[[i]])) {
        formulas[[i]] <- rep(formulas[[i]], nQ+dZ+i-1)
      }
      else {
        formulas[[i]] <- rep(list(formulas[[i]]), nQ+dZ+i-1)
      }
      lhs(formulas[[i]]) <- c(quans, LHS_Z[rank(ord[seq_len(dZ)])],
                              LHS_Y[rank(ord[dZ+dX+seq_len(i-1)])])
      ## update to allow some Zs to come after Ys
    }
  }
  names(formulas) <- LHS_Y[rank(ord[dZ+dX+seq_len(dY)])]


  ## get copula families in right format
  if (!is.list(family) || length(family) != dY) {
    if (is.list(family)) {
      family <- rep(family, dY)
    }
    if (!is.list(family)) {
      if (length(family) == dY) {
        family <- as.list(family)
      }
      else  {
        family <- rep(as.list(family), dY)
      }
    }
  }
  if (any(lengths(family) != lengths(formulas))) {
    family <- mapply(function(x, y) rep(x, y), family, nQ+dZ+seq_len(dY)-1, SIMPLIFY = FALSE)
  }
  # if (length(dims) == 4) return(list(formulas=formulas, family=family))


  ## get parameters in right format
  if (!setequal(names(pars), LHS_Y)) {
    if ("beta" %in% names(pars)) {
      pars <- rep(list(list(list(beta=pars$beta))), length(LHS_Y))
      names(pars) <- LHS_Y
      ordY <- rank(ord[dZ+dX+seq_len(dY)])
      if (nQ+dZ > 1 || dY > 1) {
        for (i in seq_len(dY)) {
          vnm <- LHS_Y[[ordY[i]]]
          pars[[vnm]] <- rep(pars[[vnm]], nQ+dZ+i-1)
          names(pars[[vnm]]) <- c(quans, LHS_Z, LHS_Y[ord[seq_len(i-1)]])
        }
      }
    }
    else {
      nrep <- setdiff(LHS_Y, names(pars))
      if (length(nrep) > 0) stop(paste0("Variable ", paste(nrep, collapse=", "),
                                        " not represented in the copula parameters list"))
      rep <- setdiff(names(pars), LHS_Y)
      if (length(rep) > 0) stop(paste0("Variable ", paste(rep, collapse=", "),
                                       " represented in copula parameters list but not a response variable"))
      stop("Shouldn't get here")
    }
  }
  else if (!all(sapply(pars, is.list))) {
    ## put some code in here to correct fact that some pairs aren't represented
  }
  # else stop("'method' should be \"inversion\" or \"rejection\"")


  return(list(formulas=formulas, family=family, pars=pars))
}

##' Obtain quantiles for prespecified variables
##'
##' @param dat data frame containing variables
##' @param prespec character vector of prespecified variables in `dat`
##'
##' @details Currently takes the rank of each entry, and subtracts 1/2.  If
##' there are \eqn{k} ties they are randomly sorted with a uniform random variable
##' in the symmetric interval around the rank of width \eqn{k/n}.
##'
##' @export
process_prespecified <- function (dat, prespec) {
  quantiles <- dat[prespec]
  n <- nrow(dat)

  ## go through each variable and put into [0,1] in approx uniform manner
  for (i in seq_along(prespec)) {
    var_nm <- prespec[i]
    quan <- (rank(dat[[var_nm]])-1/2)/n
    if (any(duplicated(quan))) {
      cts <- table(quan)
      vals <- as.numeric(names(cts))
      for (j in which(cts > 1)) {
        wh_j <- which(quan == vals[j])
        quan[wh_j] <- quan[wh_j] + (runif(cts[j])-1/2)*cts[j]/n
      }
    }
    quantiles[prespec] <- quan
  }

  return(quantiles)
}
