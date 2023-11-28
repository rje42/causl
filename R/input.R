##' Process formulas, families and parameters
##'
##' @param formulas list of lists of formulas
##' @param pars list of lists of parameters
##' @param family families for Z,X,Y and copula
##' @param link list of link functions
##' @param kwd keyword for copula
##' @param ordering logical: should an ordering of variables be computed?
##' @param ... dots from rfrugalParam
##'
process_inputs <- function (formulas, pars, family, link, kwd, ordering=FALSE, ...) {

  # ## check list of formula
  # if (ordering && !is.numeric(family[[1]]) && ("formula" %in% class(formulas[[1]]))) {
  #   process_inputs2(formulas, pars, family, link, kwd)
  # }

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
    if (!all(lengths(family[1:3]) == lengths(formulas[1:3]))) stop("Mismatch in family and formulae specifications")
  }
  else if (length(family) == 4) {
    if (sum(lengths(formulas[1:3])) > 3) stop("Mismatch in family and formulae specification")
    family <- as.list(family)
  }
  else stop("family should be a list, or vector of length 4")

  add_args <- list(...)
  if (length(add_args) > 0) {
    nms <- names(add_args)

  }

  ## check families are valid
  # if (!all(unlist(family[1:3]) %in% familyVals$val)) stop("Invalid family specification")
  if (!all(unlist(family[[4]]) %in% copulaVals$val)) stop("Invalid copula specification")

  ## check that supplied parameters are sufficient
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

  for (i in seq_along(formulas[[1]])) {
    npar <- length(attr(formsZ[[i]], "term.labels")) + attr(formsZ[[i]], "intercept")
    if (length(pars[[LHS_Z[i]]]$beta) != npar) stop(paste0("dimension of model matrix for ", LHS_Z[i], " does not match number of coefficients provided"))
  }
  for (i in seq_along(formulas[[2]])) {
    npar <- length(attr(formsX[[i]], "term.labels")) + attr(formsX[[i]], "intercept")
    if (length(pars[[LHS_X[i]]]$beta) != npar) stop(paste0("dimension of model matrix for ", LHS_X[i], " does not match number of coefficients provided"))
  }
  for (i in seq_along(formulas[[3]])) {
    npar <- length(attr(formsY[[i]], "term.labels")) + attr(formsY[[i]], "intercept")
    if (length(pars[[LHS_Y[i]]]$beta) != npar) stop(paste0("dimension of model matrix for ", LHS_Y[i], " does not match number of coefficients provided"))
  }

  ## check which variables are already in the formula for another, and set
  ## corresponding copula correlation to zero
  # formsZ <- merge_formulas(formulas[[1]])$old_forms
  # setmatch(LHS_Z, rhs_vars(formsZ))
  # formsY <- merge_formulas(formulas[[3]])$old_forms
  if (ordering) {
    ord_mat <- matrix(0, nrow=sum(dims), ncol=sum(dims))

    ## get adjacency matrix of dependencies
    vars <- c(LHS_Z, LHS_X, LHS_Y)
    ord_mat[seq_len(dZ), ] <- 1*t(sapply(formsZ, function(x) vars %in% attr(x, "term.labels")))
    ord_mat[dZ + seq_len(dX), ] <- 1*t(sapply(formsX, function(x) vars %in% attr(x, "term.labels")))
    ord_mat[dZ + dX + seq_len(dY), ] <- 1*t(sapply(formsY, function(x) vars %in% attr(x, "term.labels")))

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
    if (any(unlist(lapply(rhs_vars(formulas[[3]]),
                          function(x) any(x %in% LHS_Z))))) stop("Covariates cannot be direct predictors for outcomes")
    if (any(unlist(lapply(rhs_vars(formulas[[1]]),
                          function(x) any(x %in% LHS_Y))))) stop("Outcomes cannot be direct predictors for covariates")
    if (any(unlist(lapply(rhs_vars(formulas[[3]]),  # could relax this
                          function(x) any(x %in% LHS_Y))))) stop("Outcomes cannot directly predict one another")
    order <- NULL
  }
  #
  # if (any(sapply(formsZ,
  #        function(x) any(attr(x, "term.labels")[attr(x, "order") == 1] %in% LHS_Y)))) stop("Outcomes cannot be direct predictors for covariates")
  # if (any(sapply(formsY,
  #                function(x) any(attr(x, "term.labels")[attr(x, "order") == 1] %in% LHS_Z)))) stop("Covariates cannot be direct predictors for outcomes")
  # if (any(sapply(formsY,  # could weaken this?
  #                function(x) any(attr(x, "term.labels")[attr(x, "order") == 1] %in% LHS_Y)))) stop("Outcomes cannot directly predict one another")

  ## check if covariates predict each other
  if (length(formulas[[1]]) > 1) {
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

  ## process family variable inputs
  family <- process_family(family)

  LHSs <- list(LHS_Z, LHS_X, LHS_Y)

  ## check that variables in families 1,2,3,4 have a dispersion parameter
  for (j in 1:3) {
    if (is.numeric(family[[j]])) {
      for (i in seq_along(family[[j]])) {
        if (family[[j]][i] <= 4 && is.null(pars[[LHSs[[j]][i]]]$phi))
          stop(paste(LHSs[[j]][i], "needs a phi parameter"))
        if (family[[j]][i] == 2 && is.null(pars[[LHSs[[j]][i]]]$par2))
          stop(paste(LHSs[[j]][i], "needs a par2 parameter for degrees of freedom"))
      }
    }
    else if (length(family[[j]]) > 0 && is(family[[j]][1], "causl_family")) {
      for (i in seq_along(family[[j]])) {
        if (!all(family[[j]][i]$pars %in% names(pars[[LHSs[[j]][i]]]))) {
          miss <- family[[j]][i]$pars[!family[[j]][i]$pars %in% names(pars[[LHSs[[j]][i]]])]
          stop(paste0("Parameters ", paste(miss, collapse=", "), " missing for variable ", LHS[[j]][i]))
        }
      }
    }
  }
  for (i in seq_along(famX))
    if (famX[i] <= 4 && is.null(pars[[LHS_X[i]]]$phi))
      stop(paste(LHS_X[i], "needs a phi parameter"))
  for (i in seq_along(famY))
    if (famY[i] <= 4 && is.null(pars[[LHS_Y[i]]]$phi))
      stop(paste(LHS_Y[i], "needs a phi parameter"))


  if (ordering) {

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
          formulas[[4]][[i]] <- rep(formulas[[4]][[i]], dZ+i-1)
        }
        else {
          formulas[[4]][[i]] <- rep(list(formulas[[4]][[i]]), dZ+i-1)
        }
        lhs(formulas[[4]][[i]]) <- c(LHS_Z[rank(order[seq_len(dZ)])], LHS_Y[rank(order[dZ+dX+seq_len(i-1)])])
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
      family[[4]] <- mapply(function(x, y) rep(x, y), family[[4]], dZ+seq_len(dY)-1)
    }

    ## get parameters in right format
    if (!setequal(names(pars[[kwd]]), LHS_Y)) {
      if ("beta" %in% names(pars[[kwd]])) {
        pars[[kwd]] <- rep(list(list(list(beta=pars[[kwd]]$beta))), length(LHS_Y))
        names(pars[[kwd]]) <- LHS_Y
        if (dZ > 1 || dY > 1) {
          for (i in seq_len(dY)) {
            vnm <- vars[order[dZ+dX+i]]
            pars[[kwd]][[vnm]] <- rep(pars[[kwd]][[vnm]], dZ+i-1)
            names(pars[[kwd]][[vnm]]) <- c(LHS_Z, LHS_Y[rank(order[dZ+dX+seq_len(i-1)])])
          }
        }
      }
      else {
        nrep <- setdiff(LHS_Y, names(pars[[kwd]]))
        if (length(nrep) > 0) stop(paste0("Variable ", paste(nrep, collapse=", "),
                                          "not represented in the copula parameters list"))
        rep <- setdiff(names(pars[[kwd]]), LHS_Y)
        if (length(rep) > 0) stop(paste0("Variable ", paste(nrep, collapse=", "),
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
    if (any(c(LHS_Z,LHS_Y) %in% rhs_vars(formulas[[4]])[[1]])) stop("copula cannot depend upon Z or Y variables")
  }

  ## set up link functions
  link <- link_setup(link, family[1:3], vars=list(LHS_Z,LHS_X,LHS_Y))

  ## other useful variables
  output <- c(LHS_Z, LHS_Y)
  vars <- c(LHS_Z, LHS_X, LHS_Y)

  return(list(formulas=formulas, pars=pars, family=family, link=link,
              LHSs=list(LHS_Z=LHS_Z, LHS_X=LHS_X, LHS_Y=LHS_Y), kwd=kwd,
              dim=dims, in_form=in_form, vars=vars, output=output,
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
