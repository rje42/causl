##' Tools for manipulating formulas
##'
##' @param formulas list of formulae
##'
##' @details `lhs` returns a character vector containing left-hand sides of a
##' list of formulae.  If `surv=TRUE` then two responses are returned in the
##' event of the left-hand side being a valid `Surv` object.  `lhs<-` allows
##' one to assign the left-hand sides of variables in the obvious way.
##'
##' `tidy_formulas` ensures that all formulae in a list have a left hand side,
##' by giving them names of the form `Vn` where `n` is some positive integer. The
##' prefix `V` can be changed using the argument `prefix`.
##'
##' `rhs_vars` extracts all the variables used on the right-hand sides of a list
##' of formulas.
##'
##' @name formula_tools
NULL


## To do: put in check for invalid name and then extract multiple variables if a Surv object
##

##' @describeIn formula_tools Obtain left-hand sides from list of formulas
##' @param surv logical indicating whether to treat as survey data
##' @export
lhs <- function (formulas, surv=FALSE) {
  if (!is.list(formulas)) formulas <- list(formulas)
  term <- lapply(formulas, terms)

  ## get list of variables
  vars <- lapply(term, function(x) attr(x, which="variables"))
  resp <- rep(NA_character_, length(vars)) #character(length(vars))
  if (any(lengths(vars) >= 2)) {
    mis <- (sapply(term, attr, which="response") != 1)
    if (!all(mis)) {
      wh <- !mis & lengths(vars) >= 2
      new_vals <- lapply(vars[wh], function(x) as.character(x[[2]]))
      resp[wh & lengths(new_vals) == 1] <- unlist(new_vals[lengths(new_vals) == 1])

      if (any(lengths(new_vals) != 1 & lengths(new_vals) != 3)) stop("Not a valid left-hand side")

      if (any(lengths(new_vals) == 3)) {
        resp <- as.list(resp)
        wh3 <- which(lengths(new_vals) == 3)
        vals1 <- sapply(new_vals, function(x) x[[1]])
        if (any(vals1 != "Surv")) stop("Non-standard syntax must be a 'Surv' object")
        vals_out <- lapply(new_vals[wh3], function(x) x[2:3])
        resp[wh3] <- vals_out
      }
    }

    #
    # if (any(na.omit(grepl("[$&+:;=?@|'<>^*%!-]", resp))) || (any(na.omit(grepl("[,()]", resp))) && !surv)) stop("Invalid variables in formula")
    # if (surv && any(na.omit(grepl("[,()]", resp)))) {
    #   wh <- which(grepl("[,()]", resp))
    #   vals <- lapply(resp, function(x) strsplit(x, "[(,)]")[[1]])
    #   lVs <- lengths(vals)
    #   vals1 <- sapply(vals, function(x) x[[1]])
    #   if (any(vals1 != "Surv")) stop("Non-standard syntax must be a 'Surv' object")
    #   if (any(lVs != 3)) stop("'Surv' objects must have exactly two arguments")
    #   resp <- as.list(resp)
    #   vals_out <- mapply(c, vals[[2]], vals[[3]])
    #   resp[wh] <- vals
    # }
  }
  # mis <- (sapply(term, attr, which="response") != 1)
  # resp[mis] <- rep()

  return(resp)
}

##' @describeIn formula_tools Assign left-hand sides to list of formulas
##' @param value character vector to assign
##' @export
`lhs<-` <- function (formulas, value) {
  if (length(formulas) < length(value)) stop("Must be as many formulas as LHSs")

  if (!is.list(formulas)) formulas <- list(formulas)
  formulas <- mapply(function(x,y) update.formula(x, paste0(y, " ~ .")), formulas, value, SIMPLIFY = FALSE)

  return(formulas)
}

##' @describeIn formula_tools Extract variables from right-hand sides
##' @export
rhs_vars <- function (formulas) {
  if (!is.list(formulas)) formulas <- list(formulas)
  term <- lapply(formulas, terms)
  wres <- (sapply(term, attr, which="response") == 1)
  nores <- (sapply(term, attr, which="response") == 0)

  vars <- lapply(term, function(x) attr(x, which="variables"))

  ## get list of variables
  vars[wres] <- lapply(vars[wres], function(x) as.character(x[-(1:2)]))
  vars[nores] <- lapply(vars[nores], function(x) as.character(x[-1]))

  # if (any(lengths(vars) >= 2)) resp[lengths(vars) >= 3] <- lapply(vars[lengths(vars) >= 3], function(x) as.character(x))
  # resp[mis] <- rep()

  vars
}

##' Combine multiple formulas
##'
##' Take collection of formulae and create one formula with all variables on the
##' right-hand side of any of the originals.
##'
##' @param formulas list of formulas to merge
##'
##' @export
merge_formulas <- function (formulas) {
  if (!is.list(formulas)) formulas <- list(formulas)
  formulas <- unlist(formulas)
  term <- lapply(formulas, terms)
  LHSs <- lhs(formulas)
  # if (any(is.na(LHSs))) warning("some formulae did not have a left-hand side, output will be unlabelled")

  elems_lst <- lapply(term, attr, which="term.labels")
  elems <- unlist(elems_lst)
  ords <- unlist(lapply(term, attr, which="order"))
  intcpt <- sapply(term, attr, which="intercept")

  if (length(elems) != length(ords)) stop("length of labels and orders do not match")

  ## remove duplicates and reorder
  ords <- ords[!duplicated(elems)]
  elems <- elems[!duplicated(elems)]

  elems <- elems[order(ords)]
  ords <- ords[order(ords)]

  wh <- vector(mode="list", length=length(formulas))
  if (!any(is.na(LHSs))) names(wh) <- LHSs

  reforms <- formulas

  ## now get columns relevant for each formula
  for (i in seq_along(formulas)) {
    wh[[i]] <- match(elems_lst[[i]], elems) + 1
    if (intcpt[i] == 1) wh[[i]] <- c(1, wh[[i]])

    ## reorder original formulae to have same order as merged one
    reforms[[i]] <- formula(paste0(LHSs[i], " ~ ", ifelse(intcpt[i]==1, "", "0 + "),
                                   ifelse(intcpt[i]==1 && length(wh[[i]])==1, "1", ""),
                                   paste(elems[sort.int(wh[[i]])-1], collapse=" + ")))
  }
  if (length(elems) == 0) elems = "1"

  full_form <- formula(paste("~", paste(unlist(elems), collapse=" + ")))
  list(formula=full_form, wh=wh, reforms=reforms, old_forms=term)
}

## Get topological order from an adjacency matrix
topOrd <- function (A, from_top = FALSE) {
  ord <- integer(0)
  actv <- seq_len(nrow(A))

  if (from_top) {
    while(length(actv) > 0) {
      npa <- rowSums(A[actv,,drop=FALSE])
      wh0 <- which(npa == 0)

      if (length(wh0) == 0) return(NA)
      ord <- c(ord, actv[wh0])
      A[,actv[wh0]] <- 0
      actv <- actv[-wh0]
    }
  }
  else {
    while(length(actv) > 0) {
      npa <- colSums(A[,actv,drop=FALSE])
      wh0 <- which(npa == 0)

      if (length(wh0) == 0) return(NA)
      ord <- c(actv[wh0], ord)
      A[actv[wh0],] <- 0
      actv <- actv[-wh0]
    }
  }

  ord
}

## @inheritParams fitCausal
##' @param kwd string used to denote copula
##' @param prefix string to begin each new variable name
##' @describeIn formula_tools Tidy up formulae
##' @export
tidy_formulas <- function(formulas, kwd, prefix="V") {
  forms <- formulas
  nf <- length(forms)

  ## make sure all the univariate formulae have left-hand sides
  wh <- which(lengths(forms) < 3)

  if (any(lengths(forms) < 3)) {
    if (last(wh) == nf) wh2 <- wh[-length(wh)]
    else wh2 <- wh
    for (i in seq_along(wh2)) {
      tmp_form <- formula(paste(prefix, i, " ~ .", sep=""))
      forms[[wh[i]]] <- update.formula(forms[[wh[i]]], tmp_form)
    }
  }

  ## give copula formula the keyword as its left-hand side
  if (length(wh) > 0 && last(wh) == nf) {
    tmp_form <- formula(paste(kwd, "~ ."))
    forms[[nf]] <- update.formula(forms[[nf]], tmp_form)
  }
  else if (lhs(last(forms)) != kwd) stop("Error: keyword does not match left-hand side of copula formula")

  return(forms)
}
