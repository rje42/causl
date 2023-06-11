##' Get LHSs of list of formulae
##'
##' @param formulas list of formulae
##'
##' @details Returns character vector containing left-hand sides of a list of
##' formulae.
##'
##' @export
lhs <- function (formulas) {
  if (!is.list(formulas)) formulas <- list(formulas)
  term <- lapply(formulas, terms)

  ## get list of variables
  vars <- lapply(term, function(x) attr(x, which="variables"))
  resp <- rep(NA_character_, length(vars)) #character(length(vars))
  if (any(lengths(vars) >= 2)) {
    mis <- (sapply(term, attr, which="response") != 1)
    if (!all(mis)) resp[!mis & lengths(vars) >= 2] <- sapply(vars[!mis& lengths(vars) >= 2], function(x) as.character(x[[2]]))
  }
  # mis <- (sapply(term, attr, which="response") != 1)
  # resp[mis] <- rep()

  resp
}


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
topOrd <- function (A) {
  ord <- integer(0)
  actv <- seq_len(nrow(A))

  while(length(actv) > 0) {
    npa <- rowSums(A[actv,,drop=FALSE])
    wh0 <- which(npa == 0)

    if (length(wh0) == 0) return(NA)
    ord <- c(ord, actv[wh0])
    A[,actv[wh0]] <- 0
    actv <- actv[-wh0]
  }

  ord
}
