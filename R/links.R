##' Set up link functions
##'
##' @param link the input given to causalSamp()
##' @param family the list of families for Z,X and Y variables
##' @param vars a list of vectors of variable names with the same structure as \code{family}
##'
##' @export
link_setup <- function(link, family, vars) {

  if (!missing(vars) && !all(lengths(vars) == lengths(family))) stop("length of variable names vector does not match number of families provided")

  ## set up output
  link_out <- lapply(family, function(x) {
    if (length(x) > 0 && is(x[[1]], "causl_family")) sapply(x, function(y) y$name)
    else familyVals$family[match(x, familyVals$val)]
  })
  if (any(is.na(unlist(link_out)))) stop("Invalid family variable specified")
  link_out <- lapply(link_out, function(x) sapply(x, function(y) linksList[[y]][1]))

  if (!missing(vars)) {
    ## set names to match vars
    nms_list <- relist(unlist(vars), family)
    for (i in seq_along(link_out)) names(link_out[[i]]) <- nms_list[[i]]
  }

  ## if no link argument supplied, then just return this list
  if (missing(link) || is.null(link)) return(link_out)

  tmp <- unlist(linksList)

  ## now add in any modifications made in 'link'
  if (is.list(link)) {
    if (all(lengths(link) == lengths(family))) {
      ## if lengths the same, assume in same position as family variables
      for (i in seq_along(link_out)) {
        link_out[[i]][] <- tmp[pmatch(link[[i]], tmp)]
        if (any(is.na(link_out[[i]]))) stop("link not properly matched")
      }
      return(link_out)
    }
  }

  ## otherwise try to use names to deduce which is which
  if (!missing(vars)) {
    link <- unlist(link)
    if (is.null(names(link))) stop("names must be supplied for links in order to match with them")
    nms <- names(link)
    wh_set <- subsetmatch(as.list(nms), vars)
    if (any(is.na(wh_set))) stop("names must be supplied for links in order to match with them")

    ## now get particular entry in vector
    wh_val <- NA*wh_set
    for (i in seq_along(wh_set)) {
      wh_val[i] <- match(nms[i], vars[[wh_set[i]]])
      if (is.na(wh_val[i])) stop("some links not matched")
      link_out[[wh_set[i]]][wh_val[i]] <- tmp[pmatch(link[i], tmp)]
    }
    if (any(is.na(unlist(link_out)))) stop("some links not matched")
    #
    #
    # for (j in seq_along(vars)) {
    #   mask <- wh_set == j
    #
    #   for (i in seq_along(linksList)) {
    #     lki <- link[mask && fam_nm == nms[i]]
    #     linknm <- pmatch(lki, linksList[[i]])
    #     lki2 <- linksList[[i]][linknm]
    #     if (any(is.na(lki2))) stop("link ", paste(lki[which(is.na(lki2))], sep=", "), " not supported")
    #     lki[]
    #   }
    #   chk <- wh_set
    # }
  }
  else stop("variable names must be provided to match using them")

  # ## matching
  # fam_nm <- familyVals$family[familyVals$val==fams]
  # nms <- names(linksList)
  # for (i in seq_along(linksList)) {
  #   lki <- link[fam_nm == nms[i]]
  #   linknm <- pmatch(lki, linksList[[i]])
  #   lki2 <- linksList[[i]][linknm]
  #   if (any(is.na(lki2))) stop("link ", paste(lki[which(is.na(lki2))], sep=", "), " not supported")
  #

  return(link_out)
}

##' List of links available for each parametric family
##'
##' This is a named list whose entries are character vectors of allowed link
##' functions.
##'
##' @export
linksList <- list(
  gaussian = c("identity", "inverse", "log"),
  t = c("identity", "inverse", "log"),
  Gamma = c("log", "inverse", "identity"),
  beta = c("logit", "probit"),
  binomial = c("logit", "probit", "log"),
  lognormal = c("exp", "identity")
)
