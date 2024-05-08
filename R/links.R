##' Set up link functions
##'
##' @param link input given to `msm_samp()` or `causalSamp()`
##' @param family the list of families for `Z`,`X` and `Y` variables
##' @param vars a list of vectors of variable names with the same structure as `family`
##' @param sources list of links for parametric families
##' @param fam_list list of data frames in the same format as `family_vals`
##'
##' @export
link_setup <- function(link, family, vars, sources=links_list,
                       fam_list=list(family_vals)) {

  if (!missing(vars) && !all(lengths(vars) == lengths(family))) stop("length of variable names vector does not match number of families provided")

  ## concatenate list of courses
  lk_lsts <- sources # do.call(c, sources)
  fm_lsts <- do.call(c, fam_list)

  ## set up output
  link_out <- lapply(family, function(x) {
    if (length(x) > 0 && is(x[[1]], "causl_family")) sapply(x, function(y) y$name)
    else fm_lsts$family[match(x, fm_lsts$val)]
  })
  if (any(is.na(unlist(link_out)))) stop("Invalid family variable specified")
  link_out <- lapply(link_out, function(x) sapply(x, function(y) lk_lsts[[y]][1]))

  if (any(is.null(unlist(link_out)))) stop("Error in family specification for link functions")

  if (!missing(vars)) {
    ## set names to match vars
    for (i in seq_along(link_out)) names(link_out[[i]]) <- vars[[i]]
  }

  ## if no link argument supplied, then just return this list
  if (missing(link) || is.null(link) || is.null(link[[1]])) return(link_out)

  tmp <- unlist(lk_lsts)

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
    #   for (i in seq_along(lk_lsts)) {
    #     lki <- link[mask && fam_nm == nms[i]]
    #     linknm <- pmatch(lki, lk_lsts[[i]])
    #     lki2 <- lk_lsts[[i]][linknm]
    #     if (any(is.na(lki2))) stop("link ", paste(lki[which(is.na(lki2))], sep=", "), " not supported")
    #     lki[]
    #   }
    #   chk <- wh_set
    # }
  }
  else stop("variable names must be provided to match using them")

  # ## matching
  # fam_nm <- fm_lsts$family[fm_lsts$val==fams]
  # nms <- names(lk_lsts)
  # for (i in seq_along(lk_lsts)) {
  #   lki <- link[fam_nm == nms[i]]
  #   linknm <- pmatch(lki, lk_lsts[[i]])
  #   lki2 <- lk_lsts[[i]][linknm]
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
links_list <- list(
  gaussian = c("identity", "inverse", "log"),
  t = c("identity", "inverse", "log"),
  Gamma = c("log", "inverse", "identity"),
  beta = c("logit", "probit"),
  binomial = c("logit", "probit", "log"),
  lognormal = c("exp", "identity"),
  categorical = c("logit"),
  ordinal = c("logit")
)

##' @describeIn links_list Old name
##' @format `linksList` is the old name for `links_list`
linksList <- links_list
