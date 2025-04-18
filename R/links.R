##' Set up link functions
##'
##' @param link input from `causl_model` object or similar
##' @param family the list of families for random variables
##' @param vars a list of vectors of variable names with the same structure as `family`
##' @param sources list of links for parametric families
##' @param fam_list list of data frames in the same format as `family_vals`
##'
##' @description
##' `family` and `vars` should have the same structure if `vars` is
##' specified.
##'
##' @export
link_setup <- function(link, family, vars, sources=links_list,
                       fam_list=list(family_vals)) {

  if (!missing(vars)) {
    if (!all(lengths(vars) == lengths(family))) stop("length of variable names vector does not match number of families provided")
    if (!missing(link) && isTRUE(any(!(names(unlist(link)) %in% unlist(vars))))) stop("Link named for variable not in the model")
  }

  ## concatenate list of sources
  lk_lsts <- sources # do.call(c, sources)
  all_lks <- unlist(lk_lsts)
  fm_lsts <- do.call(c, fam_list)

  ## set up output
  is_cf <- sapply(family, function(x) (length(x) > 0 && is(x[[1]], "causl_family")))
  link_out <- lapply(family, function(x) {
    if (length(x) > 0 && is(x[[1]], "causl_family")) sapply(x, function(y) y$link)
    else fm_lsts$family[match(x, fm_lsts$val)]
  })
  ## check that links are valid for causl_family specifications
  if (any(is_cf)) {
    cf_links <- link_out[is_cf]
    cfs <- family[is_cf]
    vld <- lapply(cf_links, function (x) x %in% all_lks)
    vld <- mapply(function (x,y,z)
      mapply(function (a,b) b %in% names(a$custom_links), a=x,b=y) | z,
      cfs, cf_links, vld)
    if (!all(unlist(vld))) stop("Some links not valid")
  }
  if (any(is.na(unlist(link_out)))) stop("Invalid family variable specified")
  link_out[!is_cf] <- lapply(link_out[!is_cf], function(x) {
    sapply(x, function(y) lk_lsts[[y]][1])
  })

  if (any(is.null(unlist(link_out)))) stop("Error in family specification for link functions")

  if (!missing(vars)) {
    ## set names to match vars
    for (i in seq_along(link_out)) names(link_out[[i]]) <- vars[[i]]
  }

  ## if no link argument supplied, then just return this list
  if (missing(link) || is.null(link)) return(link_out)

  tmp <- unlist(lk_lsts)


  if (is.null(names(unlist(link)))) {
    if (!all(lengths(family) == lengths(link))) stop("If no names provided, lengths of 'link' must match those of 'family'")
    ## if lengths the same, assume in same position as family variables
    for (i in seq_along(link_out)[lengths(family) > 0]) {
      if (is(family[[i]][[1]], "causl_family")) {
        link_out[[i]] <- sapply(family[[i]], link)
        next
      }
      else {
        link_out[[i]][] <- tmp[pmatch(link[[i]], tmp)]
        valid <- mapply(
          function(x,y) {
            fms <- family_vals$family[match(y,family_vals$val)]
            match(x, lk_lsts[[fms]])
          },
          link_out[[i]], family[[i]])
        if (any(is.na(valid))) stop(paste0("Invalid link supplied for families: ",
                                           paste(family[[i]][is.na(valid)], collapse=", ")))
      }
      if (any(is.na(link_out[[i]]))) stop("link not properly matched")
    }
    return(link_out)
  }
  else {
    ## unlist provided links
    link <- unlist(link)

    ## now add in any modifications made in 'link'
    ## see if names can be used
    if (!missing(vars)) {
      if (!isTRUE(all(names(link) %in% unlist(vars)))) stop("Some names not contained in 'vars'")
    }

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
    for (i in seq_along(family)) {
      if (is(family[[i]][[1]], "causl_family")) {
        for (j in seq_along(family[[i]])) {
          if (link(family[[i]][[j]]) != link_out[[i]][j]) {
            warning("Provided link does not match that in family. Using family link")
            link_out[[i]][[j]] <- link(family[[i]][[j]])
          }
        }
      }
      else {
        valid <- mapply(
          function(x,y) {
            fms <- family_vals$family[match(y,family_vals$val)]
            match(x, lk_lsts[[fms]])
          },
          link_out[[i]], family[[i]])
        if (any(is.na(valid))) stop(paste0("Invalid link supplied for families: ",
                                           paste(family[[i]][is.na(valid)], collapse=", ")))
      }
    }
    if (any(is.na(unlist(link_out)))) stop("some links not matched")
  }


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
