##' New function to sort out links and families simultaneously
##'
##' @inheritParams process_family
##' @inheritParams link_setup
##'
##' Given a list for families, another for links, and a vector of dimensions,
##' will put the families and links into a standard format.
##'
##' @export
process_family_link <- function (family, link, dims, func_return=get_family) {

  nU <- length(family) - 1
  if (length(family) != length(dims)) stop("Should be a family entry for each set of formulas")

  ## check family is a list, or get into list form if not
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
  else if (length(family) == nU+1) {  ## always true given dfn of nU
    if (sum(dims[seq_len(nU)]) > nU) stop("Mismatch in family and formulae specification")
    family <- as.list(family)
    lens <- rep(1, nU+1)
  }
  else stop(paste0("'family' should be a list, or vector of length ", nU+1))

  ## now check link is in list format, or else modify to be a list
  if (missing(link) || is.null(link)) {
    link <- vector("list", nU)
    for (i in seq_len(nU)) {
      link[[i]] <- character(length(family[[i]]))
      if (length(family[[i]]) > 0) {
        if (is(family[[i]][[1]], "causl_family")) link[[i]] <- sapply(family[[i]], function(x) x$link)
        else if (is.numeric(family[[i]])) {
          link[[i]] <- sapply(family[[i]], function (x) {
            wh <- which(x == family_vals[,1])
            if (length(wh) != 1 || is.na(wh)) stop(paste0("Invalid family ", x))
            links_list[[family_vals[wh,2]]][1]
          })
        }
        else if (is.character(family[[i]])) {
          link[[i]] <- sapply(family[[i]], function (x) links_list[[x]][1])
        }
      }
    }
  }
  else if (is.list(link)) {
    if (!all(lengths(link) == lengths(family)[seq_len(nU)])) stop("Lengths of links do not match those of families")
  }
  else cat(paste0("'link' is vector of length ", length(link)))

  ## deal with family names and causl_fam functions
  for (i in seq_len(nU)) {
    if (lens[i] > 0) {
      family[[i]] <- fam_chk(family[[i]], lens[i], link[[i]], func_return = func_return)
    }
  }

  ## check copula families are valid
  # if (!all(unlist(family[1:3]) %in% family_vals$val)) stop("Invalid family specification")
  if (!all(unlist(family[[nU+1]]) %in% copula_vals$val)) stop("Invalid copula specification")

  return(list(family=family, link=link))
}
