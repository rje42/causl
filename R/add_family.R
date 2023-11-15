##' Add a new parametric family
##'
##' Include a new parametric family by specifying its density, quantiles and a function to simulate with
##'
##' @param name string giving name for family
##' @param dens,quan,sim density, quantile and simulation functions (see details)
##' @param pars named vector of parameters used to evaluate function
##' @param links character vector of link functions to use with family
##' @param number optional number to specify
##'
## @importFrom rlang get_env
##'
## @export
new_family <- function (name, ddist, qdist, rdist, pdist, pars=c(x=0, mu=0, phi=1, par2=1), links) {
  # if (missing(number)) {
  #   number <- max(familyVals$val) + 1
  # }
  # if (number %in% familyVals$val) {
  #   number <- max(familyVals$val) + 1
  #   message(paste0("Number already taken, so assigning ", number))
  # }
  if (name %in% familyVals$family) stop(paste0("Family '", name, "' already defined"))

  # print(rlang::env_parent())
  # return(invisible(NULL))
  # this_pkg <- rlang::get_env(new_family)
  # labs <- paste(c("d","q","r"), name, sep="")

  # unlockBinding(sym = "familyVals", env = this_pkg)
  # familyVals <<- rbind(familyVals, data.frame(val=number, family=name))
  # lockBinding(sym = "familyVals", env = this_pkg)
  # print(familyVals)
  # return()
  # linksList[[name]] <<- links

  labs <- names(pars)
  rlabs <- labs %in% c("mu", "phi", "par2")
  dlabs <- labs %in% "x" | rlabs
  qlabs <- labs %in% "p" | rlabs
  nec_pars <- names(labs[rlabs])

  ##
  dens <- do.call(ddist, pars[dlabs])
  ldens <- do.call(ddist, c(pars[dlabs], list(log=TRUE)))
  quan <- do.call(qdist, pars[qlabs])
  sim <- do.call(rdist, c(list(n=10), pars[rlabs]))
  probs <- do.call(pdist, pars[dlabs])
  if (length(dens) != 1 || is.na(dens) || is.null(dens)) stop("Density function not valid")
  if (length(ldens) != 1 || !isTRUE(all.equal(exp(ldens), dens))) stop("Log-density function not evalusted correctly")
  if (length(quan) != 1 || is.na(quan) || is.null(quan)) stop("Quantile function not valid")
  if (length(sim) != 10 || any(is.na(sim)) || any(is.null(sim))) stop("Simulation function not valid")
  if (length(probs) != 1 || is.na(probs) || is.null(probs)) stop("CDF function not valid")

  out <- list(name=name, ddist=ddist, qdist=qdist, rdist=rdist, pdist=pdist, pars=nec_pars)
  class(out) <- c("cust_family", "causl_family")

  return(out)
}

print.causl_family <- function (x, ...) {
  custom <- isTRUE("cust_family" %in% class(x))
  cat(ifelse(custom, "Custom family", "Family"), ": ", x$name, "\n", sep="")
  cat(paste0("Parameters: ", paste(x$pars, collapse=", ")), "\n")

  invisible(x)
}
