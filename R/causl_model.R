##' Define a `causl_model` object
##'
##' This defines a `causl_model` object, that can be used either for simulation
##' or inference.
##'
##' @inheritParams rfrugalParam
##' @param kwd word used for copula formula and parameters
##'
##' @details
##' The components `formulas` and `family` must both be specified, and have
##' matching lengths.  If `pars` is specified, then the model can be used for
##' simulation and inference, if not then only for inference.  `link` is optional,
##' and if not specified default links will be used.
##'
##' @export
causl_model <- function (formulas, family, pars, link, dat=NULL, method="inversion",
                         kwd="cop", control=list()) {


  out <- process_inputs(formulas=formulas, family=family, pars=pars,
                        link=link, dat=dat, kwd=kwd, method=method)

  class(out) <- "causl_model"

  return(out)
}


##' Display output from `causl_model`
##'
##' @param x an object of class `causl_model`
##' @param ... additional arguments (not used)
##'
##' @exportS3Method base::print
print.causl_model <- function (x, ...) {
  forms <- unlist(x$formulas)
  cat("Frugal model with variables ", paste(x$vars, collapse=", "), "\n")

  invisible(x)
}

##' Generic for modify
##'
##' @param x object to be modified
##' @param ... other arguments
##'
##' @examples
##' cm <- causl_model(formulas = list(Z ~ 1, Y ~ 1, ~ 1),
##'                   family = list(3, 1, 1),
##'                   pars = list(Z = list(beta = 1, phi = 1),
##'                               Y = list(beta = 0, phi = 1),
##'                               cop = list(beta = 1)))
##' modify(cm, pars = list(cop = list(beta = 1.5)))
##' cm$pars
##'
##' @export
modify <- function (x, ...) {
  NextMethod("modify")
}

modify.default <- function (x, ...) {
  args <- list(...)
  for (i in seq_along(args)) x[[names(args)[i]]] <- args[[i]]
  return(x)
}


##' Modify `causl_model` object
##'
##' Change one or more components of a `causl_model` object.
##'
##' @inheritParams causl_model
##' @param x an object of class `causl_model`
##' @param over logical: should components be added/modified or entirely over-written?
##'
##' This function can be used to modify
##'
##' @export
modify.causl_model <- function (x, over=FALSE, formulas, family, pars, link, dat, method,
                                kwd) {
  if (!is(x, "causl_model")) stop("Must include an object of class 'causl_model'")

  if (missing(formulas) && missing(family) && missing(pars) && missing(link) &&
      missing(dat) && missing(method) && missing(kwd)) return(x)
  if (missing(formulas)) formulas <- x$formulas
  if (missing(family)) family <- x$family
  if (missing(pars)) pars <- x$pars
  else {
    cpars <- x$pars
    cpars[names(pars)] <- pars
    pars <- cpars
  }
  if (missing(link)) link <- x$link
  if (missing(dat)) dat <- x$dat
  if (missing(method)) method <- x$method
  if (missing(kwd)) kwd <- x$kwd

  out <- process_inputs(formulas=formulas, family=family, pars=pars,
                        link=link, dat=dat, method=method, kwd=kwd)

  return(out)
}

# ##' @export
# modify_trt_prob.causl_model <- function (cm, prob=0) {
#   if (length(cm$formulas[[2]]) > 1) stop("Must be a single treatment")
#   if (!(isTRUE(all.equal(cm$family[[2]], 5)) ||
#         isTRUE(all.equal(cm$family[[2]], 0)) ||
#         isTRUE(all.equal(cm$family[[2]], "binomial")) ||
#         isTRUE(all.equal(cm$family[[2]], binomial_causl_fam())))) stop("Must be a single binary treatment")
#
#   cm0 <- cm
#   cm0$formulas[[2]][[1]] <- update.formula(cm0$formulas[[2]][[1]], ". ~ 1")
#
#   if (prob == 0) bd <- logit(1/(1e3*n))
#   else if (prob == 1) bd <- -logit(1/(1e3*n))
#   else stop("Probability should be 0 or 1")
#
#   ## set up new parameter values
#   new_pars0 <- list(list(beta = bd))
#   t_var <- lhs(cm$formulas[[2]])
#   names(new_pars0) <- t_var
#   cm0 <- modify.causl_model(cm0, pars = new_pars0)
#
#   return(cm0)
# }
