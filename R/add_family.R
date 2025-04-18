##' @export
print.causl_family <- function (x, ...) {
  custom <- isTRUE("cust_family" %in% class(x))
  cat(ifelse(custom, "Custom family", "Family"), ": ", x$name, "\n", sep="")
  cat(paste0("Parameters: ", paste(x$pars, collapse=", ")), "\n")
  cat(paste0("Link: ", paste(link_name(x), collapse=", ")), "\n")

  invisible(x)
}

##' @export
print.causl_copula <- function (x, ...) {
  custom <- isTRUE("cust_family" %in% class(x))
  cat(ifelse(custom, "Custom family", "Family"), ": ", x$name, " copula\n", sep="")
  cat(paste0("Parameters: ", paste(x$pars, collapse=", ")), "\n")
  cat(paste0("Link: ", paste(link_name(x), collapse=", ")), "\n")

  invisible(x)
}
