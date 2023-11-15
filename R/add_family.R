##' @export
print.causl_family <- function (x, ...) {
  custom <- isTRUE("cust_family" %in% class(x))
  cat(ifelse(custom, "Custom family", "Family"), ": ", x$name, "\n", sep="")
  cat(paste0("Parameters: ", paste(x$pars, collapse=", ")), "\n")

  invisible(x)
}
