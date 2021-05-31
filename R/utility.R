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
