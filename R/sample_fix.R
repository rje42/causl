##' Get maximum weight for each segment of a distribution
##'
##' @param pars list with all regression parameters
##' @param forms_X formulae for treatments
##' @param fam_X,fam_Z vector of families for treatments and covariates
##' @param LHS_Z variables in covariates
##' @param qden density of proposals
##' @param ranges matrix of ranges of segments
##' @param link link functions for treatments
##' @param ... not currently used
##'
get_max_weights <- function (pars, forms_X, fam_X, qden, fam_Z, LHS_Z, ranges, link, ...) {
  # if (!missing(link)) warning("Link argument is ignored")

  ## get unique list of segments
  # m_wts <- unique(Z0q)
  # m_wts$max_wt <- NA

  LHS_X <- lhs(forms_X)
  dX <- length(LHS_X)
  full_form <- merge_formulas(forms_X)

  ## get states for discrete treatments
  if (length(fam_X) > 0 && is(fam_X[[1]], "causl_fam")) {
    famChr <- sapply(fam_X, function(x) x$name)
    disX <- sum(famChr %in% "binomial")
    # disZ <- sum(fam_Z == 5)
    combs <- expand.grid(rep(list(0:1), disX))
    names(combs) <- names(LHS_X[famChr %in% "binomial"])
  }
  else {
    disX <- sum(fam_X == 5)
    # disZ <- sum(fam_Z == 5)
    combs <- expand.grid(rep(list(0:1), disX))
    names(combs) <- names(LHS_X[fam_X == 5])
  }

  ## build dummy data frame
  Zvars <- setdiff(unique(unlist(rhs_vars(forms_X))), LHS_X)
  fam_Z <- fam_Z[match(LHS_Z, Zvars, nomatch = 0L)]  # get only relevant Z variables
  nms <- c(LHS_X, Zvars)
  args <- rep(list(NA), length(nms))
  names(args) <- nms
  out <- do.call(data.frame, args)


  ## get matrix form for parameters
  msks <- initializeParams2(formulas=forms_X, family=fam_X, full_form=full_form,
                            only_masks=TRUE, nc=0, inc_cop=FALSE)
  beta <- beta_m <- msks$beta_m
  phi <- phi_m <- msks$phi_m

  for (i in seq_len(ncol(beta_m))) {
    beta[full_form$wh[[i]],i] <- pars[[LHS_X[i]]]$beta
    if (fam_X[i] %in% c(1:3,6)) phi[i] <- pars[[LHS_X[i]]]$phi
  }
  if (any(beta != 0 & beta_m == 0) ||
      any(phi != 0 & phi_m == 0)) stop("Mask problem")

  disc <- c(fam_X==5, fam_Z==5)

  func <- function (x, dis_val=integer(0)) {
    ## put variables in data frame
    ZX <- out
    ZX[1, !disc] <- x
    ZX[1, disc] <- dis_val
    ZX2 <- ZX[,seq_along(fam_X),drop=FALSE]
    # X <- x[seq_len(fam_X)]
    # Z <- x[-seq_len(fam_X)]
    # beta

    mm <- model.matrix(full_form$formula, ZX)
    eta <- mm %*% beta

    out <- get_X_density(dat = ZX2, eta = eta, phi = phi, par2=0*phi,
                         qden = qden, family=fam_X, link = link)

    -out
  }

  # if (length(ranges) > 1) rg <-
  # else if (length(ranges) == 1) rg <- list(ranges[[1]][1],ranges[[1]][2])
  # else rg <- list(numeric(0), numeric(0))

  ## need to make move in several directions
  best <- list(value=func(rep(0.1,dX+length(fam_Z)-sum(disc)),rep(0,sum(disc))))
  n_strts <- 10
  for (i in seq_len(n_strts)) {
    strt <- runif(nrow(ranges), min = ranges[,1], max = ranges[,2])
    strt <- c(rnorm(dX-disX)/10, strt)
    if (nrow(combs) > 0) {
      for (dx in seq_len(nrow(combs))) {
        opt <- tryCatch(optim(par=strt, fn=func, dis_val = combs[dx,], method = "L-BFGS-B",
                              lower=c(rep(-Inf, dX-disX),ranges[,1]), upper=c(rep(Inf, dX-disX),ranges[,2])),
                        error=function(e) return(list(value=Inf)))
        if (opt$value < best$value) best <- opt
      }
    }
    else {
      opt <- tryCatch(optim(par=strt, fn=func, dis_val = integer(0), method = "L-BFGS-B",
                            lower=c(rep(-Inf, dX),ranges[,1]), upper=c(rep(Inf, dX),ranges[,2])),
                      error=function(e) return(list(value=Inf)))
      if (opt$value < best$value) best <- opt
    }

  }
  # print(-best$value)
  # print(ranges)
  return(-best$value)

  # ### get maximum weights in each segment
  # for (i in m_wts) {
  #   rg
  #   opt <- optim(par)
  #
  # }
}

