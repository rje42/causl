# initializeParams <- function(dat, formulas, family, init, LHS, wh) {
#
#   d <- length(formulas) - 1
#   # fam_y <- family[1]
#   # fam_z <- family[1+seq_len(d)]
#   fam_cop <- last(family)
#
#   ## initialization of parameters
#   beta_start = rep(0, last(unlist(last(wh))))
#
#   ## get parameters for y variable
#   # if (fam_y <= 2) {
#   #   if (init) {
#   #     form <- update.formula(forms[[1]], paste0(". ~ . + ", paste(LHS[-1], collapse=" + ")))
#   #     lm_y <- summary(lm(form, data=dat))
#   #     beta_start[wh[[1]]] <- lm_y$coefficients[wh[[1]],1]
#   #     beta_start[wh[[2]]] <- lm_y$sigma
#   #   }
#   #   else beta_start[wh[[2]]] = sd(dat[[LHS[1]]])
#   # }
#   # else if (fam_y == 3) {
#   #   glm_y <- glm(forms[[1]], family = Gamma(link = "log"),
#   #                data = dat)
#   #   beta_start[wh[[2]]] = summary(glm_y)$dispersion
#   # }
#   # else stop("fam_y must be 1, 2 or 3")
#
#   if (init) lm_z <- vector(mode="list", length=d)
#
#   ## get parameters for z variable(s)
#   for (i in seq_len(d)) {
#     if (family[i] <= 2 || family[i] == 5) {
#       if (init) {
#         lm_z[[i]] <- summary(lm(formulas[[i]], data=dat))
#         beta_start[wh[[2*i - 1]]] <- lm_z[[i]]$coefficients[,1]
#         beta_start[wh[[2*i]]] <- lm_z[[i]]$sigma
#       }
#       else beta_start[wh[[2*i]]] = sd(dat[[LHS[i]]])
#     }
#     else if (family[i] == 3) {
#       if (init) {
#         lm_z[[i]] <- glm(formulas[[i]], family = Gamma(link = "log"),
#                    data = dat)
#         beta_start[wh[[2*i - 1]]] = lm_z[[i]]$residuals
#         beta_start[wh[[2*i]]] = summary(lm_z[[i]])$dispersion
#       }
#       else beta_start[wh[[2*i]]] = sd(dat[[LHS[i]]])
#     }
#     else stop("fam_z must be 1, 2, 3 or 5")
#   }
#
#   if (init && fam_cop <= 2) {
#     resids <- lapply(lm_z, function(x) x$residuals)
#     tmp <- cor(do.call(data.frame, resids))
#     tmp <- tmp[upper.tri(tmp)]
#     for (i in seq_len(choose(d,2))) {
#       beta_start[wh[[2*d + i]][1]] <- logit((tmp[i]+1)/2)
#     }
#   }
#
#   return(beta_start)
# }

initializeParams2 <- function(dat, formulas, family=rep(1,nv), link, init=FALSE,
                              full_form, kwd, notInCop, inc_cop=TRUE, nc, only_masks=FALSE) {

  # d <- ncol(dat)
  # fam_y <- family[1]
  # fam_z <- family[1+seq_len(d)]
  # fam_cop <- last(family)
  if (missing(dat)) {
    if (!only_masks) stop("Need data frame to initialize parameters")
    else if (missing(nc)) stop("Need data frame to initialize parameters")
  }
  else nc <- ncol(dat)
  nv <- length(formulas) - (inc_cop)

  if (length(family) != nv+(inc_cop)) stop("Must have family parameter for every variable (and copula if included)")

  ## get terms labels
  trms <- terms(full_form$formula)
  labs <- if (attr(trms, "intercept") == 1) c("(intercept)")
  else character(0)
  labs <- c(labs, attr(trms, "term.labels"))

  wh <- full_form$wh

  ## define output matrix/vector for beta,phi
  LHS <- lhs(formulas)
  if (inc_cop) LHS <- LHS[-match(kwd, LHS, nomatch=0L)]

  if (inc_cop) {
    ## get column names for beta/beta_m
    if (missing(notInCop)) LHSc <- LHS
    else LHSc <- setdiff(LHS, notInCop)
    nc <- length(LHSc)
    c2 <- combn(nc, 2)
    colnms <- c(LHS, mapply(function(x,y) paste0(kwd,"_",LHSc[x],"_",LHSc[y]), c2[1,], c2[2,]))

    ## check which variables are already in the formula for another, and set
    ## corresponding copula correlation to zero
    wh_in_cop <- match(LHSc, LHS)
    idx <- lapply(full_form$old_forms[wh_in_cop],
                  function(x) match(intersect(LHSc, attr(x, "term.labels")[attr(x, "order") == 1]), LHSc))
    ignr <- do.call(cbind, mapply(function(x,y) if (length(x) > 0) rbind(x,y) else matrix(NA, 2, 0), idx, seq_along(LHSc)))
    if (length(ignr) > 0) {
      wh_swt <- ignr[1,] > ignr[2,]
      ignr[,wh_swt] <- ignr[2:1, wh_swt]
      in_form <- setmatch(apply(ignr, 2, c, simplify = FALSE),
                          apply(c2, 2, c, simplify = FALSE))
    }
    else in_form <- numeric(0)
  }
  else colnms <- LHS

  ## now set up blank masks
  # get_masks(formula, nv, nc)

  beta <- beta_m <- matrix(0, nrow=max(unlist(wh)), ncol=nv+choose(nc,2),
                           dimnames = list(labs,colnms))
  phi <- phi_m <- numeric(nv)

  # FIGURE OUT HOW TO MAKE THIS WORK WITH NOT EVERYTHING IN THE COPULA

  # LHSs <- lhs(formulas)

  # if (init) {
  #   wts <- ipw_weights(dat, formulas[-length(formulas)])
  # }
  # else wts <- rep(1,nrow(dat))
  # dat <- cbind(dat, wts=wts)

  ## intialize parameters for each variable
  for (i in seq_along(phi)) {
    beta_m[wh[[i]],i] <- 1

    if (!only_masks && init) {
      fam <- switch(family[i], "1"=gaussian, "2"=gaussian, "3"=Gamma,
                    stop("family should be Gaussian, t or Gamma"))
      # mod_fit <- survey::svyglm(formula=formulas[[i]], family=do.call(fam, list(link=link[[i]])), design=survey::svydesign(~1, data=dat,weights=ipw))
      mod_fit <- glm(formula=formulas[[i]], data=dat, family=do.call(fam, list(link=link[[i]])))
      beta[wh[[i]],i] <- c(mod_fit$coef[1], mod_fit$coef[-1]/2)
      if (family[i] < 5) {
        phi[i] <- summary(mod_fit)$dispersion
        phi_m[i] <- 1
      }
      next
    }
    if (only_masks) {
      phi_m[family[-length(family)] %in% 1:3] <- 1
    }
    else {
      ## pick mean and sd...
      if (family[i] <= 2) {
        if (link[i] == "identity") beta[1,i] <- mean(dat[[LHS[i]]])
        else if (link[i] == "inverse") beta[1,i] <- 1/mean(dat[[LHS[i]]])
        else if (link[i] == "log") beta[1,i] <- log(exp(mean(dat[[LHS[i]]])))

        phi[i] <- var(dat[[LHS[i]]])
        phi_m[i] <- 1
      }
      else if (family[i] == 3) {
        if (link[i] == "identity") beta[1,i] <- mean(dat[[LHS[i]]])
        else if (link[i] == "inverse") beta[1,i] <- 1/mean(dat[[LHS[i]]])
        else if (link[i] == "log") beta[1,i] <- log(exp(mean(dat[[LHS[i]]])))

        phi[i] <- var(dat[[LHS[i]]])
        phi_m[i] <- 1
      }
      else phi[i] <- NA
    }
  }

  if (inc_cop) {
    ## initialize copula parameters
    beta_m[wh[[nv+1]], nv+seq_len(choose(nc,2))] <- 1
    beta_m[wh[[nv+1]], nv+in_form] <- 0
    # beta[-wh[[i]],i] <- 0
  }

  if (only_masks) out <- list(beta_m=beta_m, phi_m=phi_m)
  else out <- list(beta=beta, phi=phi, beta_m=beta_m, phi_m=phi_m)

  return(out)
}

##' Get parameter masks for regression parameters
##'
##' @param formulas formulas to create mask for
##' @param family vector or list of families
##' @param full_form (optionally) merged list of `formulas`
##
par_masks <- function(formulas, family=rep(1,nv), full_form) { #, kwd, only_masks=FALSE) {

  # d <- ncol(dat)
  # fam_y <- family[1]
  # fam_z <- family[1+seq_len(d)]
  # fam_cop <- last(family)
  nv <- length(formulas)
  fam <- family_list(family, func_return = get_family)

  if (length(fam) != nv) stop("Must have family parameter for every variable")

  ## get terms and outcome labels
  if (missing(full_form)) full_form <- merge_formulas(formulas)
  trms <- terms(full_form$formula)
  labs <- if (attr(trms, "intercept") == 1) c("(intercept)")
  else character(0)
  labs <- c(labs, attr(trms, "term.labels"))
  LHSs <- lhs(formulas)

  wh <- full_form$wh
  beta_m <- matrix(0, nrow=max(unlist(wh)), ncol=nv, dimnames = list(labs, LHSs))
  phi_m <- rep(0, nv)
  names(phi_m) <- LHSs

  ## define output matrix/vector for beta,phi
  LHS <- lhs(formulas)

  ## add in indicators for each variable
  for (i in seq_len(nv)) {
    beta_m[wh[[i]],i] <- 1
    phi_m[i] <- match("phi", names(formals(fam[[i]]()$ddist)), nomatch = 0L) > 0
  }

  out <- list(beta_m=beta_m, phi_m=phi_m)

  return(out)
}


masks <- function(formulas, family=rep(1,nc+1), wh, LHS, cp) {

  if (is.list(family)) family <- unlist(family[1:3])
  formulas <- unlist(formulas)
  nc <- length(formulas)-1

  phi_m <- numeric(nc)
  if (missing(cp)) cp <- length(phi_m)
  beta_m <- matrix(0, nrow=max(unlist(wh)), ncol=nc+choose(cp,2))

  # LHSs <- lhs(formulas)

  ## intialize parameters for each variable
  for (i in seq_along(phi_m)) {
    if (family[i] >= 1 && family[i] <= 3) {
      phi_m[i] <- 1
    }

    beta_m[wh[[i]],i] <- 1
  }

  ## initialize copula parameters
  beta_m[wh[[nc+1]], nc+seq_len(choose(cp,2))] <- 1

  return(list(beta_m=beta_m, phi_m=phi_m))
}

pars2mask <- function(pars, masks) {
  out <- masks

  wh <- pmatch(c("beta", "phi"), names(masks))

  if (!is.na(wh[1])) {
    names(out)[wh[1]] <- "beta"
    for (i in seq_along(pars)) {
      out$beta[masks[[wh[1]]][,i] > 0,i] <- pars[[i]]$beta
    }
  }
  if (!is.na(wh[2])) {
    names(out)[wh[2]] <- "phi"
    for (i in seq_along(out$phi)) {
      if (masks[[wh[2]]][i] > 0) out$phi[i] <- pars[[i]]$phi
    }
  }

  return(out)
}

## Get theta vector from pars
get_theta <- function (pars, formulas, full_form, kwd="cop") {

  ## get terms labels
  # trms <- terms(full_form$formula)
  # labs <- if (attr(trms, "intercept") == 1) c("(intercept)")
  # else character(0)
  # labs <- c(labs, attr(trms, "term.labels"))

  nv <- length(formulas) - 1
  wh <- full_form$wh

  ## define output matrix/vector for beta,phi
  LHS <- lhs(formulas)
  # LHS <- LHS[-match(kwd, LHS)]
  # if (length(LHS) != nv) stop("Formula for copula causing problems!")

  wh_par <- match(LHS, names(pars))

  # c2 <- combn(length(LHS), 2)
  # colnms <- c(LHS, mapply(function(x,y) paste0(kwd,"_",LHS[x],"_",LHS[y]), c2[1,],c2[2,]))

  beta <- beta_m <- matrix(0, nrow=max(unlist(wh)), ncol=nv+choose(nv,2))
  phi <- phi_m <- numeric(nv)

  # FIGURE OUT HOW TO MAKE THIS WORK WITH NOT EVERYTHING IN THE COPULA

  # LHSs <- lhs(formulas)

  ## intialize parameters for each variable
  for (i in seq_along(phi)) {
    beta_m[wh[[i]],i] <- 1
    beta[wh[[i]],i] <- pars[[wh_par[i]]]$beta

    if (!is.null(pars[[wh_par[i]]]$phi)) {
      phi_m[i] <- 1
      phi[i] <- pars[[wh_par[i]]]$phi
    }
  }

  ## get copula parameters
  beta[wh[[length(phi)+1]], nv+seq_len(choose(nv,2))] <- pars[[kwd]]$beta
  beta_m[wh[[length(phi)+1]], nv+seq_len(choose(nv,2))] <- 1

  theta <- c(beta[beta_m > 0], phi[phi_m > 0])

  return(theta)
}

ipw_weights <- function(dat, formulas, min_wt=0.1) {

  if (min_wt >= 1) stop("Minimum weight must be less than 1")

  ## get variables to inverse weight by
  covs <- lapply(formulas, rhs_vars)
  all_covs <- unique.default(unlist(covs))
  LHS <- lhs(formulas)
  out <- rep(1,nrow(dat))

  for (i in seq_along(all_covs)) {
    # now, for each variable that appears, reweight by its dependence on
    # variables whose regression it is _not_ in
    wh <- sapply(covs, function (x) all_covs[i] %in% x)
    if (all(wh)) next

    form <- as.formula(paste0(all_covs[i], " ~ ", paste(LHS[!wh], collapse=" + ")))
    mod <- lm(form, data=dat)
    wt <- dnorm(dat[[all_covs[i]]], cbind(1, as.matrix(dat[LHS[!wh]])) %*% mod$coefficients, sd=sd(mod$residuals))
    out <- out*wt
  }

  ## now normalize weights and enforce the minimum value
  out <- out/max(out)
  out[out < min_wt] <- min_wt

  out
}
