initializeParams <- function(dat, formulas, family, init, LHS, wh) {

  d <- length(formulas) - 1
  # fam_y <- family[1]
  # fam_z <- family[1+seq_len(d)]
  fam_cop <- last(family)

  ## initialization of parameters
  beta_start = rep(0, last(unlist(last(wh))))

  ## get parameters for y variable
  # if (fam_y <= 2) {
  #   if (init) {
  #     form <- update.formula(forms[[1]], paste0(". ~ . + ", paste(LHS[-1], collapse=" + ")))
  #     lm_y <- summary(lm(form, data=dat))
  #     beta_start[wh[[1]]] <- lm_y$coefficients[wh[[1]],1]
  #     beta_start[wh[[2]]] <- lm_y$sigma
  #   }
  #   else beta_start[wh[[2]]] = sd(dat[[LHS[1]]])
  # }
  # else if (fam_y == 3) {
  #   glm_y <- glm(forms[[1]], family = Gamma(link = "log"),
  #                data = dat)
  #   beta_start[wh[[2]]] = summary(glm_y)$dispersion
  # }
  # else stop("fam_y must be 1, 2 or 3")

  if (init) lm_z <- vector(mode="list", length=d)

  ## get parameters for z variable(s)
  for (i in seq_len(d)) {
    if (family[i] <= 2 || family[i] == 5) {
      if (init) {
        lm_z[[i]] <- summary(lm(formulas[[i]], data=dat))
        beta_start[wh[[2*i - 1]]] <- lm_z[[i]]$coefficients[,1]
        beta_start[wh[[2*i]]] <- lm_z[[i]]$sigma
      }
      else beta_start[wh[[2*i]]] = sd(dat[[LHS[i]]])
    }
    else if (family[i] == 3) {
      if (init) {
        lm_z[[i]] <- glm(formulas[[i]], family = Gamma(link = "log"),
                   data = dat)
        beta_start[wh[[2*i - 1]]] = lm_z[[i]]$residuals
        beta_start[wh[[2*i]]] = summary(lm_z[[i]])$dispersion
      }
      else beta_start[wh[[2*i]]] = sd(dat[[LHS[i]]])
    }
    else stop("fam_z must be 1, 2, 3 or 5")
  }

  if (init && fam_cop <= 2) {
    resids <- lapply(lm_z, function(x) x$residuals)
    tmp <- cor(do.call(data.frame, resids))
    tmp <- tmp[upper.tri(tmp)]
    for (i in seq_len(choose(d,2))) {
      beta_start[wh[[2*d + i]][1]] <- logit((tmp[i]+1)/2)
    }
  }

  return(beta_start)
}

initializeParams2 <- function(dat, formulas, family=rep(1,nc), init, LHS, wh) {

  # d <- ncol(dat)
  # fam_y <- family[1]
  # fam_z <- family[1+seq_len(d)]
  # fam_cop <- last(family)
  nc <- ncol(dat)

  beta <- beta_m <- matrix(0, nrow=max(unlist(wh)), ncol=length(family))
  phi <- phi_m <- numeric(length(LHS))

# FIGURE OUT HOW TO MAKE THIS WORK WITH NOT EVERYTHING IN THE COPULA

  # LHSs <- lhs(formulas)

  ## intialize parameters for each variable
  for (i in seq_along(phi)) {
    if (family[i] <= 2) {
      beta[1,i] <- mean(dat[[LHS[i]]])
      phi[i] <- var(dat[[LHS[i]]])
      phi_m[i] <- 1
    }
    else if (family[i] == 3) {
      beta[1,i] <- exp(mean(dat[[LHS[i]]]))
      phi[i] <- var(dat[[LHS[i]]])
      phi_m[i] <- 1
    }
    else phi[i] <- NA

    beta_m[wh[[i]],i] <- 1
  }


  ## initialize copula parameters
  cp <- length(phi) + 1
  beta_m[wh[[cp]], cp] <- 1
  # beta[-wh[[i]],i] <- 0

  # ## get parameters for z variable(s)
  # for (i in seq_len(d)) {
  #   if (family[i] <= 2 || family[i] == 5) {
  #     if (init) {
  #       lm_z[[i]] <- summary(lm(formulas[[i]], data=dat))
  #       beta_start[wh[[2*i - 1]]] <- lm_z[[i]]$coefficients[,1]
  #       beta_start[wh[[2*i]]] <- lm_z[[i]]$sigma
  #     }
  #     else beta_start[wh[[2*i]]] = sd(dat[[LHS[i]]])
  #   }
  #   else if (family[i] == 3) {
  #     if (init) {
  #       lm_z[[i]] <- glm(formulas[[i]], family = Gamma(link = "log"),
  #                        data = dat)
  #       beta_start[wh[[2*i - 1]]] = lm_z[[i]]$residuals
  #       beta_start[wh[[2*i]]] = summary(lm_z[[i]])$dispersion
  #     }
  #     else beta_start[wh[[2*i]]] = sd(dat[[LHS[i]]])
  #   }
  #   else stop("fam_z must be 1, 2, 3 or 5")
  # }

  # if (init && fam_cop <= 2) {
  #   resids <- lapply(lm_z, function(x) x$residuals)
  #   tmp <- cor(do.call(data.frame, resids))
  #   tmp <- tmp[upper.tri(tmp)]
  #   for (i in seq_len(choose(d,2))) {
  #     beta_start[wh[[2*d + i]][1]] <- logit((tmp[i]+1)/2)
  #   }
  # }

  return(list(beta=beta, phi=phi, beta_m=beta_m, phi_m=phi_m))
}
