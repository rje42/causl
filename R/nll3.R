nll3 <- function(dat, formulas=list(y~x, z~1, ~x),
                 family=rep(1,length(formulas)), link, par2,
                 inCop, kwd="cop") {

  ## may need to fix this to allow more flexibility in copula
  d <- length(formulas) - 1
  if (length(family) != length(formulas)) {
    stop("Should be a family variable for each formula")
  }

  ## tidy up the formulae
  forms <- tidy_formulas(formulas, kwd=kwd)
  fam_cop <- last(family)
  family <- family[-length(family)]
  link <- link_setup(link, family = family)

  LHS <- lhs(forms[-length(forms)])
  if (missing(inCop)) inCop <- match(LHS, names(dat))
  if (any(is.na(inCop))) stop("Names from formulae not matched")

  ## now merge formulae
  full_form <- merge_formulas(forms)
  wh <- full_form$wh
  mm <- model.matrix(full_form$formula, data=dat)
  # mms = lapply(forms, model.matrix, data=dat)
  dat <- dat[,inCop]

  ## set secondary parameter to 4 if in a t-Copula model
  if (missing(par2)) {
    if (fam_cop == 2) {
      par2 <- 4
      message("par2 set to 4\n")
    }
    else par2 = 0
  }

  ## get masks and initial parameter values
  beta_start <- initializeParams2(dat, formulas=forms, family=c(family, fam_cop),
                                  link=link, full_form=full_form, kwd=kwd)
  beta <- beta_start$beta_m
  phi <- beta_start$phi_m

  np <- sum(beta > 0)
  nv <- length(phi)
  nc <- ncol(dat)

  ncv <- length(inCop)
  par <- vector(mode="list", length=choose(ncv,2))

  nll <- function(theta) {
    beta[beta > 0] <- theta[seq_len(np)]
    phi[phi > 0] <- theta[-seq_len(np)]

    if (any(phi < 0)) return(-Inf)

    ## compute etas for each variable
    eta <- mm %*% beta

    ## get the densities for each observation
    log_den <- dat_u <- matrix(NA, nrow(dat), ncv)
    # if (length(family) < nc) family <- rep_len(family, nc)

    ## get univariate densities
    for (i in which(family != 5)) {
      tmp <- univarDens(dat[,i], eta[,i], phi=phi[i], family=family[i])
      log_den[,i] <- tmp$ld
      dat_u[,i] <- tmp$u
    }
    # wh_trunc = 0
    ## deal with discrete variables separately
    for (i in which(family == 5)) {
      # wh_trunc <- wh_trunc + 1
      tmp <- univarDens(dat[,i], eta[,i], family=family[i])
      log_den[,i] <- tmp$ld
      dat_u[,i] <- tmp$u
    }

    ## now set up copula parameters
    for (i in seq_len(choose(ncv,2))) {
      if (fam_cop <= 2 || fam_cop == 11) {
        par[[i]] <- pmin(pmax(2*expit(eta[,i+nv])-1, -1+1e-10), 1-1e-10)
      }
      else if (fam_cop == 3) {
        par[[i]] <- exp(eta[,i+nv])
      }
      else if (fam_cop == 4 || fam_cop == 6) {
        par[[i]] <- exp(eta[,i+nv])+1
      }
      else if (fam_cop == 5) {
        par[[i]] <- eta[,i+nv]
      }
    }

    ## define appropriate matrices
    if (fam_cop != 11) {
      if (ncv > 2 || fam_cop <= 2) {
        Sigma <- rep(diag(ncv), length(par[[1]]))
        dim(Sigma) <- c(ncv,ncv,length(par[[1]]))
        k = 1
        for (j in seq_len(ncv)[-1]) for (i in seq_len(j-1)) {
          Sigma[i,j,] <- Sigma[j,i,] <- par[[k]]
          k <- k+1
        }

        ## deal with Gaussian and t-copulas
        if (fam_cop == 1) {
          # if (any(family == 5)) {
          #   new_ord <- order(family[inCop])
          #   dat_u2 <- dat_u[,inCop,drop=FALSE][,new_ord,drop=FALSE]
          #   Sigma <- Sigma[new_ord,new_ord,,drop=FALSE]
          #   cop <- dGaussDiscCop(dat_u2, Sigma=Sigma, trunc=mms$trunc, log=TRUE, useC=useC)
          # }
          # else cop <- dGaussCop(dat_u[,inCop,drop=FALSE], Sigma=Sigma, log=TRUE, useC=useC)

          cop <- dGaussCop(dat_u, Sigma=Sigma, log=TRUE, useC=TRUE)
          # cop <- dGaussCop(dat_u[,inCop,drop=FALSE], Sigma=Sigma, log=TRUE, useC=useC)
        }
        else if (fam_cop == 2) {
          #q_dat <- qt(as.matrix(dat_u), df=par2)
          # cop <- dtCop(dat_u[,inCop,drop=FALSE], Sigma=Sigma, df=par2, log=TRUE)
          cop <- dtCop(dat_u, Sigma=Sigma, df=par2, log=TRUE)
        }
        else stop("Only Gaussian and t-copulas implemented for more than two dimensions")
      }
      else if (nc == 2) {
        ### MODIFY THIS TO USE copula PACKAGE!
        cop <- log(VineCopula::BiCopPDF(dat_u[,1], dat_u[,2], family=fam_cop, par=par[[1]], par2=par2))
      }
      else stop("should have that nc is an integer >= 2")
    }
    else {
      cop <- log(causl::dfgmCopula(dat_u[,1], dat_u[,2], alpha=par[[1]]))
    }

    out <- cop + rowSums(log_den)

    -sum(out)
  }
}
