## check that dGDCop2_sig works in the way that is intended
forms <- list(Z ~ 1, X ~ Z, Y ~ X, ~ X)
pars <- list(Z = list(beta = 0, phi = 1),
             X = list(beta = c(0,0.5)),
             Y = list(beta = c(-0.25,0.5)),
             cop = list(beta=c(1,0.25)))

set.seed(123)
dat <- rfrugalParam(1e4, formulas = forms, family = c(1,5,5,1), pars = pars,
                    method = "inversion", control=list(quiet=TRUE))

modX <- glm(X ~ Z, data=dat, family=binomial)
ps <- predict(modX, type = "response")
wt <- dat$X/ps + (1-dat$X)/(1-ps)
modY <- suppressWarnings(svyglm(Y ~ X, family=binomial, design=svydesign(~1, data=dat, weights = wt)))
modYn <- suppressWarnings(svyglm(Y ~ X, family=binomial, design=svydesign(~1, data=dat, weights = ~1)))
smodX <- summary(modX)
smodY <- summary(modY)
smodYn <- summary(modYn)

Xerr <- sum((smodX$coef[,1] - pars$X$beta)^2/smodX$coef[,2]^2)
Yerr <- sum((smodY$coef[,1] - pars$Y$beta)^2/smodY$coef[,2]^2)
Ynerr <- sum((smodYn$coef[,1] - pars$Y$beta)^2/smodYn$coef[,2]^2)

test_that("distribution is correct", {
  expect_lt(Xerr, qchisq(0.99, df=2))
  expect_lt(Yerr, qchisq(0.99, df=2))
  expect_gt(Ynerr, qchisq(0.99, df=2))
})

## try fitting model by ML
out <- fit_causl(dat[1:1e3,], formulas = list(Y ~ X, Z ~ 1, ~ X), family=c(5,1,1))

Yerr <- sum((out$pars$Y$beta - pars$Y$beta)^2/out$pars$cop$beta_sandwich^2)
coperr <- sum((out$pars$cop$beta - pars$cop$beta)^2/out$pars$cop$beta_sandwich^2)

test_that("fitting works", {
  expect_lt(Yerr, qchisq(0.99, df=2))
  expect_lt(abs(out$pars$Z$beta - 0), 2.5*out$pars$Z$beta_sandwich)
  expect_lt(abs(out$pars$Z$phi - 1), 2.5*out$pars$Z$phi_sandwich)
  expect_lt(coperr, qchisq(0.99, df=2))
})

