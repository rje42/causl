## check that dGDCop2_sig works in the way that is intended
forms <- list(Z ~ 1, X ~ Z, Y ~ X, ~ X)
pars <- list(Z = list(beta = 0, phi = 1),
             X = list(beta = c(0,0.5)),
             Y = list(beta = c(-0.25,0.5)),
             cop = list(beta=c(1,0.25)))

set.seed(123)
dat <- suppressMessages(rfrugalParam(1e4, formulas = forms, family = c(1,5,5,1), pars = pars,
                    method = "inversion"))

modX <- glm(X ~ Z, data=dat, family=binomial)
ps <- predict(modX, type = "response")
wt <- dat$X/ps + (1-dat$X)/(1-ps)
modY <- suppressWarnings(glm(Y ~ X, data=dat, family=binomial, weights = wt))

test_that("distribution is correct", {
  expect_true(sqrt(sum(modX$coef - c(0,0.5))^2) < 0.01)
  expect_true(sqrt(sum(modY$coef - c(-0.25,0.5))^2) < 0.05)
})

## try fitting model by ML
out <- fit_causl(dat[1:1e3,], formulas = list(Y ~ X, Z ~ 1, ~ X), family=c(5,1,1))

test_that("fitting works", {
  expect_true(sqrt(sum(out$pars$Y$beta - c(-0.25,0.5))^2) < 0.1)
  expect_true(abs(out$pars$Z$beta - 0) < 0.05)
  expect_true(abs(out$pars$Z$phi - 1) < 0.05)
  expect_true(sqrt(sum(out$pars$cop$beta - c(1,0.25))^2) < 0.1)
})
