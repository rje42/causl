fam <- rep(1,4)
pars <- list(z = list(beta=0, phi=1),
             x = list(beta=c(0,0.5), phi=1),
             y = list(beta=c(0,0.5), phi=0.5),
             cop = list(beta=matrix(1, 1, 1)))

# test_that("No families message", {
#   expect_message(dat <- causalSamp(10, par=pars),
#                  "No families provided, so assuming all variables are Gaussian")
# })

set.seed(123)
dat <- suppressMessages(causalSamp(1e3, family=c(1,1,1,1), par=pars))
z_p <- ks.test(dat$z, pnorm)$p.value
x_p <- ks.test(dat$x - 0.5*dat$z, pnorm)$p.value
# mvn_p <- MVN::mvn(dat, mvnTest = "hz")$multivariateNormality$`p value`

test_that("simulation works 1", {
  expect_gt(z_p, 0.02)
  expect_gt(x_p, 0.02)
  # expect_gt(mvn_p, 0.02)
})

fam <- list(c(1,3),5,3,c(1,1,1))
forms <- list(list(z ~ 1, u ~ 1), list(x ~ z), list(y ~ x), ~ x)
pars <- list(z = list(beta=0, phi=1),
             u = list(beta=0, phi=1),
             x = list(beta=c(0,0.5)),
             y = list(beta=c(0,0.5), phi=0.5),
             cop = list(beta=matrix(c(0,0.5,
                                      -0.5,0.2,
                                      0,0.5), nrow=2)))

set.seed(124)
dat <- causalSamp(1e4, formulas = forms, family=fam, par=pars)

z_p <- ks.test(dat$z, pnorm)$p.value
x_p <- glm(x ~ z, family=binomial, data=dat)$coef

test_that("simulation works 2", {
  expect_gt(z_p, 0.02)
  expect_lt(abs(x_p[2]-0.5), 0.1)
})

## check that binomial works OK
fam <- c(5,1,1,1)
forms <- list(z ~ 1, x ~ z, y ~ x, ~ 1)
pars <- list(z = list(beta=0),
             x = list(beta=c(0,0.5), phi=1),
             y = list(beta=c(0,0.5), phi=0.5),
             cop = list(beta=matrix(1, 1, 1)))

set.seed(125)
dat <- causalSamp(1e4, formulas = forms, family=fam, par=pars)

z_p <- chisq.test(table(dat$z), p=c(0.5,0.5))$p.value
x_p <- lm(x ~ z, data=dat)$coef

test_that("simulation works 3a", {
  expect_gt(z_p, 0.02)
  expect_lt(abs(x_p[2]-0.5), 0.01)
})

fam <- c(5,1,1,1)

set.seed(125)
dat <- causalSamp(1e4, formulas = forms, family=fam, par=pars)

z_p <- chisq.test(table(dat$z), p=c(0.5,0.5))$p.value
x_p <- lm(x ~ z, data=dat)$coef

test_that("simulation works 3b", {
  expect_gt(z_p, 0.02)
  expect_lt(abs(x_p[2]-0.5), 0.01)
})

set.seed(124)
fam <- c(1,5,1,1)
pars <- list(z = list(beta=0, phi=1),
             x = list(beta=c(0,0.5)),
             y = list(beta=c(0,0.5), phi=0.5),
             cop = list(y=list(z=list(beta=1))))
cm <- causl_model(formulas = forms, family=fam, par=pars, method="inversion")
dat <- rfrugal(n=1e4, causl_model = cm, control=list(quiet=TRUE))
# dat <- rfrugalParam(1e4, formulas = forms, family=fam, par=pars,
                    # method="inversion", control=list(quiet=TRUE))

z_p <- ks.test(pnorm(dat$z), runif(1e4))$p.value
x_p <- glm(x ~ z, family=binomial, data=dat)$coef

test_that("simulation (inv) works 1", {
  expect_gt(z_p, 0.05)
  expect_lt(abs(x_p[2]-0.5), 0.04)
})

## test log link for binomial
# set.seed(124)
# fam <- c(1,5,1,1)
# pars <- list(z = list(beta=0, phi=1),
#              x = list(beta=c(-5,0.5)),
#              y = list(beta=c(0,0.5), phi=0.5),
#              cop = list(y=list(z=list(beta=1))))
# dat <- suppressMessages(rfrugalParam(1e4, formulas = forms, family=fam, par=pars,
#                                      link = list("identity", "log", "identity")))

set.seed(124)
fam <- list(1,5,5,1)
forms <- list(X ~ 1,
              list(A ~ X),
              Y ~ A,
              ~ A)

pars <- list(X = list(beta = 0,phi = 1),
             A = list(beta = c(-1,2)),
             Y = list(beta = c(-1,0.2)),
             cop = list(beta = c(0,0.5)))

set.seed(123)
dat <- suppressMessages(rfrugalParam(1e5, formulas = forms, pars=pars, family = fam,
                    link = list("identity","logit","log")))


modX <- glm(A ~ X, family=binomial, data=dat)
ps <- fitted(modX)
wt <- dat$A/ps + (1-dat$A)/(1-ps)
modY <- suppressWarnings(glm(Y ~ A, data=dat, family=binomial("log"), weights = wt))


test_that("distribution is correct", {
  expect_true(sqrt(sum(modX$coef - c(-1,2))^2) < 0.05)
  expect_true(sqrt(sum(modY$coef - c(-1,0.2))^2) < 0.05)
})


cm2 <- modify.causl_model(cm, pars=list(x=list(beta=c(-5,0.5))),
                          link = list("identity", "log", "identity"))
set.seed(124)
dat <- rfrugal(n=1e5, causl_model=cm2, control=list(quiet=TRUE))

modx <- summary(glm(x ~ z, family=binomial(log), data=dat))
ztest <- abs(modx$coefficients[,1] - c(-5,0.5))/modx$coefficients[,2]

test_that("log-link works", {
  expect_lt(ztest[1], 2.5)
  expect_lt(ztest[2], 2.5)
})

