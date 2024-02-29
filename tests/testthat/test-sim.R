fam <- rep(1,4)
pars <- list(z = list(beta=0, phi=1),
             x = list(beta=c(0,0.5), phi=1),
             y = list(beta=c(0,0.5), phi=0.5),
             cop = list(beta=matrix(1, 1, 1)))

set.seed(123)
dat <- causalSamp(1e3, par=pars)

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
dat <- rfrugalParam(1e4, formulas = forms, family=fam, par=pars,
                    method="inversion")

z_p <- ks.test(pnorm(dat$z), runif(1e4))$p.value
x_p <- glm(x ~ z, family=binomial, data=dat)$coef

test_that("simulation (inv) works 1", {
  expect_gt(z_p, 0.05)
  expect_lt(abs(x_p[2]-0.5), 0.02)
})
