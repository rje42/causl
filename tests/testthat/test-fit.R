fam <- rep(1,4)
pars <- list(z = list(beta=0, phi=1),
             x = list(beta=c(0,0.5), phi=1),
             y = list(beta=c(0,0.5), phi=0.5),
             cop = list(beta=matrix(1, 1, 1)))

set.seed(123)
dat <- causalSamp(1e2, par=pars, family=fam)
out <- fit_causl(dat, family=c(1,1,1))

set.seed(124)
fam2 <- c(1,5,1,1)
dat2 <- causalSamp(5e2, par=pars, family=fam2)
out2 <- fit_causl(dat2, family=c(1,1,1))


test_that("fitting works as expected", {
  testthat::expect_lt(sum((out$pars$y$beta - c(0.004957232, 0.548474548))^2), 1e-10)
  testthat::expect_lt(sum((out2$pars$y$beta - c(-0.01425083, 0.42604663))^2), 1e-10)
})

# nll1 <- nll3(dat, family=rep(1,3), inCop = c("y","z"))
# out3 <- optim(nll1, par=c(0,0,0,0,0,1,1))
#
# test_that("fitting3 works as expected", {
#   testthat::expect_lt(sum((out3$pars - c(-0.00216897817065516, 0.591080756662059, -0.105764362140315,
#                                          0.714666339742363, -0.10159280239843, 0.444723304635891,
#                                          0.763841648563355))^2), 1e-6)
# })

fam[1] <- 3

cm <- causl_model(formulas = list(z ~ 1, x ~ z, y ~ x, ~ 1),
                  family = fam,
                  pars = pars,
                  link = c("inverse", "identity", "identity"))
cm$pars$z$beta <- 2

set.seed(123)
dat <- suppressMessages(rfrugal(1e2, causl_model = cm))
out <- fit_causl(dat, list(z ~ 1, y ~ x, cop ~ 1), family=fam[-2], link = c("inverse", "identity"))

test_that("fitting works as expected with different links", {
  testthat::expect_lt(sum((out$pars$z$beta - c(2.21898240437888))^2), 1e-12)
})
