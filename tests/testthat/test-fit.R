fam <- rep(1,4)
pars <- list(z = list(beta=0, phi=1),
             x = list(beta=c(0,0.5), phi=1),
             y = list(beta=c(0,0.5), phi=0.5),
             cop = list(beta=matrix(1, 1, 1)))

pars2 <- list(z1 = list(beta=0, phi=1),
              z2 = list(beta=c(0,0.25), phi=1),
             x = list(beta=c(0,0.5), phi=1),
             y = list(beta=c(0,0.5), phi=0.5),
             cop = list(beta=matrix(c(0,0.5,0.5), 1, 3)))
forms2 <- list(list(z1~1, z2~z1), x~z1, y~x, ~1)

set.seed(123)
dat <- causalSamp(1e2, par=pars, family=fam)
out <- fitCausal(dat, family=c(1,1,1))

set.seed(124)
fam2 <- c(1,5,1,1)
dat2 <- causalSamp(5e2, par=pars, family=fam2)
out2 <- fitCausal(dat2, family=c(1,1,1))

set.seed(125)
fam3 <- list(c(1,1),5,1,1)
dat3 <- causalSamp(1e2, formulas=forms2, par=pars2, family=fam3)
out3 <- fitCausal(dat3, formulas=unlist(forms2[-2]), family=c(1,1,1,1))


test_that("fitting works as expected", {
  testthat::expect_lt(sum((out$pars$y$beta - c(0.004957232, 0.548474548))^2), 1e-10)
  testthat::expect_lt(sum((out2$pars$y$beta - c(-0.01425083, 0.42604663))^2), 1e-10)
  testthat::expect_lt(sum((out3$pars$y$beta - c(-0.018912648, 0.60378023))^2), 1e-10)
})

nll1 <- nll3(dat, family=rep(1,3), inCop = c("y","z"))
out3 <- optim(nll1, par=c(0,0,0,0,0,1,1))

test_that("fitting3 works as expected", {
  testthat::expect_lt(sum((out3$pars - c(-0.00216897817065516, 0.591080756662059, -0.105764362140315,
                                         0.714666339742363, -0.10159280239843, 0.444723304635891,
                                         0.763841648563355))^2), 1e-6)
})
