suppressMessages(library(survey, quietly = TRUE))

set.seed(130)
forms <- list(list(Z1 ~ 1, Z2 ~ Z1), list(X1 ~ Z1 + Z2, X2 ~ Z1*Z2 + X1),
              list(Y1 ~ X1*X2, Y2 ~ X1*X2), ~ 1)
fam <- list(c(1,1), c(5,5), c(3,1), 1)
fam2 <- list(c("gaussian","gaussian"), c("binomial","binomial"), c("Gamma","gaussian"), 1)

pars <- list(Z1 = list(beta=0, phi=1),
             Z2 = list(beta=c(0,0.5), phi=0.25),
             X1 = list(beta=c(-0.2,0.1,0.4)),
             X2 = list(beta=c(0,0.1,0.4,0.4,-0.05)),
             Y1 = list(beta=c(-0.4,0.3,0.4,-0.2), phi=1),
             Y2 = list(beta=c(-0.5,0.5,0.2,-0.2), phi=0.75),
             cop = list(beta=matrix(c(0,rep(1,5)), nrow=1)))
pars2 <- pars
pars2$cop <- list(beta=1)

n <- 1e5
dat0 <- suppressMessages(rfrugalParam(n, formulas = forms, family = fam, pars = pars,
                    method = "rejection"))
dat2 <- suppressMessages(rfrugalParam(n, formulas = forms, family = fam2, pars = pars2,
                     method = "inversion"))
dats <- list(dat0, dat2)

for (dat in dats) {
  lmZ1 <- summary(lm(Z1 ~ 1, data=dat))
  lmZ2 <- summary(lm(Z2 ~ Z1, data=dat))
  lmX1 <- glm(X1 ~ Z1 + Z2, family=binomial, data=dat)
  slmX1 <- summary(lmX1)
  lmX2 <- glm(X2 ~ Z1*Z2 + X1, family=binomial, data=dat)
  slmX2 <- summary(lmX2)
  zZ1 <- (lmZ1$coef[,1] - pars$Z1$beta)/lmZ1$coef[,2]
  zZ2 <- (lmZ2$coef[,1] - pars$Z2$beta)/lmZ2$coef[,2]
  zX1 <- (slmX1$coef[,1] - pars$X1$beta)/slmX1$coef[,2]
  zX2 <- (slmX2$coef[,1] - pars$X2$beta)/slmX2$coef[,2]
  ps1 <- predict(lmX1, type = "response")
  ps2 <- predict(lmX2, type = "response")
  wts <- (dat$X1/ps1 + (1-dat$X1)/(1-ps1))*(dat$X2/ps2 + (1-dat$X2)/(1-ps2))

  lmY1 <- summary(svyglm(Y1 ~ X1*X2, family=Gamma(link=log),
                         design = svydesign(~1, weights = ~wts, data=dat)))
  lmY2 <- summary(svyglm(Y2 ~ X1*X2, design = svydesign(~1, weights = ~wts,
                                                        data=dat)))
  zY1 <- (lmY1$coef[,1] - pars$Y1$beta)/lmY1$coef[,2]
  zY2 <- (lmY2$coef[,1] - pars$Y2$beta)/lmY2$coef[,2]

  test_that("simulation with multiple variables works", {
    expect_lt(max(abs(zZ1)), 3)
    expect_lt(max(abs(zZ2)), 3)
    expect_lt(max(abs(zX1)), 3)
    expect_lt(max(abs(zX2)), 3)
    expect_lt(max(abs(zY1)), 3)
    expect_lt(max(abs(zY2)), 3)

    expect_gt(ks.test(pnorm(c(zZ1,zZ2,zX1,zX2,zY1,zY2)), (1:18)/19)$p.value, 0.05)
  })
}

