suppressMessages(library(survey, quietly = TRUE))

forms <- list(list(Z1 ~ 1, Z2 ~ 1),
              X ~ Z1 + Z2,
              Y ~ X,
              ~ 1)
fams <- list(c(1,1), 5, 1, 1)
pars <- list(Z1 = list(beta=0, phi=1),
             Z2 = list(beta=0, phi=1),
             X = list(beta=c(0,1,1)),
             Y = list(beta=c(1,2), phi=1),
             cop = list(beta=matrix(c(0.1,0.5,0.6), nrow=1)))

n <- 1e4
dat <- rfrugalParam(n, formulas = forms, family = fams, pars = pars,
                    method="inversion_mv", control=list(quiet=TRUE))

modX <- glm(X ~ Z1 + Z2, family=binomial, data=dat)
coefsX <- summary(modX)$coefficients
ps <- predict(modX, type="response")

wt <- dat$X/ps + (1-dat$X)/(1-ps)


modY <- svyglm(Y ~ X, design = svydesign(~ 1, weights=wt, data=dat))
coefsY <- summary(modY)$coefficients

test_that("multivariate copula method works", {
  expect_lt(max(abs((coefsX[,1] - c(0,1,1))/coefsX[,2])), 2.5)
  expect_lt(max(abs((coefsY[,1] - c(1,2))/coefsY[,2])), 2.5)
})

