suppressMessages(library(survey))
set.seed(123)
n <- 1e4

## check that plasmode approach works.  First simulate a dataset:
df <- data.frame(Z1=rnorm(n), Z2=factor(sample(3, size=n, replace=TRUE)))

forms <- list(list(),
              list(A ~ Z1 + Z2),
              list(Y ~ A),
              ~ A)
fams <- list(integer(0), 5, 1, 1)
pars <- list(A = list(beta=c(-0.5,0.5,0.25,-0.5)),
             Y = list(beta=c(-0.5,1), phi=1),
             cop = list(Y = list(Z1 = list(beta=c(0.5,0.5)),
                                 Z2 = list(beta=c(0.75,0.2)))))

dat <- rfrugalParam(formulas = forms, family = fams, pars = pars, dat = df)

test_that("plasmode simulation works", {
  expect_equal(dim(dat), c(n, 4L))
  expect_false(any(is.na(dat)))
})

## obtain parameter estimates
glmA <- glm(A ~ Z1 + Z2, data=dat, family = binomial)
ps <- predict(glmA, type="response")
wt <- dat$A/ps + (1-dat$A)/(1-ps)

glmY <- svyglm(Y ~ A, design=svydesign(~1, weights=wt, data=dat))

coefA <- summary(glmA)$coef
coefY <- summary(glmY)$coef

test_that("plasmode simulation generates correct data", {
  expect_lt(max(abs(coefA[,1] - pars$A$beta)/coefA[,2]), 2.5)
  expect_lt(max(abs(coefY[,1] - pars$Y$beta)/coefY[,2]), 2.5)
})


## check that works with 1 prespecified variable
forms <- list(list(),
              list(A ~ U),
              list(Y ~ A),
              list(~ 1))
fams <- list(integer(0), 5, 1, 1)
pars <- list(A = list(beta=c(0,1)),
             Y = list(beta=c(1, 1), phi=1),
             cop = list(beta=-1))

dat <- data.frame(U = rep(1,100))
dat <- rfrugalParam(formulas=forms, family=fams, pars=pars, dat=dat)
