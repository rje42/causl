library(survey)

forms <- list(list(L0 ~ A0, L1 ~ A0),
              list(A0 ~ 1, A1 ~ A0*L0),
              Y ~ A0*A1,
              ~ A0)
pars <- list(A0 = list(beta = 0),
             L0 = list(beta = c(0.3,-0.2,0.4,0.1), nlevel=3),
             L1 = list(beta = c(0.3,-0.2,0.4,0.1), nlevel=3),
             A1 = list(beta = 2*c(-0.3,0.4,0.3,0.5,0.3,0.5)),
             Y = list(beta = c(-0.5,0.2,0.3,0), phi=1),
             cop = list(beta=c(2,0.5)))

set.seed(123)
n <- 5e4
dat <- rfrugalParam(n, formulas = forms, pars=pars,
                    family = list(c("categorical","categorical"),c(5,5),1,1),
                    method="inversion")

ps <- predict(glm(A1 ~ A0*L0, family=binomial, data=dat), type="response")
wt <- dat$A1/ps + (1-dat$A1)/(1-ps)
hist(wt, breaks=25)

cor_tab <- summary(svyglm(Y ~ A0*A1, design = svydesign(~1,data=dat,weights=wt)))$coef
wrg_tab <- summary(svyglm(Y ~ A0*A1, design = svydesign(~1,data=dat,weights=~1)))$coef

testthat::test_that("categorical variables simulated correctly", {
  expect_true(max(abs(cor_tab[,1]-c(-0.5,0.2,0.3,0))/cor_tab[,2]) < 2.5)
  expect_true(max(abs(wrg_tab[,1]-c(-0.5,0.2,0.3,0))/wrg_tab[,2]) > 10)
})
