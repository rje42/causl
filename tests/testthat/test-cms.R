forms <- list(list(Z1 ~ 1, Z2 ~ Z1),
              X ~ Z1 + Z2,
              Y ~ X,
              ~ 1)
fam <- list(c(1,3), 5, 1, 1)
pars <- list(Z1 = list(beta=0, phi=1),
             Z2 = list(beta=c(-0.2,0.5), phi=1),
             X = list(beta=c(-0.5,0.2,0.4)),
             Y = list(beta=c(-0.5,0.75), phi=1),
             cop = list(beta=matrix(c(0,0.3,0.5), nrow=1))
             )

cm <- causl_model(forms, fam, pars=pars)

test_that("causl_model() behaves as expected", {
  expect_s3_class(cm, "causl_model")
  expect_error(causl_model(forms, fam, pars[-1]), regexp = "Variable Z1")
  expect_error(causl_model(forms, fam[-4], pars), regexp = "Should be a family entry for each set of formulas")
})


