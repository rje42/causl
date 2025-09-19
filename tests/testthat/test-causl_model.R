cm <- causl_model(formulas = list(Z ~ 1, X~ Z, Y ~ 1, ~ 1),
                  family = list(3, 1, 1, 1),
                  pars = list(Z = list(beta = 1, phi = 1),
                              X = list(beta = c(0.5, 1), phi = 1),
                              Y = list(beta = 0, phi = 1),
                              cop = list(beta = 1)))

test_that("causl_model works", {
  expect_equal(cm$pars$X$beta, c(0.5, 1))
})


cm2 <- modify(cm, pars = list(cop = list(beta = 1.5)))

test_that("modify.causl_model works", {
  expect_equal(cm2$pars$cop$Y$Z$beta, 1.5)
})
