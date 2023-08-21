## check that dGDCop2_sig works in the way that is intended
Sigma <- matrix(c(1,0.3,0.3,1), 2, 2)
forms <- list(Z ~ 1, X ~ Z, Y ~ X, ~ X)
pars <- list(Z = list(beta = 0, phi = 1),
             X = list(beta = c(0,0.5)),
             Y = list(beta = c(-0.25,0.5)),
             cop = list(list(beta=1)))

rfr

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
