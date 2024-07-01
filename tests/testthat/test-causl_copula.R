set.seed(123)

n <- 10
rho <- matrix(c(1,0.2,0.1,
                0.2,1,0.05,
                0.1,0.05,1), 3, 3)
samp <- rGaussCop(n=n, Sigma=rho)

dflt <- dGaussCop(samp, Sigma=rho)
gcc <- gaussian_causl_cop()

test_that("causl_copula functions work", {
  expect_equal(gcc$ddist(samp, Sigma=rho), dflt)
})
