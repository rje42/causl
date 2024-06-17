## check that glm_dens gives the same answers with a causl_fam as a number
set.seed(123)

x <- glm_sim(1, eta = rep(0,100), phi = 1, link = "identity")
q0 <- glm_dens(x=x, eta = rep(0,100), phi = 1)
q1 <- glm_dens(x=x, eta = rep(0,100), family= gaussian_causl_fam(), phi = 1)

test_that("Gaussian distribution works", {
  expect_equal(q0, q1)
})


x <- factor(1:3)
q1 <- glm_dens(x=x, eta = cbind(rep(0,3), rep(0.25,3)), family=categorical_causl_fam(),
               other_pars = list(nlevel=3))
p <- c(theta_to_p_cat(c(0,0.25)))
cp <- cumsum(p)
q0 <- list(u=cp, ld=log(p))

test_that("Categorical distribution works", {
  expect_equal(q0, q1)
})
