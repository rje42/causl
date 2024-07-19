fams <- list(c(1,3), list(binomial_causl_fam()), c(1,5))
vars <- list(c("Z1","Z2"), "X1", c("Y1", "Y2"))

link1 <- list(Z1="log", Y2="probit")

out1 <- list(c(Z1="identity", Z2="log"), c(X1="logit"), c(Y1="identity", Y2="logit"))
out2 <- list(c(Z1="log", Z2="log"), c(X1="logit"), c(Y1="identity", Y2="probit"))

test_that("link_setup() works for different inputs", {
  expect_equal(link_setup(family=fams, vars=vars), out1)
  expect_equal(link_setup(link=link1, family=fams, vars=vars), out2)
})

link2 <- list(X1="probit")

test_that("link_setup() returns right warnings", {
  expect_warning(link_setup(link=link2, family=fams, vars=vars), "Provided link does not match that in family. Using family link")
})

link3 <- list(X2="probit")

test_that("link_setup() returns right errors", {
  expect_error(link_setup(link=link3, family=fams, vars=vars), "Link named for variable not in the model")
})
