formulas <- list(list(C ~ 1,
                      Z ~ C),
                 X ~ Z + C,
                 Y ~ X + C,
                 ~1)
family <- list(c(5,1),5,1,1) # binomial, normal, binomial, normal, normal
link <- list(list("logit", "identity"), "logit", "identity")
pars <- list(C = list(beta=0),
             Z = list(beta = c(-1/2,0.25), phi=0.5),
             X = list(beta = c(0,1/2,1/10)),
             Y = list(beta = c(0.5, 1, 2), phi = 1),
             cop = list(Y = list(Z = list(beta=0.8472979), 
                                 C = list(beta = 0))))
cm <- causl_model(formulas = formulas,
                  family = family, link = link, pars = pars)
estimator <- list(outcome_model = 
                    list(model = survey::svyglm,
                    formula = Y~X,
                    family = "gaussian",
                    se = NULL),
              
                  ipw_model = 
                    list(model = glm,
                         family = "binomial",
                         formula = X ~ Z + C))
power_table1 <- estimate_power(cm, n = 50, estimator = estimator)
power_table2 <- estimate_power(cm, n = 100, estimator = estimator)
test_that("power function works and works as expected", {
  expect_lt(power_table1$Power, power_table2$Power) 
})
  

