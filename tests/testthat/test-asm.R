test_that("test empty model", {


})


test_that("gaussian covariate, gaussian error", {
  n = 100
  p = 2
  X = matrix(rnorm(n*p), n, p)
  beta0 = rep(0, p)
  noise = rnorm(n)
  y = X %*% beta0 + noise
  mydf = data.frame(y, X)
  colnames(mydf) = c("y", paste0("X", 1:p))
  asm_res = asm(y ~ ., data = mydf)

  betahat = coef(asm_res)
  expect_equal(length(betahat), p+1)

  ## expect that the true beta is close to the estimated beta
  expect_true(all(abs(betahat[2:(p+1)] - beta0) < 1))
})



