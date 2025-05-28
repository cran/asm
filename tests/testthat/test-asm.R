
test_that("test boundary cases", {

  # test no covariate, only intercept
  n = 10
  Y = rnorm(n)
  asm_fit = asm(Y ~ 1)
  expect_true(length(asm_fit$betahat) == 1)
  expect_true(abs(asm_fit$betahat) < 1)

  asm_fit = asm(Y ~ 1, symmetric=FALSE)
  expect_true(abs(asm_fit$betahat - mean(Y)) < 1e-5)

  # test null model
  asm_fit = asm(Y ~ - 1)
  expect_true(is.null(asm_fit$betahat))

})



test_that("gaussian covariate, gaussian error: basic functionality", {
  n = 200
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
  expect_true(abs(betahat[1]) < 1)

  asm_res2 = asm.fit(cbind(rep(1, n), X), y)
  betahat2 = asm_res2$betahat
  expect_true(all(abs(betahat - betahat2) < 1e-6))
})




test_that("cauchy error: prediction interval", {
  n = 600
  p = 6

  X = matrix(rnorm(p*n, -5, 2), n, p)
  beta0 = rnorm(p)

  Y = X %*% beta0 + rcauchy(n)
  train_idx = sample(1:n, size = round(0.6*n))
  X_train = X[train_idx, ]; Y_train = Y[train_idx]
  X_test = X[-train_idx, ]; Y_test = Y[-train_idx]

  df = data.frame(Y = Y_train, X = X_train)
  df_test = data.frame(X = X_test)

  asm_fit = asm(Y ~ ., data=df)

  preds = predict(asm_fit, newdata=df_test, interval="prediction", level=0.9)

  expect_true(nrow(preds) == length(Y_test) && ncol(preds) == 3)

  coverage = mean(Y_test < preds[, 3] & Y_test > preds[, 2])

  expect_true(mean(coverage) > 0.6)
})



test_that("Cauchy error: asymmetric loss", {
  n = 200
  p = 5
  X = matrix(rnorm(p*n, 1, 1), n, p)
  beta0 = rep(1, p)
  Y = X %*% beta0 + rcauchy(n)

  df = data.frame(X=X, Y=Y)
  asm_fit = asm(Y ~ X, data=df, symmetric=FALSE)

  expect_true(length(asm_fit$betahat) == p+1)
  expect_true(all(abs(asm_fit$betahat[-1] - 1) < 1))

  asm_fit2 = asm(Y ~ X, data=df, symmetric=FALSE, intercept.selection="median")

  expect_true(length(asm_fit2$betahat) == p+1)
  expect_true(all(abs(asm_fit2$betahat[-1] - asm_fit$betahat[-1]) < 1e-5))

  expect_true(abs(asm_fit2$betahat[1]) < 1)
})


test_that("large p", {

  n = 50
  p = 40

  X = matrix(rnorm(p*n, 0, 1), n, p)
  beta0 = rep(1/p, p)
  Y = X %*% beta0 + rcauchy(n)

  df = data.frame(X=X, Y=Y)
  asm_fit = asm(Y ~ ., data=df)

})


