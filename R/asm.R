#' @import fdrtool
#' @import pracma
#' @import MASS
#' @import quantreg
#' @import Iso
#' @import stats
#' @import grDevices
#' @import graphics

library(fdrtool)
library(pracma)
library(Iso)
library(MASS)
library(quantreg)

## IN: psi
## OUT: negative log-density (normalized) corresponding to psi
##
getPhiPDF <- function(pts_evals, psi_evals) {
  K <- length(pts_evals)
  phi <- c(0, cumsum((psi_evals[-1] + psi_evals[-K]) / 2 * diff(pts_evals)))
  phi = phi - max(phi)
  pdf <- exp(phi)
  pdf <- pdf / sum((pdf[-1] + pdf[-K]) * diff(pts_evals) / 2)
  phi <- log(pdf)
  return(phi)
}

## LCM(J) = least concave majorant of J
##
## Computes
## psi* = LCM(densities circ Q)^{(R)} circ F
##
## where Q, F are the quantile and distribution functions respectively corresponding to
## "densities", represented as evaluations on "grid"
##
## INPUT: densities -- vector of density values
##        grid -- points at which the density is evaluated
##        k -- number of discrete points between [0,1] at which
##             to evaluate Jhat and (psi^*)
##
## OUTPUT: psi -- function  R -> R
##         psi_deriv -- function  R -> R
##
isotonize_score_given_density <- function(grid, densities, J_eval = NULL, k = 3000,
                                          return_deriv = TRUE,
                                          debug = FALSE,
                                          lcm = FALSE, k1=1, k2=k, ptm) {
  if (debug) browser()
  K <- length(grid)

  # Interpolate between density values on the grid
  f <- approxfun(grid, densities, method = "linear", rule = 2)
  masses <- (densities[-1] + densities[-K])/2 * diff(grid)
  Fx <- c(0, cumsum(masses))
  masses <- masses / Fx[K]
  Fx <- Fx / Fx[K]

  Q <- approxfun(Fx, grid, method = "linear", rule = 2, ties = max)

  u <- 0:k / k # customise this for light-tailed densities for plotting purposes?
  J <- f(Q(u))

  # Evaluate J at a given point J_eval (used to estimate the intercept standard
  # error when intercept.selection = "quantile" in asm.fit)
  if (!is.null(J_eval)) J_eval <- f(Q(J_eval))

  ## plotting LCM of the density quantile function
  if (lcm) {
    lcm <- gcmlcm(u, J, type = "lcm")
    if (debug) {
      plot(u, J, type = 'l', xlab = 'u', ylab = 'J', xlim = c(0, 1))
      lines(lcm$x.knots, lcm$y.knots, col = 2)
    }
  }

  init_scores <- diff(J) * k # psihat circ Qhat(0, 1/k, 2/k, ..., 1)

  ## init_scores_sub = init_scores[k1:k2]
  domain = (Q(u[-1]) + Q(u[-(k + 1)])) / 2

  fn_vals = pava(init_scores, decreasing = TRUE,
                 long.out = FALSE, stepfun = FALSE)

  fn_vals = fn_vals[k1:k2]
  domain = domain[k1:k2]

  jumps = diff(fn_vals) != 0
  jumps_lshift = c(jumps, TRUE)
  jumps_rshift = c(TRUE, jumps)
  non_redundant_index <- jumps_lshift | jumps_rshift

  domain <- domain[non_redundant_index]
  fn_vals <- fn_vals[non_redundant_index]

  yleft <- max(fn_vals[1], 1e-7)
  yright <- min(fn_vals[length(fn_vals)], -1e-7)
  psi <- approxfun(domain, fn_vals, method = "linear",
                   yleft = yleft, yright = yright)

  if (return_deriv) {
    deriv_vals <- diff(fn_vals) / diff(domain)
    ## Constant extension (instead of 0) outside the support of grid
    ## for numerical stability in Newton's method later
    # deriv_vals <- c(deriv_vals, 0)
    deriv_vals <- c(deriv_vals, deriv_vals[length(deriv_vals)])
    psi_deriv <- approxfun(domain, deriv_vals,
                           method = "constant", f = 0,
                           # yleft = 0, yright = 0,
                           rule = 2)
    if (debug) {
      plot(domain, fn_vals, type = 'l')
      plot(domain, deriv_vals, type = 'l')
    }
    if (lcm) {
      return(list(psi = psi, psi_deriv = psi_deriv, J = J, J_eval = J_eval,
                  Jhat = lcm))
    }
    return(list(psi = psi, psi_deriv = psi_deriv, J_eval = J_eval))
  } else {
    return(psi)
  }
}


## Kernel-based antitonic projected score estimates
## truncation not implemented
## ... are passed to density();
## see density() for more details (kernel, bw, etc.)
kde_decr_score_est <- function(residuals, quantile = NULL, k = 3000, kernel_pts = 2^15, ...) {

  residuals <- sort(residuals)

  supp1 = range(residuals) + 3*bw.nrd0(residuals) * c(-1, 1)
  supp2 = median(residuals) + IQR(residuals) * c(-20, 20)
  supp = c(max(supp1[1], supp2[1]), min(supp1[2], supp2[2]))

  kde <- density(residuals, n = kernel_pts,
                from=supp[1], to=supp[2], ...)

  n = length(residuals)

  k = max(k, 2 * n)
  if (n > 2) {
    skip = 2 * floor(k / n)
  } else {
    skip = 0
  }

  iso_res = isotonize_score_given_density(kde$x, kde$y, J_eval = quantile,
            k=k, k1=1+skip, k2=k-skip)

  return(iso_res)
}


### General root finding ####
# Gives the betahat that satisfies sum_i Xi psi(Yi-<Xi,betahat>) = 0
linear_regression_newton_fn <- function(betainit, X, Y, fn_psi, fn_psi_deriv,
                                        max_iter = 55, verbose=FALSE) {
  betahat <- betainit
  d <- ncol(X)
  N <- nrow(X)
  CONV_THRESH <- (d * N)^0.5 * 1e-14

  hess_cond_tol = 1e6

  if (d == 0){
    return(list(betahat = matrix(0, 0, 1), conv=1))
  }

  residuals <- Y - X %*% betahat

  psi_vec <- fn_psi(residuals)
  psi_deriv_vec <- fn_psi_deriv(residuals)
  grad <- - t(X) %*% psi_vec

  Xpsi <- X * as.vector(psi_deriv_vec)
  H <- t(Xpsi) %*% X

  for (l in 1:max_iter) {

    alpha <- 1

    cond_num = kappa(H)

    if (cond_num < hess_cond_tol)
      update <- solve(H - diag(nrow(H)), grad)
    else
      update <- grad

    while (TRUE) {
      betatemp <- betahat - alpha * update

      resid_temp <- Y - X %*% betatemp
      psi_vec_temp <- fn_psi(resid_temp)
      grad_temp <- - t(X) %*% psi_vec_temp

      if (sum(grad_temp^2) <= sum(grad^2) || alpha < 1e-4) {
        betahat <- betatemp
        grad <- grad_temp
        break
      } else
        alpha <- 0.8 * alpha
    }

    if (verbose){
      print(sqrt(sum(grad^2)))
    }

    if (sqrt(sum(grad^2)) < CONV_THRESH) {
      return(list(betahat = betahat, conv = 1))
    }

    psi_deriv_vec = fn_psi_deriv(resid_temp)

    Xpsi <- X * as.vector(psi_deriv_vec)
    H <- t(Xpsi) %*% X

  }
  return(list(betahat = betahat, conv = 0))
}


## Jointly computes theta_hat and mu_hat
## NOT used currently
asm_regression_test <- function(betapilot, residuals, X, Y, quantile = NULL,
                                k = 3000, kernel_pts=2^15,
                                est_info = FALSE, est_score_obj = FALSE,
                                max_iter = 65, return_fn = FALSE,
                                bw = "nrd0", kernel = "gaussian") {

  res = kde_decr_score_est(residuals, quantile = quantile, k = k, kernel_pts = kernel_pts,
                           bw = bw, kernel = kernel)

  fn_psi_deriv = res$psi_deriv
  fn_psi = res$psi

  n = length(Y)
  Xint = cbind(rep(1, n), X)

  root = linear_regression_newton_fn(betapilot, Xint,
                                     Y,
                                     fn_psi, fn_psi_deriv,
                                     max_iter = max_iter)

  return_list = list(thetahat = root$betahat[-1],
                     muhat = root$betahat[1],
                     conv = root$conv)

  if (return_fn) {
    return_list$psi = fn_psi
    return_list$psi_deriv = fn_psi_deriv
  }
  if (est_info || est_score_obj) {
    resids = Y - X %*% root$betahat[-1]
    temp1 = mean(fn_psi(resids)^2)
    if (!is.null(quantile)) return_list$quantile <- res$J_eval
    return_list$info_asm = temp1
    if (est_score_obj) {
      temp2 = 2 * fn_psi_deriv(resids)
      return_list$score_obj = mean(temp2) + temp1
    }
  }

  return(return_list)
}


## Estimate antitonic score psi from residuals, then
## compute
##
## min_theta sum_i ell( Yi - offset - (Xi - Xbar)' theta )
##
## where offset = betapilot[1] + Xbar' betapilot[-1]
##
## REQUIRE: betapilot[1] is the intercept;
##          betapilot has dimension ncol(X) + 1
##
asm_regression <- function(betapilot, residuals, X, Y, quantile = NULL,
                           k = 3000, kernel_pts=2^15,
                           est_info = FALSE, est_score_obj = FALSE,
                           max_iter = 65, return_fn = FALSE,
                           bw = "nrd0", kernel = "gaussian") {

  res = kde_decr_score_est(residuals, quantile = quantile, k = k, kernel_pts = kernel_pts,
                          bw = bw, kernel = kernel)

  fn_psi_deriv = res$psi_deriv
  fn_psi = res$psi

  Xbar <- colMeans(X)
  Xcentered <- sweep(X, 2, Xbar, "-")
  offset = sum(Xbar * betapilot[-1]) + betapilot[1]
  Y_offset = Y - offset

  root = linear_regression_newton_fn(betapilot[-1], Xcentered,
                                     Y_offset,
                                     fn_psi, fn_psi_deriv,
                                     max_iter = max_iter)

  return_list = list(thetahat = root$betahat, conv = root$conv)

  if (return_fn) {
    return_list$psi = fn_psi
    return_list$psi_deriv = fn_psi_deriv
  }
  if (est_info || est_score_obj) {
    resids = Y_offset - Xcentered %*% root$betahat
    temp1 = mean(fn_psi(resids)^2)
    if (!is.null(quantile)) return_list$quantile <- res$J_eval
    return_list$info_asm = temp1
    if (est_score_obj) {
      temp2 = 2 * fn_psi_deriv(resids)
      return_list$score_obj = mean(temp2) + temp1
    }
  }

  return(return_list)
}


## one iteration of asm regression with symmetric loss
##
## min_beta sum_i ell( Yi  - Xi' beta )
##
## ASSUME: betapilot has dimension ncol(X)
##
asm_regression_sym <- function(betapilot, residuals, X, Y,
                            k = 3000, kernel_pts=2^15,
                            est_info = FALSE, est_score_obj = FALSE,
                            max_iter = 100, return_fn = FALSE,
                            bw = "nrd0", kernel="gaussian") {

  res = kde_decr_score_est(residuals, k = k, kernel_pts = kernel_pts,
                           bw = bw, kernel = kernel)

  fn_psi = function(x) {
    return((res$psi(x) - res$psi(-x)) / 2)
  }
  fn_psi_deriv = function(x) {
    return((res$psi_deriv(x) + res$psi_deriv(-x)) / 2)
  }

  root = linear_regression_newton_fn(betapilot, X, Y,
                                     fn_psi, fn_psi_deriv,
                                     max_iter = max_iter)

  return_list = list(betahat = root$betahat, conv = root$conv)

  if (return_fn) {
    return_list$psi = fn_psi
    return_list$psi_deriv = fn_psi_deriv
  }
  if (est_info || est_score_obj) {
    resids = Y - X %*% root$betahat
    temp1 = mean(fn_psi(resids)^2)
    return_list$info_asm = temp1
    if (est_score_obj) {
      temp2 = 2 * fn_psi_deriv(resids)
      return_list$score_obj = mean(temp2) + temp1
    }
  }

  return(return_list)
}








#' Fit a linear regression model via antitonic score matching
#' @description Performs linear regression via M-estimation with respect to a data-driven convex loss function
#'
#' @param X design matrix
#' @param Y response vector
#' @param betapilot initial estimate of the regression coefficients:
#' can be "LAD", "OLS" or a vector of coefficients
#' @param symmetric logical; if TRUE, estimate a symmetric loss function
#' @param alt_iter number of iterations of the alternating procedure:
#' when alt_iter == 1, this function is equivalent to asm_regression
#' @param error_quantile quantile of the residuals to be returned as intercept
#' Used only if intercept.selection = "quantile"
#' If error_quantile = 0.5, then the intercept is the median of the residual
#' Ignored if symmetric = TRUE
#' @param intercept.selection mean, median, or quantile of the residuals
#' If intercept.selection = "quantile", then error_quantile specifies the quantile value
#' Ignored if symmetric == TRUE
#' @param k the density quantile function is evaluated at (0, 1/\code{k}, 2/\code{k}, ..., 1)
#' @param max_iter maximum number of iterations for the damped Newton–Raphson algorithm
#' when minimizing the convex loss function
#' @param kernel_pts number of points at which the kernel density estimate is evaluated,
#' i.e. the parameter "n" in density()
#' @param bw bandwidth for kernel density estimation
#' i.e. the parameter "bw" in density()
#' @param kernel kernel for kernel density estimation
#' i.e. the parameter "kernel" in density()
#' @param verbose logical; if TRUE, print optimization progress
#' @param ... additional arguments to ensure compatibility with generic functions
#'
#' @return \code{asm} class object containing the following components:
#' \describe{
#'  \item{\code{betahat}:}{vector of estimated coefficients}
#'  \item{\code{std_errs}:}{vector of standard errors of the estimated coefficients}
#'  \item{\code{fitted.values}:}{fitted values}
#'  \item{\code{residuals}:}{residuals}
#'  \item{\code{zvals}:}{z-values}
#'  \item{\code{sig_vals}:}{p-values}
#'  \item{\code{info_asm}:}{estimated antitonic information}
#'  \item{\code{I_mat}:}{estimated antitonic information matrix}
#'  \item{\code{Cov_mat}:}{asymptotic covariance matrix of the estimated coefficients}
#'  \item{\code{psi}:}{estimated antitonic score function}
#'  \item{\code{symmetric}:}{logical; indicating whether the loss is constrained to be symmetric}
#' }
#'
#' @export
#'
#' @examples
#' n <- 1000 ; d <- 2
#' X <- matrix(rnorm(n * d), n, d)
#' Y <- X %*% c(2, 3) + 1 + rnorm(n)
#' asm.fit(X, Y)
#'
#' Y <- X %*% c(2, 3) + rexp(n)
#' asm.fit(X, Y, symmetric=FALSE)
asm.fit <- function(X, Y, betapilot = "LAD", error_quantile = NULL,
                    symmetric = TRUE,
                    alt_iter = 2, intercept.selection = "mean",
                    k = 3000, max_iter = 200, kernel_pts = 2^15,
                    bw = "nrd0", kernel = "gaussian", verbose=FALSE, ...) {

  THRESH = 1e-6

  if (!symmetric){
    # Allow for general quantiles; q = 0.5 is the default
    q <- error_quantile
    if (is.null(q)) q <- 0.5
    if (verbose) cat("\nNote: Using asymmetric loss - estimating the intercept by the",
                     error_quantile, intercept.selection, "of the residuals\n")
  }

  X = as.matrix(X)

  if (ncol(X) == 0){
    model_res = list(null = TRUE)
    class(model_res) = "asm"
    return(model_res)
  }
  if (nrow(X) < 3){
    stop("The input must contain at least 3 rows.")
  }
  if (nrow(X) )


  if (symmetric)
    intercept.selection = "none"
  if (intercept.selection == "median"){
    q = 0.5
    intercept.selection = "quantile"
  }

  ## check if X already contains intercept as column of all 1s
  ## Xtilde used for asymmetric loss: equals X if no intercept, otherwise it is
  ##    X without the column of all 1s
  if (any(X[, 1] != 1)) {
    has_intercept = FALSE
    Xtilde = X

    if (is.null(colnames(X))) {
      colnames(X) = paste0("X", 1:ncol(X))
    }
  } else {
    has_intercept = TRUE

    Xtilde = as.matrix(X[, -1])

    if (is.null(colnames(X))) {
      if (ncol(X) == 1) {
        colnames(X) = c("(Intercept)")
      } else {
      colnames(X) = c("(Intercept)", paste0("X", 1:(ncol(X)-1)) )
      }
    }
  }
  n = nrow(X)

  if (betapilot[1] == "OLS") {
    betapilot = solve(t(X) %*% X, t(X) %*% Y)
  } else if (betapilot[1] == "LAD") {
    # Suppress 'solution may not be unique' warnings or make them more informative
    LAD_fit <- suppressWarnings(quantreg::rq.fit(x=X, y=Y, tau=0.5))
    betapilot = LAD_fit$coefficients
  }
  beta_init = betapilot

  residuals = Y - X %*% betapilot

  scoreobj = Inf
  res_prev = NULL

  ## if asymmetric and no intercept, then
  ## betapilot input to asm_regression should always
  ## be length d+1
  ##
  ## betapilot input to asm_regression_sym should always
  ## be length d
  if (!has_intercept && !symmetric) {
    betapilot = c(0, betapilot)
  }


  for (it in 1:alt_iter) {

    if (!symmetric){


      res_asm = asm_regression(betapilot, residuals, quantile = q,
                               Xtilde, Y,
                            k = k, kernel_pts = kernel_pts,
                            max_iter = max_iter,
                            bw = bw, kernel = kernel,
                            est_score_obj = TRUE, return_fn = TRUE)




      thetahat_asm = res_asm$thetahat
      resid_uncentered = Y - Xtilde %*% thetahat_asm

      if (intercept.selection == "mean" && has_intercept) {
        muhat = mean(resid_uncentered)
      } else if (intercept.selection == "quantile" && has_intercept){
        muhat = quantile(resid_uncentered, q)
      } else {
        muhat = 0
      }
      betapilot = c(muhat, thetahat_asm)

    } else {

      res_asm = asm_regression_sym(betapilot, residuals, X, Y,
                            k = k, kernel_pts = kernel_pts,
                            max_iter = max_iter,
                            bw = bw, kernel = kernel,
                            est_score_obj = TRUE, return_fn = TRUE)

      betahat_asm = res_asm$betahat
      betapilot = betahat_asm

    }

    info_asm = res_asm$info_asm

    if (!symmetric){
      if (has_intercept) betahat = c(muhat, thetahat_asm)
      if (!has_intercept) betahat = thetahat_asm
    } else {
      betahat = betahat_asm
    }

    betahat = as.vector(betahat)
    residuals = Y - X %*% betahat

    ## compare score matching objective
    scoreobj_new = res_asm$score_obj

    if (scoreobj - scoreobj_new < THRESH && it > 1) {
      if (verbose)
       print(paste("Alternation finished at iteration", it))
      break
    } else {
      ## continue to next iteration
      scoreobj = scoreobj_new
      res_prev = res_asm
    }

  }

  if (!symmetric){
    I_mat = info_asm * var(Xtilde)
  } else {
    I_mat = info_asm * (1/n) * t(X) %*% X
  }

  if (ncol(Xtilde) > 0 || symmetric){
    v22 = solve(I_mat)
  }

  # Use general formula to calculate the full asymptotic covariance matrix,
  # including the intercept
  # tau is an estimate of (\int \zeta^2 p_0)/(\int p_0 d\zeta)^2,
  # where zeta(z) = z when intercept.selection = "mean" and
  # zeta(z) = q - 1_{z < 0} when intercept.selection = "quantile" and
  # error_quantile = q
  if (has_intercept && !symmetric) {
    Xbar <- colMeans(Xtilde)
    if (intercept.selection == "mean") {
      tau <- var(Y - Xtilde %*% thetahat_asm)
    } else if (intercept.selection == "quantile") {
      tau <- q*(1 - q) / res_asm$quantile^2
    }
    tau <- as.numeric(tau)

    if (ncol(Xtilde) == 0){
      Cov_mat = matrix(tau, 1, 1)
      std_errs = rep(tau, 1)
      I_full = matrix(1/tau, 1, 1)
    } else {
      I_full <- info_asm * var(X) + outer(colMeans(X), colMeans(X)) / tau

      v11 = tau + Xbar %*% v22 %*% Xbar
      v12 = Xbar %*% v22
      Cov_mat <- I_full
      Cov_mat[1, 1] = v11
      Cov_mat[-1, -1] = v22
      Cov_mat[1, -1] = -v12
      Cov_mat[-1, 1] = -t(v12)

      stopifnot( max(abs(solve(I_full) - Cov_mat)) < 1e-5 )

      std_errs = sqrt(diag(Cov_mat) / n)
    }
  } else {
    Cov_mat = v22
    std_errs = sqrt(diag(Cov_mat) / n)
    I_full <- I_mat
  }


  names(std_errs) = colnames(X)
  names(betahat) = colnames(X)

  zvals = abs(betahat) / std_errs
  sig_vals = pnorm(zvals, lower.tail = FALSE) * 2

  fitted.values = X %*% betahat
  residuals = Y - fitted.values

  model_res <- list(
    betahat = betahat,
    beta_init = beta_init,
    std_errs = std_errs,
    fitted.values = fitted.values,
    residuals = residuals,
    zvals = zvals,
    sig_vals = sig_vals,
    info_asm = info_asm,
    I_mat = I_mat,
    I_full = I_full,
    Cov_mat = Cov_mat,
    psi = res_asm$psi,
    conv = res_asm$conv,
    symmetric = symmetric,
    xbar = apply(Xtilde, 2, mean),
    intercept.selection = intercept.selection
  )

  class(model_res) = "asm"

  return(model_res)

}


#' Linear regression via antitonic score matching
#'
#' @description Performs linear regression with a data-driven convex loss function
#'
#' @param formula regression formula
#' @param data input data frame
#' @param ... additional arguments for asm.fit
#'
#' @return \code{asm} class object containing the following components:
#' \describe{
#'  \item{\code{betahat}:}{vector of estimated coefficients}
#'  \item{\code{std_errs}:}{vector of standard errors of the estimated coefficients}
#'  \item{\code{fitted.values}:}{fitted values}
#'  \item{\code{residuals}:}{residuals}
#'  \item{\code{zvals}:}{z-values}
#'  \item{\code{sig_vals}:}{p-values}
#'  \item{\code{info_asm}:}{antitonic information}
#'  \item{\code{I_mat}:}{estimated antitonic information matrix}
#'  \item{\code{Cov_mat}:}{covariance matrix of the estimated coefficients}
#'  \item{\code{psi}:}{estimated antitonic score function}
#' }
#'
#'
#' @export
#'
#' @examples
#' asm(mpg ~ cyl + hp + disp, data=mtcars)
#'
#' asm(mpg ~ cyl + hp + disp, data=mtcars, symmetric=FALSE)
#'
#' n <- 1000 ; d <- 2
#' X <- matrix(rnorm(n * d), n, d)
#' Y <- X %*% c(2, 3) + 1 + rnorm(n)
#' asm(Y ~ X - 1)
#'
#' Y <- X %*% c(2, 3) + rchisq(n, 6) - qchisq(0.4, 6)
#' asm(Y ~ X, symmetric=FALSE, intercept.selection="quantile", error_quantile=0.4)
#'
#' Y <- X %*% c(2, 3) + rcauchy(n)
#' asm(Y ~ X, symmetric=FALSE, intercept.selection="median")
#'
asm <- function(formula, data = NULL, ...) {
  Y = model.response(model.frame(formula, data))
  Xint = model.matrix(formula, data)

  cl = match.call()

  model = asm.fit(Xint, Y, ...)
  model$formula = formula
  model$call = cl
  return(model)
}


#' Short description of a fitted \code{asm} regression model
#'
#' @description Outputs estimated coefficients and standard errors
#'
#' @param x asm object
#' @param ... additional arguments to ensure compatibility with the generic function print()
#' @return No return value, called for its side effect
#'
#' @export
#'
#' @examples
#' model = asm(mpg ~ cyl + hp + disp, data=mtcars)
#' print(model)
print.asm <- function(x, ...) {

  cat("\nCall:", paste(deparse(x$call), sep = "\n", collapse = "\n"), sep = "")

  if ("null" %in% names(x)){
    cat("\n\nNo coefficients to fit\n")
    return(invisible(x))
  }

  cat("\nEstimates:\n")
  print(x$betahat)
  cat("\nStandard errors:\n")
  print(x$std_errs)
}


#' Summary of an \code{asm} regression model
#'
#' @description Outputs estimated coefficients, standard errors and p-values based on a fitted \code{asm} regression model
#'
#' @param object asm object
#' @param ... additional arguments to ensure compatibility with the generic function summary()
#'
#' @return \code{summary.asm} class object containing the following components:
#' \describe{
#'  \item{\code{coefficients}:}{estimated coefficients, standard errors, z-values and p-values}
#'  \item{\code{residuals}:}{residuals of the fitted model}
#'  \item{\code{call}:}{call to the \code{asm} function}
#' }
#'
#' @export
#'
#' @examples
#' model = asm(mpg ~ cyl + hp + disp, data=mtcars)
#' summary(model)
summary.asm <- function(object, ...) {

  if ("null" %in% names(object)){
    ans = list()
    ans$call = object$call
    ans$null = TRUE
    class(ans) = "summary.asm"
    return(ans)
  }

  # Extract coefficients, standard errors, z-values and p-values
  est <- object$betahat
  stderr <- object$std_errs
  zvalue <- object$zvals
  pvalue <- object$sig_vals

  ans = list()
  ans$call = object$call
  ans$residuals = object$residuals

  ans$coefficients <- cbind(Estimate = est, `Std. Error` = stderr,
        `z value` = zvalue, `Pr(>|z|)` = pvalue)

  class(ans) <- "summary.asm"
  return(ans)
}


#' Print summary of the \code{asm} regression model
#' @description Prints the summary of a fitted \code{asm} regression model
#' @param x summary.asm object
#' @param digits number of digits to print
#' @param signif.stars logical; if TRUE, 'significance stars' are printed
#' @param concise logical; if TRUE, the output is concise
#' @param ... additional arguments to ensure compatibility with the generic function print()
#' @return No return value
#'
#' @export
#' @examples
#' model = asm(mpg ~ cyl + hp + disp, data=mtcars)
#' print(summary(model))
#'
print.summary.asm <- function (x, digits = max(3L, getOption("digits") - 3L),
          signif.stars = getOption("show.signif.stars"), concise = FALSE, ...) {
  cat("\nCall:", if (!concise) "\n" else " ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
            if (!concise) "\n\n", sep = "")

  if ("null" %in% names(x)){
    cat("No coefficients to fit\n")
    return(invisible(x))
  }

  resid <- x$residuals

  if (!concise) cat("Residuals:\n")
  zz = zapsmall(quantile(resid), digits + 1L)
  nam = c("Min", "1Q", "Median", "3Q", "Max")
  rq = structure(zz, names = nam)

  if (!concise) print(rq, digits = digits, ...)

  if (!concise) cat("\nCoefficients:\n")
  coefs <- x$coefficients

  printCoefmat(coefs, digits = digits, signif.stars = signif.stars, signif.legend = (!concise),
                         na.print = "NA", eps.Pvalue = if (!concise) .Machine$double.eps else 1e-4, ...)

  cat("\n")
  invisible(x)
}


#' Coefficients of an \code{asm} regression model
#' @description Outputs the coefficients of a fitted \code{asm} regression model
#' @param object asm object
#' @param ... additional arguments to ensure compatibility with the generic function coef()
#' @return vector of coefficients of the \code{asm} regression model
#'
#' @export
#' @examples
#' model = asm(mpg ~ cyl + hp + disp, data=mtcars)
#' coef(model)
#'
coef.asm <- function(object, ...) {
  return(object$betahat)
}


#' Residuals from an \code{asm} regression model
#' @description Outputs the residuals (on the training data) from a fitted \code{asm} regression
#'              model
#' @param object asm object
#' @param ... additional arguments to ensure compatibility with the generic function residuals()
#' @return vector of residuals from the \code{asm} regression model
#' @export
#' @examples
#' model = asm(mpg ~ cyl + hp + disp, data=mtcars)
#' residuals(model)
#'
residuals.asm <- function(object, ...) {
  return(object$residuals)
}

#' Confidence intervals for coefficients
#' in an \code{asm} regression model
#' @description Computes confidence intervals for individual regression coefficients based on a fitted \code{asm} regression model
#' @param object asm object
#' @param parm parameters to calculate confidence intervals
#' @param level confidence level
#' @param ... additional arguments to ensure compatibility with the generic function confint()
#' @return matrix of confidence intervals for the regression coefficients
#'
#' @export
#' @examples
#' model = asm(mpg ~ cyl + hp + disp, data=mtcars)
#' confint(model)
#'
confint.asm <- function(object, parm, level = 0.95, ...) {
  ses = object$std_errs
  cf = object$betahat
  pnames <- names(cf)

  if (missing(parm))
    parm <- pnames
  else if (is.numeric(parm))
    parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a)

  a_str = c(sprintf("%.1f%%", a[1]*100), sprintf("%.1f%%", a[2]*100))
  fac <- qnorm(a)
  # pct <- .format_perc(a, 3)
  ci <- array(NA_real_, dim = c(length(parm), 2L),
              dimnames = list(parm, a_str))
  ci[] <- cf[parm] + ses[parm] %o% fac
  ci
}

#' Predict new responses using an \code{asm} regression model.
#'
#' @description Outputs predictions on new test data based on a fitted \code{asm} regression model. Also returns a confidence interval around
#' the conditional mean (if interval = "confidence") or predicted response (if interval = "prediction").
#'
#' @param object asm object
#' @param newdata new data frame
#' @param interval type of interval calculation,
#' either "none", "confidence" or "prediction". Default is "none".
#' @param level confidence level
#' @param debug boolean; enables debug mode
#' @param ... additional arguments to ensure compatibility with the generic function predict()
#'
#' @return matrix of predicted values
#' * if interval = "none", the matrix has one column of predicted values
#' * if interval = "confidence" or "prediction", the matrix has three columns: predicted value, lower bound, and upper bound
#'
#' @export
#'
#' @examples
#' model = asm(mpg ~ cyl + hp + disp, data=mtcars)
#' predict(model, newdata = data.frame(cyl = 4, hp = 121, disp = 80), interval = "prediction")
#'
#' n <- 1000
#' X <- rnorm(n)
#' beta <- 2
#' Y <- beta*X + rt(n,df=3)
#' asm_model <- asm(Y ~ X)
#' predict(asm_model, newdata = data.frame(X = 1), interval = "prediction")
#'
predict.asm <- function(object, newdata = NULL, interval = "none",
                         level = 0.95, debug = FALSE, ...) {
  if (is.null(newdata))
    return(object$fitted.values)

  if (debug) browser()

  terms <- setdiff(names(object$betahat), "(Intercept)")

  if (is.data.frame(newdata)) {
    if (all(terms %in% names(newdata)))
      X <- newdata[, terms]
    else stop("One or more model variables missing from newdata")
    if (!all(names(newdata) %in% terms)) warning("Ignoring some variables in newdata")
    X <- as.matrix(X, ncol = length(terms))
  }
  else {
    if (length(newdata) %% length(terms)) {
      stop("newdata must be a data frame or matrix, with the number of columns
         equal to the number of non-intercept coefficients of model$betahat")
    }
    X <- matrix(newdata, ncol = length(terms))
  }

  if (length(object$betahat) > length(terms)) {
    Xint = cbind(1, X)
    pred = Xint %*% object$betahat
    beta_has_intercept = TRUE
    theta = object$betahat[-1]
    resids = object$residuals + object$betahat[1]
  } else {
    pred = X %*% object$betahat
    beta_has_intercept = FALSE
    theta = object$betahat
    resids = object$residuals
  }


  if (interval == "none") {
    return(pred)
  }

  if (interval == "prediction"){
    cov = object$Cov_mat
    xbar = object$xbar

    if (beta_has_intercept){
      cov = cov[-1, -1]
    }

    n = length(resids)

    delta = (1 - level)/10
    upp_quant = level + (1 - level)/2 + delta
    low_quant = (1 - level)/2 - delta

    pred_low = X %*% theta + quantile(resids, low_quant)
    pred_upp = X %*% theta + quantile(resids, upp_quant)

    suppressWarnings({
    supp1 = range(resids) + 3*bw.nrd0(resids) * c(-1, 1)
    supp2 = median(resids) + IQR(resids) * c(-20, 20)
    supp = c(max(supp1[1], supp2[1]), min(supp1[2], supp2[2]))
    resid_kde = density(resids, n = 2^14, from=supp[1], to=supp[2])
    f = approxfun(resid_kde$x, resid_kde$y, method="linear", rule=2)
    K = length(resid_kde$x)
    masses = (resid_kde$y[-1] + resid_kde$y[-K])/2 * diff(resid_kde$x)
    Fx = c(0, cumsum(masses))
    masses = masses / Fx[K]
    Fx = Fx/Fx[K]
    Q = approxfun(Fx, resid_kde$x, method="linear", rule=2)
    })

    tau_lo = low_quant * (1-low_quant) / (f(Q(low_quant))^2)
    tau_upp = upp_quant * (1-upp_quant) / (f(Q(upp_quant))^2)

    tmp = diag(X %*% cov %*% t(X)) - 2 * X %*% cov %*% xbar
    var_common = tmp + as.numeric(t(xbar) %*% cov %*% xbar)

    se_upp = sqrt(var_common + tau_upp)
    se_low = sqrt(var_common + tau_lo)

    adj = qnorm(1-delta/2)

    predictor = cbind(pred, pred_low - se_low * adj/sqrt(n),
                      pred_upp + se_upp*adj/sqrt(n))
    colnames(predictor) <- c("fit", "lwr", "upr")
    return(predictor)
  }
}

#' Generate diagnostic plots for an \code{asm} regression model
#'
#' @description Generates plots of residuals vs fitted values, and the estimated convex loss
#' and antitonic score functions based on a fitted \code{asm} regression model
#'
#' @param x asm object
#' @param which a subset of the plots to be displayed
#' @param caption a list of captions for the plots
#' @param extend.ylim.f factor to extend the y-axis limits
#' for the residuals vs fitted plot
#' @param id.n number of residuals to label in the residuals vs fitted plot
#' @param labels.id labels for the residuals in the residuals vs fitted plot
#' @param label.pos position of the labels in the residuals vs fitted plot
#' @param ext.xlim.f factor to extend the x-axis limits for the convex loss
#' and antitonic score function plots
#' @param grid.length.f the number of grid points for the convex loss plot
#' is defined as grid.length.f * length(x$residuals)
#' @param ask logical; if TRUE, the user is asked before each plot
#' @param ... additional arguments to ensure compatibility with the generic function plot()
#'
#' @return No return value
#'
#' @export
#'
#' @examples
#' model = asm(mpg ~ cyl + hp + disp, data=mtcars)
#' plot(model)
#'
plot.asm <- function(x, which = c(1, 2, 3),
                     caption = list("Residuals vs fitted",
                                    "Convex loss function",
                                    "Antitonic score function"),
                     extend.ylim.f = 0.08, id.n = 3,
                     labels.id = rownames(x$residuals),
                     label.pos = c(4, 2), ext.xlim.f = 0.08,
                     grid.length.f = 10,
                     ask = prod(par("mfcol")) < length(which) && dev.interactive(),
                     ...){
  show <- rep(FALSE, 6)
  show[which] <- TRUE
  yh = x$fitted.values
  r = x$residuals
  psi = x$psi
  n <- length(r)
  # plot the residuals vs fitted values
  if (id.n > 0L) {
    if (is.null(labels.id))
      labels.id <- paste(1L:n)
    iid <- 1L:id.n
    show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
    text.id <- function(x, y, ind, adj.x = TRUE, usr = par("usr")) {
      labpos <- if (adj.x)
        label.pos[(x > mean(usr[1:2])) + 1L]
      else 3
      text(x, y, labels.id[ind], cex = 0.75, xpd = TRUE,
           pos = labpos, offset = 0.25)
    }
  }
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  panel = function(x, y, ...) panel.smooth(x, y, iter = 3, ...)
  if (show[1L]) {
    ylim <- range(r, na.rm = TRUE)
    if (id.n > 0)
      ylim <- extendrange(r = ylim, f = extend.ylim.f)
    dev.hold()
    plot(yh, r, xlab = "Fitted values", ylab = "Residuals",
         main = "Residuals vs Fitted",
         ylim = ylim)
    panel(yh, r)
    if (id.n > 0) {
      y.id <- r[show.r]
      y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
      text.id(yh[show.r], y.id, show.r)
    }
    abline(h = 0, lty = 3, col = "gray")
    dev.flush()
  }

  if (show[2L]) {
    ## plot the convex loss function
    xlim = range(r)
    xlim <- extendrange(r = xlim, f = ext.xlim.f)
    grid = seq(xlim[1], xlim[2], length.out = grid.length.f * n)
    evals = psi(grid)
    phi <- getPhiPDF(grid, evals)
    plot(grid, -phi, type = "l", main = "Convex loss function",
         xlab = "x", ylab = "loss")
  }

  if (show[3L]) {
    ## plot the antitonic score function
    xlim = range(r)
    xlim <- extendrange(r = xlim, f = ext.xlim.f)
    dev.hold()
    plot(psi, main = "Antitonic score function",
         ylab = "score",
         xlim = xlim)
    points(r, rep(0, length(yh)), col = "red", pch = 4, cex = 0.5)
    dev.flush()
  }
  invisible()
}
