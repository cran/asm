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
##        grid -- domain of "densities"
##        k -- number of discrete points between [0,1] at which
##             to evaluate Jhat and (psi^*)?
##
## OUTPUT: psi -- function  R -> R
##         psi_deriv -- function  R -> R
##
isotonize_score_given_density <- function(grid, densities, k = 3000,
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
      return(list(psi = psi, psi_deriv = psi_deriv, J = J, Jhat = lcm))
    }
    return(list(psi = psi, psi_deriv = psi_deriv))
  } else {
    return(psi)
  }
}


## Kernel-based antitonic projected score estimates
## truncation not implemented
## ... are passed to density();
## see density() for more details (kernel, bw, etc.)
kde_decr_score_est <- function(residuals, k = 3000, kernel_pts = 2^15,
                               #  truncation_lim = NULL, set_to_zero = FALSE
                               is_sorted = FALSE,  ...) {
  if (!is_sorted) {
    residuals <- sort(residuals)
  }

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

  iso_res = isotonize_score_given_density(kde$x, kde$y,
            k = k, k1 = 1 + skip, k2 = k - skip)

  return(iso_res)
}


### General root finding ####
# Gives the betahat that satisfies sum_i Xi psi(Yi-<Xi,betahat>) = 0
linear_regression_newton_fn <- function(betainit, X, Y, fn_psi, fn_psi_deriv,
                                        max_iter = 55) {
  betahat <- betainit
  d <- ncol(X)
  N <- nrow(X)
  CONV_THRESH <- (d * N)^0.5 * 1e-14

  residuals <- Y - X %*% betahat

  psi_vec <- fn_psi(residuals)
  psi_deriv_vec <- fn_psi_deriv(residuals)
  grad <- - t(X) %*% psi_vec

  Xpsi <- X * psi_deriv_vec
  H <- t(Xpsi) %*% X

  for (l in 1:max_iter) {

    alpha <- 1
    update <- solve(H + diag(nrow(H)), grad)
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

    if (sqrt(sum(grad^2)) < CONV_THRESH) {
      return(list(betahat = betahat, conv = 1))
    }

    psi_deriv_vec = fn_psi_deriv(resid_temp)

    Xpsi <- X * psi_deriv_vec
    H <- t(Xpsi) %*% X

  }
  return(list(betahat = betahat, conv = 0))
}

### without cross fitting, return theta_hat, i.e. no intercept term
asm_regression <- function(betapilot, residuals, X, Y,
                           k = 3000, kernel_pts=2^15,
                           est_info = FALSE, est_score_obj = FALSE,
                           max_iter = 65, return_fn = FALSE,
                           bw = "nrd0", kernel="gaussian") {

  res = kde_decr_score_est(residuals, k = k, kernel_pts = kernel_pts,
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
    return_list$info_asm = temp1
    if (est_score_obj) {
      temp2 = 2 * fn_psi_deriv(resids)
      return_list$score_obj = mean(temp2) + temp1
    }
  }
  # if (est_info) {
  #   return_list$info_asm = mean(fn_psi(Y - X %*% root$betahat-betapilot[1])^2)
  # }

  return(return_list)
}










#' Fit a linear regression model via antitonic score matching
#' @description Performs linear regression via M-estimation with respect to a data-driven convex loss function
#'
#' @param X design matrix
#' @param Y response vector
#' @param betapilot initial estimate of the regression coefficients:
#' can be "LAD", "OLS" or a vector of coefficients
#' @param alt_iter number of iterations of the alternating procedure:
#' when alt_iter == 1, this function is equivalent to asm_regression
#' @param intercept.selection mean or median of the residuals
#' if intercept.selection == "median",
#' then the standard error of the intercept estimate is set to NA
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
#'  \item{\code{info_asm}:}{antitonic information}
#'  \item{\code{I_mat}:}{estimated antitonic information matrix}
#'  \item{\code{Cov_mat}:}{asymptotic covariance matrix of the estimated coefficients}
#'  \item{\code{psi}:}{estimated antitonic score function}
#' }
#'
#' @export
#'
#' @examples
#' n <- 1000 ; d <- 2
#' X <- matrix(rnorm(n * d), n, d)
#' Y <- X %*% c(2, 3) + rnorm(n) # no intercept!
#' asm.fit(X,Y)
asm.fit <- function(X, Y, betapilot = "OLS",
                    alt_iter = 1, intercept.selection = "mean",
                    k = 3000, max_iter = 65, kernel_pts = 2^15,
                    bw = "nrd0", kernel = "gaussian", verbose=FALSE, ...) {

  THRESH = 1e-6

  X = as.matrix(X)

  has_intercept = TRUE

  ## if the first column of X is not the all-ones vector, add a column of 1s
  if (any(X[, 1] != 1)) {
    has_intercept = FALSE
    Xtilde = X

    if (is.null(colnames(X))) {
      colnames(X) = paste0("X", 1:ncol(X))
    }
  } else {
    Xtilde = as.matrix(X[, -1])

    if (is.null(colnames(X))) {
      if (ncol(X) == 1) {
        colnames(X) = c("(Intercept)")
      } else {
      colnames(X) = c("(Intercept)", paste0("X", 2:ncol(X)))
      }
    }
  }
  n = nrow(X)

  if (betapilot[1]=="OLS") {
    betapilot = solve(t(X) %*% X, t(X) %*% Y)
  } else if (betapilot[1]=="LAD") {
    LAD_fit = quantreg::rq.fit(x=X, y=Y, tau=0.5)
    betapilot = LAD_fit$coefficients
  }

  residuals = Y - X %*% betapilot

  scoreobj = Inf

  res_prev = NULL

  if (!has_intercept) {
    betapilot = c(0, betapilot)
  }

  for (it in 1:alt_iter) {
    res_asm = asm_regression(betapilot, residuals, Xtilde, Y,
                            k = k, kernel_pts = kernel_pts,
                            max_iter = max_iter,
                            bw = bw, kernel = kernel,
                            est_score_obj = TRUE, return_fn = TRUE)

    ## compare score matching objective
    scoreobj_new = res_asm$score_obj

    if (scoreobj - scoreobj_new < THRESH && it > 1) {
      if (verbose)
       cat(paste("Alternation finished at iteration", it))
      res_asm = res_prev
      break
    } else {
      ## continue to next iteration
      scoreobj = scoreobj_new
      res_prev = res_asm
    }

    thetahat_asm = res_asm$thetahat
    resid_with_mean = Y - Xtilde %*% thetahat_asm

    if (intercept.selection == "mean" && has_intercept) {
      muhat = mean(resid_with_mean)
    } else if (intercept.selection == "median" && has_intercept){
      muhat = median(resid_with_mean)
    } else
      muhat = 0

    betapilot = c(muhat, thetahat_asm)

  }

  info_asm = res_asm$info_asm
  thetahat_asm = res_asm$thetahat
  I_mat = info_asm * var(Xtilde)
  v22 = solve(I_mat)

  resid_with_mean = Y - Xtilde %*% thetahat_asm

  if (intercept.selection == "mean" && has_intercept) {
      muhat = mean(resid_with_mean)
      noise_var = var(resid_with_mean)
      Xbar <- colMeans(Xtilde)
      v11 = noise_var + Xbar %*% v22 %*% Xbar
      v12 = Xbar %*% v22
      Cov_mat = matrix(0, ncol(X), ncol(X))
      Cov_mat[1, 1] = v11
      Cov_mat[-1, -1] = v22
      Cov_mat[1, -1] = v12
      Cov_mat[-1, 1] = t(v12)
      std_errs = sqrt(diag(Cov_mat) / n)
    } else if (intercept.selection == "median" && has_intercept) {
      muhat = median(resid_with_mean)
      ## needs to be fixed
      std_errs = c(NA, sqrt(diag(v22) / n))
    } else {
      std_errs = sqrt(diag(v22) / n )
      Cov_mat = v22
    }

  if (has_intercept) {
    betahat = c(muhat, thetahat_asm)
  } else {
    betahat = thetahat_asm
  }

  betahat = as.vector(betahat)

  names(std_errs) = colnames(X)
  names(betahat) = colnames(X)

  zvals = abs(betahat) / std_errs
  sig_vals = pnorm(zvals, lower.tail = FALSE) * 2

  fitted.values = X %*% betahat
  residuals = Y - fitted.values

  model_res <- list(
    betahat = betahat,
    std_errs = std_errs,
    fitted.values = fitted.values,
    residuals = residuals,
    zvals = zvals,
    sig_vals = sig_vals,
    info_asm = info_asm,
    I_mat = I_mat,
    Cov_mat = Cov_mat,
    psi = res_asm$psi
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
  cat("Estimates:\n")
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
#' Mod <- asm(Y ~ X)
#' predict(Mod, newdata = data.frame(X = 1), interval = "prediction")
#'
predict.asm <- function(object, newdata = NULL, interval = "none",
                        level = 0.95, ...) {
  if (is.null(newdata))
    return(object$fitted.values)

  if (is.data.frame(newdata)){
    tt = terms(object$formula)
    Terms = delete.response(tt)
    m = model.frame(Terms, newdata, na.action = na.pass)
    X = model.matrix(Terms, m)
  } else {
    X = as.matrix(newdata)
  }

  if (any(X[, 1] != 1)) {
    Xint = cbind(1, X)
  } else {
    Xint = X
    X = as.matrix(X[, -1])
  }
  if (interval == "none") {
    return(Xint %*% object$betahat)
  }
  cov = object$Cov_mat

  if (ncol(cov) != ncol(Xint))
    stop("covariance matrix is not of the correct dimension")

  n = length(object$residuals)
  var_pred_mean = colSums(t(Xint) * (cov %*% t(Xint))) / n

  if (interval == "prediction") {
    var_pred_mean = var_pred_mean + as.numeric(var(object$residuals))
  }

  se = sqrt(var_pred_mean)

  a <- (1 - level)/2
  a <- c(a, 1 - a)
  fac <- qnorm(a)
  predictor = Xint %*% object$betahat
  predictor = cbind(predictor, predictor + se * fac[1],
                    predictor + se * fac[2])
  colnames(predictor) <- c("fit", "lwr", "upr")
  return(predictor)
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
