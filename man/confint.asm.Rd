% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/asm.R
\name{confint.asm}
\alias{confint.asm}
\title{Confidence intervals for coefficients
in an \code{asm} regression model}
\usage{
\method{confint}{asm}(object, parm, level = 0.95, ...)
}
\arguments{
\item{object}{asm object}

\item{parm}{parameters to calculate confidence intervals}

\item{level}{confidence level}

\item{...}{additional arguments to ensure compatibility with the generic function confint()}
}
\value{
matrix of confidence intervals for the regression coefficients
}
\description{
Computes confidence intervals for individual regression coefficients based on a fitted \code{asm} regression model
}
\examples{
model = asm(mpg ~ cyl + hp + disp, data=mtcars)
confint(model)

}
