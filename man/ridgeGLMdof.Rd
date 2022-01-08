\name{ridgeGLMdof}
\alias{ridgeGLMdof}
\title{
Degrees of freedom of the generalized ridge estimator.
}
\description{
Function that evaluates the degrees of freedom of the generalized ridge estimator of the regression parameter of generalized linear models. 
}
\usage{
ridgeGLMdof(X, U=matrix(ncol=0, nrow=nrow(X)), lambda, 
         lambdaG, Dg=matrix(0, ncol=ncol(X), nrow=ncol(X)), 
         model="linear", linPred=rep(0,nrow(X)))
}
\arguments{
\item{X}{           The design \code{matrix} of the penalized covariates. The number of rows should match the number of elements of \code{Y}. }
\item{U}{           The design \code{matrix} of the unpenalized covariates. The number of rows should match the number of elements of \code{Y}. }
\item{lambda}{      A positive \code{numeric} that is the ridge penalty parameter. }
\item{lambdaG}{     A positive \code{numeric} that is the generalized ridge penalty parameter. }
\item{Dg}{          A non-negative definite \code{matrix} of the unscaled generalized ridge penalty. }
\item{model}{       A \code{character}, either \code{"linear"} and \code{"logistic"} (a reference to the models currently implemented), indicating for which generalized linear model model instance the degrees of freedom is to be evaluated. }
\item{linPred}{ A \code{numeric}, the linear predictor associated with the provided \code{X}, \code{U}, \code{lambda}, \code{lambdaG}, and \code{Dg}. The number of elements of \code{linPred} should match the number of rows of \code{X} and \code{U}. }
}
\details{
The degrees of freedom of the regular ridge regression estimator is usually defined the trace of the ridge hat matrix: \eqn{\mbox{tr} [ \mathbf{X} (\mathbf{X}^{\top} \mathbf{X} + \lambda \mathbf{I}_{pp})^{-1} \mathbf{X}^{\top}]}. That of the regular ridge logistic regression estimator is defined analoguously by Park, Hastie (2008). Lettink et al. (2022) translates these definitions to the generalized ridge (logistic) regression case.
}
\value{
A \code{numeric}, the degrees of freedom consumed by the (generalized) ridge (logistic) regression estimator.
} 
\references{
Park, M. Y., & Hastie, T. (2008). Penalized logistic regression for detecting gene interactions. \emph{Biostatistics}, 9(1), 30-50.

Lettink, A., Chinapaw, M.J.M., van Wieringen, W.N. (2022), "Two-dimensional fused targeted ridge regression for health indicator prediction from accelerometer data", \emph{submitted}.
}
\author{
W.N. van Wieringen.
}
\examples{
# set the sample size
n <- 50

# set the true parameter
betas <- (c(0:100) - 50) / 20

# generate covariate data
X <- matrix(rnorm(length(betas)*n), nrow=n)

# set the penalty parameter
lambda <- 3

# estimate the logistic regression parameter
dofs <- ridgeGLMdof(X, lambda=lambda, lambdaG=0, 
                    model="logistic", 
                    linPred=tcrossprod(X, t(betas)))
}
