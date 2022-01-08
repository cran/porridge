\name{ridgeGLM}
\alias{ridgeGLM}
\title{
Ridge estimation of generalized linear models.
}
\description{
Function that evaluates the targeted ridge estimator of the regression parameter of generalized linear models. 
}
\usage{
ridgeGLM(Y, X, U=matrix(ncol=0, nrow=length(Y)), lambda, 
         lambdaG=0, Dg=matrix(0, ncol=ncol(X), nrow=ncol(X)), 
         target=rep(0, ncol(X)), model="linear", 
         minSuccDiff=10^(-10), maxIter=100)
}
\arguments{
\item{Y}{           A \code{numeric} being the response vector. }
\item{X}{           The design \code{matrix} of the penalized covariates. The number of rows should match the number of elements of \code{Y}. }
\item{U}{           The design \code{matrix} of the unpenalized covariates. The number of rows should match the number of elements of \code{Y}. }
\item{lambda}{      A positive \code{numeric} that is the ridge penalty parameter. }
\item{lambdaG}{     A positive \code{numeric} that is the generalized ridge penalty parameter. }
\item{Dg}{          A non-negative definite \code{matrix} of the unscaled generalized ridge penalty. }
\item{target}{      A \code{numeric} towards which the estimate is shrunken. }
\item{model}{       A \code{character}, either \code{"linear"} and \code{"logistic"} (a reference to the models currently implemented), indicating which generalized linear model model instance is to be fitted. }
\item{minSuccDiff}{ A \code{numeric}, the minimum distance between the loglikelihoods of two successive iterations to be achieved. Used only if \code{model="logistic"}. }
\item{maxIter}{     A \code{numeric} specifying the maximum number of iterations. Used only if \code{model="logistic"}. }
}
\details{
This function finds the maximizer of the following penalized loglikelihood: \eqn{ \mathcal{L}( \mathbf{Y}, \mathbf{X}, \mathbf{U}; \boldsymbol{\beta}, \boldsymbol{\gamma}) - \frac{1}{2} \lambda \| \boldsymbol{\beta}  - \boldsymbol{\beta}_0 \|_2^2 - \frac{1}{2} \lambda_g ( \boldsymbol{\beta}  - \boldsymbol{\beta}_0 )^\top \Delta_g ( \boldsymbol{\beta}  - \boldsymbol{\beta}_0 )}, with loglikelihood \eqn{\mathcal{L}( \mathbf{Y}, \mathbf{X}; \boldsymbol{\beta})}, response \eqn{\mathbf{Y}}, design matrices \eqn{\mathbf{X}} and \eqn{\mathbf{U}}, regression parameters \eqn{\boldsymbol{\beta}} and \eqn{\boldsymbol{\gamma}}, penalty parameter \eqn{\lambda}, shrinkage target \eqn{\boldsymbol{\beta}_0}, and generalized ridge penalty matrix \eqn{\Delta_g}. For more details, see van Wieringen, Binder (2020) and Lettink et al. (2022).
}
\note{
The penalized IRLS (Iterative Reweighted Least Squares) algorithm for the evaluation of the generalized ridge logistic regression estimator may fail to converge for small penalty parameter values in combination with a nonzero shrinkage target. 
} 
\value{
A \code{numeric}, the generalized ridge estimate of the regression parameter. If a nonempty \eqn{\mathbf{U}} is supplied, the first few elements are the unpenalized effect estimates of the covariates that comprise this design matrix.
} 
\references{
van Wieringen, W.N. Binder, H. (2020), "Online learning of regression models from a sequence of datasets by penalized estimation", \emph{submitted}.

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

# sample the response
probs <- exp(tcrossprod(betas, X)[1,]) / (1 + exp(tcrossprod(betas, X)[1,]))
Y     <- numeric()
for (i in 1:n){
    Y <- c(Y, sample(c(0,1), 1, prob=c(1-probs[i], probs[i])))
}

# set the penalty parameter
lambda <- 3

# estimate the logistic regression parameter
bHat <- ridgeGLM(Y, X, lambda=lambda, target=betas/2, model="logistic")
}
