\name{ridgeGLM}
\alias{ridgeGLM}
\title{
Ridge estimation of generalized linear models.
}
\description{
Function that evaluates the targeted ridge estimator of the regression parameter of generalized linear models. 
}
\usage{
ridgeGLM(Y, X, lambda, target=rep(0, ncol(X)), model, 
         minSuccDiff=10^(-10), maxIter=100)
}
\arguments{
\item{Y}{           A \code{numeric} being the response vector. }
\item{X}{           The design \code{matrix}. The number of rows should match the number of elements of \code{Y}. }
\item{lambda}{      A positive \code{numeric} that is the ridge penalty parameter. }
\item{target}{      A \code{numeric} towards which the estimate is shrunken. }
\item{model}{       A \code{character}, either \code{"linear"} and \code{"logistic"} (a reference to the models currently implemented), indicating which generalized linear model model instance is to be fitted. }
\item{minSuccDiff}{ A \code{numeric}, the minimum distance between the loglikelihoods of two successive iterations to be achieved. Used only if \code{model="logistic"}. }
\item{maxIter}{     A \code{numeric} specifying the maximum number of iterations. Used only if \code{model="logistic"}. }
}
\details{
This function finds the maximizer of the following penalized loglikelihood: \eqn{ \mathcal{L}( \mathbf{Y}, \mathbf{X}; \boldsymbol{\beta}) - \frac{1}{2} \lambda \| \boldsymbol{\beta}  - \boldsymbol{\beta}_0 \|_2^2}, with loglikelihood \eqn{\mathcal{L}( \mathbf{Y}, \mathbf{X}; \boldsymbol{\beta})}, response \eqn{\mathbf{Y}}, design matrix \eqn{\mathbf{X}}, regression parameter \eqn{\boldsymbol{\beta}}, penalty parameter \eqn{\lambda}, and shrinkage target \eqn{\boldsymbol{\beta}_0}. For more details, see van Wieringen, Binder (2020).
}
\value{
The ridge estimate of the regression parameter.
} 
\references{
van Wieringen, W.N. Binder, H. (2020), "Online learning of regression models from a sequence of datasets by penalized estimation", \emph{submitted}.
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
bHat <- ridgeGLM(Y, X, lambda, target=betas/2, model="logistic")
}
