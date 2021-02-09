\name{ridgeGLMmultiT}
\alias{ridgeGLMmultiT}
\title{
Multi-targeted ridge estimation of generalized linear models.
}
\description{
Function that evaluates the multi-targeted ridge estimator of the regression parameter of generalized linear models. 
}
\usage{
ridgeGLMmultiT(Y, X, lambdas, targetMat, model, 
               minSuccDiff=10^(-10), maxIter=100)
}
\arguments{
\item{Y}{           A \code{numeric} being the response vector. }
\item{X}{           The design \code{matrix}. The number of rows should match the number of elements of \code{Y}. }
\item{lambdas}{     A positive \code{numeric}, vector of penalty parameters, one per target. }
\item{targetMat}{   A \code{matrix} with targets for the regression parameter as columns. }
\item{model}{       A \code{character}, either \code{"linear"} and \code{"logistic"} (a reference to the models currently implemented), indicating which generalized linear model model instance is to be fitted. }
\item{minSuccDiff}{ A \code{numeric}, the minimum distance between the loglikelihoods of two successive iterations to be achieved. Used only if \code{model="logistic"}. }
\item{maxIter}{     A \code{numeric} specifying the maximum number of iterations. Used only if \code{model="logistic"}.}
}
\details{
This function finds the maximizer of the following penalized loglikelihood: \eqn{ \mathcal{L}( \mathbf{Y}, \mathbf{X}; \boldsymbol{\beta}) - \frac{1}{2} \sum_{k=1}^K \lambda_k \| \boldsymbol{\beta}  - \boldsymbol{\beta}_{k,0} \|_2^2}, with loglikelihood \eqn{\mathcal{L}( \mathbf{Y}, \mathbf{X}; \boldsymbol{\beta})}, response \eqn{\mathbf{Y}}, design matrix \eqn{\mathbf{X}}, regression parameter \eqn{\boldsymbol{\beta}}, penalty parameter \eqn{\lambda}, and the \eqn{k}-th shrinkage target \eqn{\boldsymbol{\beta}_{k,0}}. For more details, see van Wieringen, Binder (2020).
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
lambdas <- c(1,3)

# estimate the logistic regression parameter
bHat <- ridgeGLMmultiT(Y, X, lambdas, model="logistic",
                       targetMat=cbind(betas/2, rnorm(length(betas))))
}
