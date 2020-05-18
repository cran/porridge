\name{ridgePgen}
\alias{ridgePgen}
\title{
Ridge estimation of the inverse covariance matrix with element-wise penalization and shrinkage.
}
\description{
Function that evaluates the generalized ridge estimator of the inverse covariance matrix with element-wise penalization and shrinkage. 
}
\usage{
ridgePgen(S, lambda, target, nInit=100, minSuccDiff=10^(-10))
}
\arguments{
\item{S}{           Sample covariance \code{matrix}. }
\item{lambda}{      A symmetric \code{matrix} with element-wise positive penalty parameters. }
\item{target}{      A semi-positive definite target \code{matrix} towards which the estimate is shrunken. }
\item{nInit}{       A \code{numeric} specifying the number of iteration. }
\item{minSuccDiff}{ A \code{numeric}: minimum distance between two succesive estimates to be achieved. }
}
\details{
This function generalizes the \code{\link[rags2ridges:ridgeP]{ridgeP}}-function in the sense that, besides element-wise shrinkage, it allows for element-wise penalization in the estimation of the precision matrix of a zero-mean multivariate normal distribution. Hence, it assumes that the data stem from \eqn{\mathcal{N}(\mathbf{0}_p, \boldsymbol{\Omega}^{-1})}. The estimator maximizes the following penalized loglikelihood:
\deqn{
\log( | \boldsymbol{\Omega} |) - \mbox{tr} ( \boldsymbol{\Omega} \mathbf{S} ) - \| \boldsymbol{\Lambda} \circ (\boldsymbol{\Omega} - \mathbf{T}) \|_F^2,
}
where \eqn{\mathbf{S}} the sample covariance matrix, \eqn{\boldsymbol{\Lambda}} a symmetric, positive matrix of penalty parameters, the \eqn{\circ}-operator represents the Hadamard or element-wise multipication, and \eqn{\mathbf{T}} the precision matrix' shrinkage target. For more details see van Wieringen (2019).
}
\value{
The function returns a regularized inverse covariance \code{matrix}.
}
\references{
van Wieringen, W.N. (2019), "The generalized ridge estimator of the inverse covariance matrix", \emph{Journal of Computational and Graphical Statistics}, 28(4), 932-942.
}
\author{
W.N. van Wieringen.
}
\seealso{
\code{\link[rags2ridges:ridgeP]{ridgeP}}, \code{\link[rags2ridges:ridgePchordal]{ridgePchordal}}.
}
\examples{
# set dimension and sample size
p <- 10
n <- 10

# penalty parameter matrix
lambda       <- matrix(1, p, p)
diag(lambda) <- 0.1

# generate precision matrix
Omega       <- matrix(0.4, p, p)
diag(Omega) <- 1
Sigma       <- solve(Omega)

# data 
Y <- mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=Sigma)
S <- cov(Y)

# unpenalized diagonal estimate
Phat <- ridgePgen(S, lambda, 0*S)
}
