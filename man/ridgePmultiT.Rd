\name{ridgePmultiT}
\alias{ridgePmultiT}
\title{
Ridge estimation of the inverse covariance matrix with multi-target shrinkage.
}
\description{
Function that evaluates the ridge estimator of the inverse covariance matrix with multi-target shrinkage. 
}
\usage{
ridgePmultiT(S, lambda, targetList)
}
\arguments{
\item{S}{            Sample covariance \code{matrix}. }
\item{lambda}{     A \code{numeric} of positive penalty parameters. Values should be specified in the same order as the target's appearance in \code{targetList}. }
\item{targetList}{ A list of semi-positive definite target matrices towards which the precision matrix is potentially shrunken. }
}
\details{
This function generalizes the \code{\link[rags2ridges:ridgeP]{ridgeP}}-function in the sense that multiple shrinkage targets can be provided in the estimation of the precision matrix of a zero-mean multivariate normal distribution. Hence, it assumes that the data stem from \eqn{\mathcal{N}(\mathbf{0}_p, \boldsymbol{\Omega}^{-1})}. The estimator maximizes the following penalized loglikelihood:
\deqn{
\log( | \boldsymbol{\Omega} |) - \mbox{tr} ( \boldsymbol{\Omega} \mathbf{S} ) - \sum\nolimits_{g=1}^G \lambda_g \| \boldsymbol{\Omega} - \mathbf{T}_g \|_F^2,
}
where \eqn{\mathbf{S}} the sample covariance matrix, \eqn{\{ \lambda_g \}_{g=1}^G} the penalty parameters of each target matrix, and the \eqn{\{ \mathbf{T}_g \}_{g=1}^G} the precision matrix' shrinkage targets. For more details see van Wieringen \emph{et al.} (2020).
}
\value{
The function returns a regularized inverse covariance \code{matrix}.
}
\references{
van Wieringen, W.N., Stam, K.A., Peeters, C.F.W., van de Wiel, M.A. (2020), "Updating of the Gaussian graphical model through targeted penalized estimation", \emph{Journal of Multivariate Analysis}, 178, Article 104621.
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

# specify vector of penalty parameters
lambda       <- c(2, 1)

# generate precision matrix
T1       <- matrix(0.7, p, p)
diag(T1) <- 1
T2       <- diag(rep(2, p))

# generate precision matrix
Omega       <- matrix(0.4, p, p)
diag(Omega) <- 2
Sigma       <- solve(Omega)

# data 
Y <- mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=Sigma)
S <- cov(Y)

# unpenalized diagonal estimate
Phat <- ridgePmultiT(S, lambda, list(T1=T1, T2=T2))
}
