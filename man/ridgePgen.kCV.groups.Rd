\name{ridgePgen.kCV.groups}
\alias{ridgePgen.kCV.groups}
\title{
K-fold cross-validated loglikelihood of ridge precision estimator with group-wise penalized variates.
}
\description{
Function that calculates of the k-fold cross-validated negative (!) loglikelihood of the generalized ridge precision estimator, assuming that variates are grouped and penalized group-wise.
}
\usage{
ridgePgen.kCV.groups(lambdaGrps, Y, fold=nrow(Y), 
                     groups, target, 
                     zeros=matrix(nrow=0, ncol=2), 
                     penalize.diag=TRUE, nInit=100, 
                     minSuccDiff=10^(-5)) 
}
\arguments{
\item{lambdaGrps}{ A \code{numeric} with penalty parameter values, one per group. Values should be specified in the same order as the first appearance of a group representative. }
\item{Y}{ Data \code{matrix} with samples as rows and variates as columns. }
\item{fold}{ A \code{numeric} or \code{integer} specifying the number of folds to apply in the cross-validation. }
\item{groups}{ A\code{numeric} indicating to which group a variate belongs. Same values indicate same group. }
\item{target}{ A semi-positive definite target \code{matrix} towards which the estimate is shrunken. }
\item{zeros}{ A two-column \code{matrix}, with the first and second column containing the row- and column-index of zero precision elements. }
\item{penalize.diag}{ A \code{logical} indicating whether the diagonal should be penalized. }    
\item{nInit}{ A \code{numeric} specifying the number of iterations. }
\item{minSuccDiff}{ A \code{numeric}: minimum successive difference (in terms of their penalized loglikelihood) between two succesive estimates to be achieved. }
}
\value{
The function returns a \code{numeric} containing the cross-validated negative loglikelihood.
}
\details{
The penalty matrix \eqn{\boldsymbol{\Lambda}} is parametrized as follows. The elements of \eqn{\boldsymbol{\Lambda}} are \eqn{(\boldsymbol{\Lambda})_{j,j'} = \frac{1}{2} (\lambda_k + \lambda_{k'})} for 
\eqn{j, j' = 1, \ldots, p} if \eqn{j} and \eqn{j'} belong to groups \eqn{k} and \eqn{k'}, respectively, where \eqn{\lambda_k} and \eqn{\lambda_{k'}} are the corresponding group-specific penalty parameters.
}
\author{
W.N. van Wieringen.
}
\references{
van Wieringen, W.N. (2019), "The generalized ridge estimator of the inverse covariance matrix", \emph{Journal of Computational and Graphical Statistics}, 28(4), 932-942.
}
\seealso{
\code{ridgePgen}
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

# find optimal penalty parameters through cross-validation
lambdaOpt <- optPenaltyPgen.kCVauto.groups(Y, rep(10^(-10), 2), rep(10^(10), 2), 
                          groups=c(rep(0, p/2), rep(1, p/2)), 
                          target=matrix(0, p, p),
                          penalize.diag=FALSE, nInit=100, 
                          minSuccDiff=10^(-5)) 

# format the penalty matrix
lambdaOptVec <- c(rep(lambdaOpt[1], p/2), rep(lambdaOpt[2], p/2))
lambdaOptMat <- outer(lambdaOptVec, lambdaOptVec, "+")

# generalized ridge precision estimate
Phat <- ridgePgen(S, lambdaOptMat, matrix(0, p, p))
}
