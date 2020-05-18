\name{optPenaltyPgen.kCVauto.banded}
\alias{optPenaltyPgen.kCVauto.banded}
\title{
Automatic search for optimal penalty parameter (generalized ridge precision).
}
\description{
Function that determines the optimal penalty parameters through maximization of the k-fold cross-validated log-likelihood score, with a penalization that encourages banded precisions.
}
\usage{
optPenaltyPgen.kCVauto.banded(Y, lambdaMin, lambdaMax, 
                          lambdaInit=(lambdaMin + lambdaMax)/2,
                          fold=nrow(Y), target, 
                          zeros=matrix(nrow=0, ncol=2), 
                          penalize.diag=TRUE, nInit=100, 
                          minSuccDiff=10^(-5)) 
}
\arguments{
\item{Y}{ Data \code{matrix} with samples as rows and variates as columns. }
\item{lambdaMin}{ A \code{numeric} giving the minimum value for the penalty parameters. One value per group. Values should be specified in the same order as the first appearance of a group representative. }
\item{lambdaMax}{ A \code{numeric} giving the maximum value for the penalty parameters. One value per group. Values should be specified in the same order as the first appearance of a group representative. }
\item{lambdaInit}{ A \code{numeric} giving the initial (starting) value for the penalty parameters. One value per group. Values should be specified in the same order as the first appearance of a group representative. }
\item{fold}{ A \code{numeric} or \code{integer} specifying the number of folds to apply in the cross-validation. }
\item{target}{ A semi-positive definite target \code{matrix} towards which the estimate is shrunken. }
\item{zeros}{ A two-column \code{matrix}, with the first and second column containing the row- and column-index of zero precision elements. }
\item{penalize.diag}{ A \code{logical} indicating whether the diagonal should be penalized. }    
\item{nInit}{ A \code{numeric} specifying the number of iterations. }
\item{minSuccDiff}{ A \code{numeric}: minimum successive difference (in terms of their penalized loglikelihood) between two succesive estimates to be achieved. }
}
\details{
The penalty matrix \eqn{\boldsymbol{\Lambda}} is parametrized as follows. The elements of \eqn{\boldsymbol{\Lambda}} are \eqn{(\boldsymbol{\Lambda})_{j,j'} = \lambda (| j - j'| + 1)} for 
\eqn{j, j' = 1, \ldots, p}.
}
\value{
The function returns a \code{numeric} containing the cross-validated optimal positive penalty parameters.
}
\references{
van Wieringen, W.N. (2019), "The generalized ridge estimator of the inverse covariance matrix", \emph{Journal of Computational and Graphical Statistics}, 28(4), 932-942.
}
\author{
W.N. van Wieringen.
}
\seealso{
\code{\link{ridgePgen}}
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
lambdaOpt <- optPenaltyPgen.kCVauto.banded(Y, 10^(-10), 10^(10),
                          target=matrix(0, p, p),
                          penalize.diag=FALSE, nInit=100, 
                          minSuccDiff=10^(-5)) 

# format the penalty matrix
lambdaOptMat <- matrix(NA, p, p)
for (j1 in 1:p){
    for (j2 in 1:p){
        lambdaOptMat[j1, j2] <- lambdaOpt * (abs(j1-j2)+1)
    }
}

# generalized ridge precision estimate
Phat <- ridgePgen(S, lambdaOptMat, matrix(0, p, p))
}
