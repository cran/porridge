\name{ridgePgen.kCV}
\alias{ridgePgen.kCV}
\title{
K-fold cross-validated loglikelihood of ridge precision estimator.
}
\description{
Function that calculates of the k-fold cross-validated negative (!) loglikelihood of the generalized ridge precision estimator.
}
\usage{
ridgePgen.kCV(lambda, Y, fold=nrow(Y), target, 
              nInit=100, minSuccDiff=10^(-5)) 
}
\arguments{
\item{lambda}{ A symmetric \code{matrix} with element-wise positive penalty parameters. }
\item{Y}{ Data \code{matrix} with samples as rows and variates as columns. }
\item{fold}{ A \code{numeric} or \code{integer} specifying the number of folds to apply in the cross-validation. }
\item{target}{ A semi-positive definite target \code{matrix} towards which the estimate is shrunken. }
\item{nInit}{ A \code{numeric} specifying the number of iterations. }
\item{minSuccDiff}{ A \code{numeric}: minimum successive difference (in terms of their penalized loglikelihood) between two succesive estimates to be achieved. }
}
\value{
The function returns a \code{numeric} containing the cross-validated negative loglikelihood.
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
