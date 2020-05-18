\name{optPenaltyPmultiT.kCVauto}
\alias{optPenaltyPmultiT.kCVauto}
\title{
Automatic search for optimal penalty parameter (ridge precision with multi-targets).
}
\description{
Function that determines the optimal penalty parameters through maximization of the k-fold cross-validated log-likelihood score, assuming that variates are grouped and penalized group-wise.
}
\usage{
optPenaltyPmultiT.kCVauto(Y, lambdaMin, lambdaMax, 
                          lambdaInit=(lambdaMin+lambdaMax)/2,
                          fold=nrow(Y), targetList) 
}
\arguments{
\item{Y}{ Data \code{matrix} with samples as rows and variates as columns. }
\item{lambdaMin}{ A \code{numeric} giving the minimum value for the penalty parameters. One value per target. Values should be specified in the same order as the target's appearance in \code{targetList}. }
\item{lambdaMax}{ A \code{numeric} giving the maximum value for the penalty parameters. One value per group. Values should be specified in the same order as the target's appearance in \code{targetList}. }
\item{lambdaInit}{ A \code{numeric} giving the initial (starting) value for the penalty parameters. One value per group. Values should be specified in the same order as the target's appearance in \code{targetList}. }
\item{fold}{ A \code{numeric} or \code{integer} specifying the number of folds to apply in the cross-validation. }
\item{targetList}{ A list of semi-positive definite target matrices towards which the precision matrix is potentially shrunken. }
}
\value{
The function returns a \code{numeric} containing the cross-validated optimal positive penalty parameters.
}
\author{
W.N. van Wieringen.
}
\references{
van Wieringen, W.N., Stam, K.A., Peeters, C.F.W., van de Wiel, M.A. (2020), "Updating of the Gaussian graphical model through targeted penalized estimation", \emph{Journal of Multivariate Analysis}, 178, Article 104621.
}
\seealso{
\code{ridgePmultiT}
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

# find optimal penalty parameters through cross-validation
lambdaOpt <- optPenaltyPmultiT.kCVauto(Y, rep(10^(-10), 2), 
                          rep(10^(10), 2), rep(1, 2), 
                          targetList=list(T1=T1, T2=T2)) 

# unpenalized diagonal estimate
Phat <- ridgePmultiT(S, lambdaOpt, list(T1=T1, T2=T2))
}
