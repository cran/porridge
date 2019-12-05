\name{optPenaltyGGMmixture.kCVauto}
\alias{optPenaltyGGMmixture.kCVauto}
\title{
Automatic search for optimal penalty parameter (mixture of GGMs).
}
\description{
Function that performs an automatic search for the optimal penalty parameter for the \code{ridgeGGMmixture} call by employing Brent's
method to the calculation of a cross-validated (negative) log-likelihood score.
}
\usage{
optPenaltyGGMmixture.kCVauto(Y, K, lambdaMin, lambdaMax, 
                lambdaInit=(lambdaMin+lambdaMax)/2, 
                fold=nrow(Y), target,               
                iWeights=matrix(sample(seq(0+1/nrow(Y), 
                                1-1/nrow(Y), by=1/(2*nrow(Y))), 
                                nrow(Y)*K, replace=TRUE), 
                                nrow=nrow(Y), ncol=K),
                nInit=100, minSuccDiff=10^(-10),
                minMixProp=0.01)
}
\arguments{
\item{Y}{ Data \code{matrix} with samples as rows and variates as columns. }
\item{K}{ A \code{numeric}, specifying the number of mixture components.}
\item{lambdaMin}{ A \code{numeric} giving the minimum value for the penalty parameter. }
\item{lambdaMax}{ A \code{numeric} giving the maximum value for the penalty parameter. }
\item{lambdaInit}{ A \code{numeric} giving the initial (starting) value for the penalty parameter. }
\item{fold}{ A \code{numeric} or \code{integer} specifying the number of folds to apply in the cross-validation. }
\item{target}{ A semi-positive definite target \code{matrix} towards which the estimate is shrunken. }
\item{iWeights}{ Sample-specific positive component weight \code{matrix}. Rows correspond to samples, while columns to components. } 
\item{nInit}{ A \code{numeric} specifying the number of iterations. }
\item{minSuccDiff}{ A \code{numeric}: minimum successive difference (in terms of their penalized loglikelihood) between two succesive estimates to be achieved. }
\item{minMixProp}{ Smallest mixing probability tolerated. }
}
\note{
The elements of \code{iWeights} may be larger than one as they are rescaled internally to sum to one.
}
\value{
The function returns a positive \code{numeric}, the cross-validated optimal penalty parameter.
}
\references{
Aflakparast, M., de Gunst, M.C.M., van Wieringen, W.N. (2018), "Reconstruction of molecular network evolution from cross-sectional omics data", \emph{Biometrical Journal}, 60(3), 547-563.
}
\author{
W.N. van Wieringen, M. Aflakparast.
}
\seealso{
\code{ridgeGGMmixture}
}
\examples{
# define mixing proportions
pis <- c(0.2, 0.3, 0.4)

# set dimension and sample size
p <- 5
n <- 100

# define population covariance matrices
diags       <- list(rep(1,    p), 
                    rep(0.5,  p-1), 
                    rep(0.25, p-2), 
                    rep(0.1,  p-3))
Omega       <- as.matrix(Matrix::bandSparse(p, 
                                            k   =-c(0:3), 
                                            diag=c(diags), 
                                            symm=TRUE))
Sigma1      <- solve(Omega)
Omega       <- matrix(0.3, p, p)
diag(Omega) <- 1
Sigma2      <- solve(Omega)
Sigma3      <- cov(matrix(rnorm(p*n), ncol=p))

# mean vectors
mean1 <- rep(0,p)
mean2 <- rexp(p)
mean3 <- rnorm(p)

# draw data data from GGM mixture
Z <- sort(sample(c(1:3), n, prob=pis, replace=TRUE))
Y <- rbind(mvtnorm::rmvnorm(sum(Z==1), mean=mean1, sigma=Sigma1),
           mvtnorm::rmvnorm(sum(Z==2), mean=mean2, sigma=Sigma2),
           mvtnorm::rmvnorm(sum(Z==3), mean=mean3, sigma=Sigma3))

# find optimal penalty parameter
\donttest{optLambda <- optPenaltyGGMmixture.kCVauto(Y,  K,           }
\donttest{                                          0.00001, 100,    }
\donttest{                                          10, fold=5,      }
\donttest{                                          target=0*Sigma1) }

# ridge penalized estimation of the GGM mixture
\donttest{ridgeGGMmixFit <- ridgeGGMmixture(Y, 3, optLambda, target=0*Sigma1) }
}
