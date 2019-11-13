\name{ridgeGGMmixture}
\alias{ridgeGGMmixture}
\title{
Ridge penalized estimation of a mixture of GGMs.
}
\description{
Function that estimates a mixture of GGMs (Gaussian graphical models) through a ridge penalized EM (Expectation-Maximization) algorithm as described in Aflakparast \emph{et al}. (2018).
}
\usage{
ridgeGGMmixture(Y, K, lambda, target,                                    
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
\item{lambda}{ A positive \code{numeric} representing the ridge penalty parameter. }
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
The function returns a regularized inverse covariance \code{list}-object with slots:
\item{mu}{ A \code{matrix} with estimated mean vectors are rows. }
\item{P}{ A \code{matrix} with estimated mixture precision matrices stacked on top of each other. }
\item{pi}{ A \code{numeric} wth estimated mixing probabilities.}
\item{weights}{ A \code{matrix} wth estimated component memberships.}
\item{penLL}{ A \code{numeric} with the penalized loglikelihood of the estimated model. }
}

\details{
The data are assumed to follow a mixture of \eqn{K} Gaussian graphical models:
\deqn{ \mathbf{Y}_i \sim \sum\nolimits_{k=1}^K \theta_k \mathcal{N}(\boldsymbol{\mu}_k, \boldsymbol{\Omega}_k^{-1}), }
where \eqn{\theta_k = P(Z_i = k)} is the probability that the \eqn{i}-th sample stems from the \eqn{k}-the component. The model parameters are estimated by ridge penalized likelihood maximization:
\deqn{ \sum\nolimits_{i=1}^n \log [ \sum\nolimits_{k=1}^K \theta_k P(\mathbf{Y}_i \, | \, Z_i = k; \boldsymbol{\mu}_k, \boldsymbol{\Omega}_k) ] + \lambda \sum\nolimits_{k=1}^K \| \boldsymbol{\Omega}_k - \mathbf{T}_k \|_F^2, }
where \eqn{\lambda} is the penalty parameter and \eqn{\mathbf{T}_k} is the shrinkage target of the \eqn{k}-th component's precision matrix. This function yields the maximizer of this penalized loglikelihood, which is found by means of a penalized EM algorithm.
}
\author{
W.N. van Wieringen, M. Aflakparast.
}
\references{
Aflakparast, M., de Gunst, M.C.M., van Wieringen, W.N. (2018), "Reconstruction of molecular network evolution from cross-sectional omics data", \emph{Biometrical Journal}, 60(3), 547-563.
}
\seealso{
\code{optPenaltyGGMmixture.kCVauto}
}
\examples{
# define mixing proportions
pis <- c(0.2, 0.3, 0.4)

# set dimension and sample size
p <- 5
n <- 100

# define population covariance matrices
diags       <- list(rep(1, p), 
                    rep(0.5, p-1), 
                    rep(0.25, p-2), 
                    rep(0.1, p-3))
Omega       <- as.matrix(Matrix::bandSparse(p, 
                                            k=-c(0:3), 
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
\donttest{optLambda <- optPenaltyGGMmixture.kCVauto(Y,  K,            }
\donttest{                                          0.00001, 100,     }
\donttest{                                          10, fold=5,       }
\donttest{                                          target=0*Sigma1)  }

# ridge penalized estimation of the GGM mixture
ridgeGGMmixFit <- ridgeGGMmixture(Y, 3, 1, target=0*Sigma1)
}
