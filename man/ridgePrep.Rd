\name{ridgePrep}
\alias{ridgePrep}
\title{
Ridge penalized estimation of the precision matrix from data with replicates.
}
\description{
Estimation of the precision matrix from data with replicates through a ridge penalized EM (Expectation-Maximization) algorithm. It assumes a simple 'signal+noise' model, both random variables are assumed to be drawn from a multivariate normal distribution with their own unstructured precision matrix. These precision matrices are estimated.
}
\usage{
ridgePrep(Y, ids, lambdaZ, lambdaE, 		
          targetZ=matrix(0, ncol(Y), ncol(Y)),
          targetE=matrix(0, ncol(Y), ncol(Y)),
          nInit=100, minSuccDiff=10^(-10))
}
\arguments{
\item{Y}{ Data \code{matrix} with samples (including the repetitions) as rows and variates as columns. }
\item{ids}{ A \code{numeric} indicating which rows of \code{Y} belong to the same individal.}
\item{lambdaZ}{ A positive \code{numeric} representing the ridge penalty parameter for the signal precision matrix estimate. }
\item{lambdaE}{ A positive \code{numeric} representing the ridge penalty parameter for the error precision matrix estimate. }
\item{targetZ}{ A semi-positive definite target \code{matrix} towards which the signal precision matrix estimate is shrunken. }
\item{targetE}{ A semi-positive definite target \code{matrix} towards which the error precision matrix estimate is shrunken. }
\item{nInit}{ A \code{numeric} specifying the number of iterations. }
\item{minSuccDiff}{ A \code{numeric}: minimum successive difference (in terms of the relative change in the absolute difference of the penalized loglikelihood) between two succesive estimates to be achieved. }
}
\value{
The function returns the regularized inverse covariance \code{list}-object with slots:
\item{Pz}{ The estimated signal precision matrix. }
\item{Pz}{ The estimated error precision matrix. }
\item{penLL}{ The penalized loglikelihood of the estimated model. }
}
\details{
Data are assumed to originate from a design with replicates. Each observation \eqn{\mathbf{Y}_{i,k_i}} with \eqn{k_i} (\eqn{k_i = 1, \ldots, K_i}) the \eqn{k_i}-th replicate of the \eqn{i}-th sample, is described by a `signal+noise' model: \eqn{\mathbf{Y}_{i,k_i} = \mathbf{Z}_i + \boldsymbol{\varepsilon}_{i,k_i}}, where \eqn{\mathbf{Z}_i} and \eqn{\boldsymbol{\varepsilon}_{i,k_i}} represent the signal and noise, respectively. Each observation \eqn{\mathbf{Y}_{i,k_i}} follows a multivariate normal law of the form
\eqn{\mathbf{Y}_{i,k_i} \sim \mathcal{N}(\mathbf{0}_p, \boldsymbol{\Omega}_z^{-1} + \boldsymbol{\Omega}_{\varepsilon}^{-1})}, which results from the distributional assumptions of the signal and the noise, \eqn{\mathbf{Z}_{i} \sim \mathcal{N}(\mathbf{0}_p, \boldsymbol{\Omega}_z^{-1})} and \eqn{\boldsymbol{\varepsilon}_{i, k_i} \sim \mathcal{N}(\mathbf{0}_p, \boldsymbol{\Omega}_{\varepsilon}^{-1})}, and their independence. The model parameters are estimated by means of a penalized EM algorithm that maximizes the loglikelihood augmented with the penalty \eqn{\lambda_z \| \boldsymbol{\Omega}_z - \mathbf{T}_z \|_F^2 + \lambda_{\varepsilon} \| \boldsymbol{\Omega}_{\varepsilon} - \mathbf{T}_{\varepsilon} \|_F^2}, in which \eqn{\mathbf{T}_z} and \eqn{\mathbf{T}_{\varepsilon}} are the shrinkage targets of the signal and noise precision matrices, respectively. For more details see van Wieringen and Chen (2019).
}
\author{
W.N. van Wieringen.
}
\references{
van Wieringen, W.N., Chen, Y. (2019), "Penalized estimation of the Gaussian graphical model from data with replicates", \emph{submitted}.
}
\seealso{
\code{optPenaltyPrep.kCVauto}
}
\examples{
# set parameters
p        <- 10
Se       <- diag(runif(p))
Sz       <- matrix(3, p, p)
diag(Sz) <- 4

# draw data
n <- 100
ids <- numeric()
Y   <- numeric()
for (i in 1:n){
     Ki <- sample(2:5, 1)
     Zi <- mvtnorm::rmvnorm(1, sigma=Sz)
     for (k in 1:Ki){
          Y   <- rbind(Y, Zi + mvtnorm::rmvnorm(1, sigma=Se))
          ids <- c(ids, i)
     }
}

# estimate
Ps <- ridgePrep(Y, ids, 1, 1) 
}
