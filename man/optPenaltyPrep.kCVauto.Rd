\name{optPenaltyPrep.kCVauto}
\alias{optPenaltyPrep.kCVauto}
\title{
Automatic search for optimal penalty parameters (for precision estimation of data with replicates).
}
\description{
Function that performs an automatic search of the optimal penalty parameter for the \code{ridgePrep} call by employing either the Nelder-Mead or quasi-Newton 
method to calculate the cross-validated (negative) log-likelihood score.
}
\usage{
optPenaltyPrep.kCVauto(Y, ids, lambdaInit, 
                       fold=nrow(Y), CVcrit, 
                       splitting="stratified",
                       targetZ=matrix(0, ncol(Y), ncol(Y)),
                       targetE=matrix(0, ncol(Y), ncol(Y)),
                       nInit=100, minSuccDiff=10^(-10))

}
\arguments{
\item{Y}{ Data \code{matrix} with samples (including the repetitions) as rows and variates as columns. }
\item{ids}{ A \code{numeric} indicating which rows of \code{Y} belong to the same individal.}
\item{lambdaInit}{ A \code{numeric} giving the initial (starting) values for the two penalty parameters. }
\item{fold}{ A \code{numeric} or \code{integer} specifying the number of folds to apply in the cross-validation. }
\item{CVcrit}{ A \code{character} with the cross-validation criterion to applied. Either \code{CVcrit="LL"} (the loglikelihood) or \code{CVcrit="Qloss"} (the quadratic loss). }
\item{splitting}{ A \code{character}, either \code{splitting="replications"}, \code{splitting="samples"}, or \code{splitting="stratified"}, specifying either how the splits are to be formed: either  replications or samples are randomly divided over the \code{fold} splits (first two options, respectively), or samples are randomly divided over the \code{fold} splits but in a stratified manner such that the total number of replicates in each group is roughly comparable. }
\item{targetZ}{ A semi-positive definite target \code{matrix} towards which the signal precision matrix estimate is shrunken. }
\item{targetE}{ A semi-positive definite target \code{matrix} towards which the error precision matrix estimate is shrunken. }
\item{nInit}{ A \code{numeric} specifying the number of iterations. }
\item{minSuccDiff}{ A \code{numeric}: minimum successive difference (in terms of the relative change in the absolute difference of the penalized loglikelihood) between two succesive estimates to be achieved. }
}
\value{
The function returns a all-positive \code{numeric}, the cross-validated optimal penalty parameters.
}
\references{
van Wieringen, W.N., Chen, Y. (2019), "Penalized estimation of the Gaussian graphical model from data with replicates", \emph{submitted}.
}
\author{
W.N. van Wieringen.
}
\seealso{
\code{ridgePrep}
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

# find optimal penalty parameters
\donttest{optLambdas <- optPenaltyPrep.kCVauto(Y, ids,             }
\donttest{                                     lambdaInit=c(1,1),  }
\donttest{                                     fold=nrow(Y),       }
\donttest{                                     CVcrit="LL")        }

# estimate the precision matrices
\donttest{Ps <- ridgePrep(Y, ids, optLambdas[1], optLambdas[2])    }
}

