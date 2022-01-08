\name{optPenaltyGLMmultiT.kCVauto}
\alias{optPenaltyGLMmultiT.kCVauto}
\title{
Automatic search for optimal penalty parameters of the targeted ridge GLM estimator.
}
\description{
Function finds the optimal penalty parameter of the targeted ridge 
regression estimator of the generalized linear model parameter. The 
optimum is defined as the minimizer of the cross-validated loss 
associated with the estimator. 
}
\usage{
optPenaltyGLMmultiT.kCVauto(Y, X, lambdaInit, model="linear", targetMat,
                      folds=makeFoldsGLMcv(min(10, length(X)), Y, model=model),
                      loss="loglik", lambdaMin=10^(-5),
                      minSuccDiff=10^(-5), maxIter=100)     
}
\arguments{
\item{Y}{           A \code{numeric} being the response vector. }
\item{X}{           The design \code{matrix}. The number of rows should match the number of elements of \code{Y}. }
\item{lambdaInit}{  A \code{numeric} giving the starting values for search of the optimal penalty parameter. }
\item{model}{       A \code{character}, either \code{"linear"} and \code{"logistic"} (a reference to the models currently implemented), indicating which generalized linear model model instance is to be fitted. }
\item{targetMat}{   A \code{matrix} with targets for the regression parameter as columns. }
\item{folds}{       A \code{list}. Each list item representing a fold. It is an \code{integer} vector indexing the samples that comprise the fold. This object can be generated with the \code{makeFoldsGLMcv} function. }
\item{loss}{        A \code{character}, either \code{loss="loglik"} or \code{"sos"}, specifying loss criterion to be used in the cross-validation. Used only if \code{model="linear"}. }
\item{lambdaMin }{ A positive \code{numeric}, the lower bound of search interval of the regular ridge penalty parameter. }
\item{minSuccDiff}{ A \code{numeric}, the minimum distance between the loglikelihoods of two successive iterations to be achieved. Used only if \code{model="logistic"}. }
\item{maxIter}{     A \code{numeric} specifying the maximum number of iterations. Used only if \code{model="logistic"}. }
}
\value{
The function returns an all-positive \code{numeric}, the cross-validated optimal penalty parameters. The average loglikelihood over the left-out samples is used as the cross-validation criterion. If \code{model="linear"}, also the average sum-of-squares over the left-out samples is offered as cross-validation criterion.
}
\references{
van Wieringen, W.N. Binder, H. (2020), "Online learning of regression models from a sequence of datasets by penalized estimation", \emph{submitted}.
}
\author{
W.N. van Wieringen.
}
\examples{
# set the sample size
n <- 50

# set the true parameter
betas <- (c(0:100) - 50) / 20

# generate covariate data
X <- matrix(rnorm(length(betas)*n), nrow=n)

# sample the response
probs <- exp(tcrossprod(betas, X)[1,]) / (1 + exp(tcrossprod(betas, X)[1,]))
Y     <- numeric()
for (i in 1:n){
    Y <- c(Y, sample(c(0,1), 1, prob=c(1-probs[i], probs[i])))
}

# create targets
targets <- cbind(betas/2, rep(0, length(betas)))

# tune the penalty parameter
### optLambdas <- optPenaltyGLMmultiT.kCVauto(Y, X, c(50,0.1), fold=5,            
###                                           targetMat=targets, model="logistic", 
###                                           minSuccDiff=10^(-3))                  

# estimate the logistic regression parameter
### bHat <- ridgeGLMmultiT(Y, X, lambdas=optLambdas, targetMat=targets, model="logistic") 
}


