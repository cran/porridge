\name{optPenaltyGLM.kCVauto}
\alias{optPenaltyGLM.kCVauto}
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
optPenaltyGLM.kCVauto(Y, X, lambdaInit, fold=nrow(X), stratified=TRUE,
                      target=rep(0, ncol(X)), model, loss="loglik",
                      minSuccDiff=10^(-5), maxIter=100)
}
\arguments{
\item{Y}{           A \code{numeric} being the response vector. }
\item{X}{           The design \code{matrix}. The number of rows should match the number of elements of \code{Y}. }
\item{lambdaInit}{ A \code{numeric} giving the initial (starting) values for the two penalty parameters. }
\item{fold}{ A \code{numeric} or \code{integer} specifying the number of folds to apply in the cross-validation. }
\item{stratified}{ A \code{logical}, should the data by stratified such that range of the response is comparable among folds? }
\item{target}{      A \code{numeric} towards which the estimate is shrunken. }
\item{model}{       A \code{character}, either \code{"linear"} and \code{"logistic"} (a reference to the models currently implemented), indicating which generalized linear model model instance is to be fitted. }
\item{loss}{        A \code{numeric}, loss criterion to be used in the cross-validation. Used only if \code{model="linear"}. }
\item{minSuccDiff}{ A \code{numeric}, the minimum distance between the loglikelihoods of two successive iterations to be achieved. Used only if \code{model="logistic"}. }
\item{maxIter}{     A \code{numeric} specifying the maximum number of iterations. Used only if \code{model="logistic"}. }
}
\value{
The function returns a all-positive \code{numeric}, the cross-validated optimal penalty parameters. The average loglikelihood over the left-out samples is used as the cross-validation criterion. If \code{model="linear"}, also the average sum-of-squares over the left-out samples is offered as cross-validation criterion.
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

# tune the penalty parameter
optLambda <- optPenaltyGLM.kCVauto(Y, X, 1, fold=5, 
                                  target=betas/2, model="logistic", 
                                  minSuccDiff=10^(-3))

# estimate the logistic regression parameter
bHat <- ridgeGLM(Y, X, optLambda, target=betas/2, model="logistic")
}

