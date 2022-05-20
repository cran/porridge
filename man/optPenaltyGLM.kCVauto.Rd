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
optPenaltyGLM.kCVauto(Y, X, U=matrix(ncol=0, nrow=length(Y)), lambdaInit, 
                      lambdaGinit=0, Dg=matrix(0, ncol=ncol(X), nrow=ncol(X)),
                      model="linear", target=rep(0, ncol(X)), 
                      folds=makeFoldsGLMcv(min(10, length(X)), Y, model=model),
                      loss="loglik", lambdaMin=10^(-5), 
                      lambdaGmin=10^(-5), minSuccDiff=10^(-5), maxIter=100, 
                      implementation="org")
}
\arguments{
\item{Y}{           A \code{numeric} being the response vector. }
\item{X}{           The design \code{matrix}. The number of rows should match the number of elements of \code{Y}. }
\item{U}{           The design \code{matrix} of the unpenalized covariates. The number of rows should match the number of elements of \code{Y}. }
\item{lambdaInit}{  A \code{numeric}, the initial (starting) values for the regular ridge penalty parameter. }
\item{lambdaGinit}{ A \code{numeric}, the initial (starting) values for the generalized ridge penalty parameter. }
\item{Dg}{          A non-negative definite \code{matrix} of the unscaled generalized ridge penalty. }
\item{model}{       A \code{character}, either \code{"linear"} and \code{"logistic"} (a reference to the models currently implemented), indicating which generalized linear model model instance is to be fitted }
\item{target}{      A \code{numeric} towards which the estimate is shrunken. }
\item{folds}{       A \code{list}. Each list item representing a fold. It is an \code{integer} vector indexing the samples that comprise the fold. This object can be generated with the \code{makeFoldsGLMcv} function. }
\item{loss}{        A \code{character}, either \code{loss="loglik"} or \code{"sos"}, specifying loss criterion to be used in the cross-validation. Used only if \code{model="linear"}. }\item{lambdaMin }{ A positive \code{numeric}, the lower bound of search interval of the regular ridge penalty parameter. }
\item{lambdaGmin }{ A positive \code{numeric} larger than \code{lambdaGmin}, the upper bound of the search interval of the generalized ridge penalty parameter. }
\item{minSuccDiff}{ A \code{numeric}, the minimum distance between the loglikelihoods of two successive iterations to be achieved. Used only if \code{model="logistic"}. }
\item{maxIter}{     A \code{numeric} specifying the maximum number of iterations. Used only if \code{model="logistic"}. }
\item{implementation}{ A \code{character}, either \code{"org"} or \code{"alt"}, specifying the implementation to be used. The implementations (should) only differ in computation efficiency. }
}
\note{
The joint selection of penalty parameters \eqn{\lambda} and \eqn{\lambda_g} through the optimization of the cross-validated loss may lead to a locally optimal choice. This is due to the fact that the penalties are to some extent communicating vessels. Both shrink towards the same target, only in slightly (dependending on the specifics of the generalized penalty matrix \eqn{\Delta}) different ways. As such, the shrinkage achieved by one penalty may be partially compensated for by the other. This may hamper the algorithm in its search for the global optimizers.

Moreover, the penalized IRLS (Iterative Reweighted Least Squares) algorithm for the evaluation of the generalized ridge logistic regression estimator and implemented in the \code{ridgeGLM}-function may fail to converge for small penalty parameter values in combination with a nonzero shrinkage target. This phenomenon propogates to the \code{optPenaltyGLM.kCVauto}-function.
} 
\value{
The function returns a all-positive \code{numeric}, the cross-validated optimal penalty parameters. The average loglikelihood over the left-out samples is used as the cross-validation criterion. If \code{model="linear"}, also the average sum-of-squares over the left-out samples is offered as cross-validation criterion.
}
\references{
van Wieringen, W.N. Binder, H. (2022), "Sequential learning of regression models by penalized estimation", \emph{accepted}.

Lettink, A., Chinapaw, M.J.M., van Wieringen, W.N. et al. (2022), "Two-dimensional fused targeted ridge regression for health indicator prediction from accelerometer data", \emph{submitted}.
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
optLambda <- optPenaltyGLM.kCVauto(Y, X, lambdaInit=1, fold=5, 
                                  target=betas/2, model="logistic", 
                                  minSuccDiff=10^(-3))

# estimate the logistic regression parameter
bHat <- ridgeGLM(Y, X, lambda=optLambda, target=betas/2, model="logistic")
}

