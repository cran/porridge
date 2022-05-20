\name{makeFoldsGLMcv}
\alias{makeFoldsGLMcv}
\title{
Generate folds for cross-validation of generalized linear models.
}
\description{
Function that evaluates the targeted ridge estimator of the regression parameter of generalized linear models. 
}
\usage{
makeFoldsGLMcv(fold, Y, stratified=TRUE, model="linear")
}
\arguments{
\item{fold}{ An \code{integer}, the number of folds to be generated. }
\item{Y}{ A \code{numeric} being the response vector. }
\item{stratified}{ A \code{logical}. If \code{stratified=TRUE}, the folds are generated such the distribution of the response \code{Y} is (roughly) the same across folds. }
\item{model}{ A \code{character}, either \code{"linear"} and \code{"logistic"} (a reference to the models currently implemented), indicative of the type of response for stratification. }
}
\value{
A \code{list} of length \code{fold}. Each list item is a fold.
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

# generate folds 
folds <- makeFoldsGLMcv(10, Y, model="logistic")
}
