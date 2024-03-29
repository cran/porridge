\docType{package}
\name{porridge-package}
\alias{porridge-package}
\alias{porridge}
\title{Ridge-Type Penalized Estimation of a Potpourri of Models.}
\description{
The following functions facilitate the ridge-type penalized estimation of various models. Currently, it includes: 
\itemize{
\item Generalized ridge estimation of the precision matrix of a Gaussian graphical model (van Wieringen, 2019) through the function \code{\link{ridgePgen}}. This function is complemented by the functions \code{\link{ridgePgen.kCV}}, \code{\link{ridgePgen.kCV.banded}}, \code{\link{ridgePgen.kCV.groups}}, 
\code{\link{optPenaltyPgen.kCVauto.banded}} and \code{\link{optPenaltyPgen.kCVauto.groups}} for penalty parameters selection through K-fold cross-validation assuming a particularly structured precision matrix. 

\item Multi-targeted ridge estimation of the precision matrix of a Gaussian graphical model (van Wieringen et al., 2020) through the functions \code{\link{ridgePmultiT}}. This function is complemented by the functions \code{\link{optPenaltyPmultiT.kCVauto}} for penalty parameters selection through K-fold cross-validation. 

\item Gaussian graphical model estimation from data with replicates in ridge penalized fashion (van Wieringen, Chen, 2021) (\code{\link{ridgePrep}} and \code{\link{ridgePrepEdiag}}). The two functions 
\code{\link{optPenaltyPrep.kCVauto}} and \code{\link{optPenaltyPrepEdiag.kCVauto}} implement the corresponding K-fold cross-validation procedures for an optimal choice of the penalty parameter.

\item Ridge penalized estimation of a mixture of Gaussian graphical models: \code{\link{ridgeGGMmixture}} and its penalty selection via K-fold cross-validation 
\code{\link{optPenaltyGGMmixture.kCVauto}}. 

\item Targeted and multi-targeted ridge estimation of the regression parameter of the generalized linear model (van Wieringen, Binder, 2022; van Wieringen, 2021; Lettink et al., 2022) through the functions \code{\link{ridgeGLM}} and \code{\link{ridgeGLMmultiT}}. This function is complemented by the functions \code{\link{optPenaltyGLM.kCVauto}} and \code{\link{optPenaltyGLMmultiT.kCVauto}} for penalty parameters selection through K-fold cross-validation, and the \code{ridgeGLMdof}-function for the evaluation of the fitted model's degrees of freedom.  
}
Future versions aim to include more ridge-type functionality.

In part the \code{porridge}-package extends/builds upon the \code{\link[rags2ridges:rags2ridges]{rags2ridges}}-packages, in which some or all functionality of the \code{porridge}-package may be absorped at some point in the future.
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.ml>
}
\references{
Aflakparast, M., de Gunst, M.C.M., van Wieringen, W.N. (2018), "Reconstruction of molecular network evolution from cross-sectional omics data", \emph{Biometrical Journal}, 60(3), 547-563.

Lettink, A., Chinapaw, M.J.M., van Wieringen, W.N. (2022), "Two-dimensional fused targeted ridge regression for health indicator prediction from accelerometer data", \emph{submitted}.

Peeters, C.F.W., Bilgrau, A.E., and van Wieringen, W.N. (2021), "\code{rags2ridges}: Ridge Estimation of Precision Matrices from  High-Dimensional Data", R package version 2.2.5. \url{https://CRAN.R-project.org/package=rags2ridges}.

van Wieringen, W.N. (2019), "The generalized ridge estimator of the inverse covariance matrix", \emph{Journal of Computational and Graphical Statistics}, 28(4), 932-942.

van Wieringen, W.N. (2021), "Lecture notes on ridge regression", \emph{Arxiv preprint}, arXiv:1509.09169.

van Wieringen W.N., Chen, Y. (2021), "Penalized estimation of the Gaussian graphical model from data with replicates", \emph{Statistics in Medicine}, 40(19), 4279-4293.

van Wieringen, W.N., Stam, K.A., Peeters, C.F.W., van de Wiel, M.A. (2020), "Updating of the Gaussian graphical model through targeted penalized estimation", \emph{Journal of Multivariate Analysis}, 178, Article 104621.

van Wieringen, W.N. Binder, H. (2022), "Sequential learning of regression models by penalized estimation", \emph{Journal of Computational and Graphical Statistics}, accepted.
}
\keyword{ package }
\seealso{
The \code{porridge}-package.
}


