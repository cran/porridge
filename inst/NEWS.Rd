\name{NEWS}
\title{News for Package \pkg{porridge}}
\newcommand{\cpkg}{\href{http://cran.r-project.org/package=#1}{\pkg{#1}}}

\section{Changes in version 0.3.0}{
 \itemize{
    \item Functionality of the \code{ridgeGLM} and \code{optPenaltyGLM.kCVauto} expanded. Among other, unpenalized covariates may now be supplied and a generalized penalty matrix provided.
    \item Bug fix: correct stopping criterion for the IRLS algorithm of the penalized ridge logistic regression estimator.
    \item Citations updated.
  }
}

\section{Changes in version 0.2.1}{
 \itemize{
    \item Faster implementation of the functions for target and multi-targeted ridge estimation of the GLM offered.
  }
}

\section{Changes in version 0.2.0}{
 \itemize{
    \item Reference added to the Description-file.
    \item Functions for target and multi-targeted ridge estimation of the GLM added.
  }
}

\section{Changes in version 0.1.1}{
 \itemize{
    \item Corrected DOIs in the Description field.
  }
}

\section{Changes in version 0.1.0}{
 \itemize{
    \item Restructuring of R and C++ files into modules.
    \item Novel module added.
    \item Some layout polish in the manual.
  }
}

\section{Changes in version 0.0.4}{
 \itemize{
    \item Modification of C++ code related to the \code{ridgePrepEdiag}- and \code{optPenaltyPrepEdiag}-functions. In case of an unbalanced design the initialization of EM algorithm no longer yields negative definite precisiom matrices.
  }
}

\section{Changes in version 0.0.3}{
 \itemize{
    \item Modification of input checks to meet CRAN's best practice.  
    \item Removed an accidently inserted section in the manual of the \code{optPenaltyGGMmixture.kCVauto}-function.  
  }
}

\section{Changes in version 0.0.2}{
 \itemize{
    \item The unnecessary dependency on the \code{rags2ridges}-packages has been removed. 
  }
}

\section{Changes in version 0.0.1}{
 \itemize{
    \item Release of the Urfassung. 
  }
}

