\name{NEWS}
\title{News for Package \pkg{porridge}}
\newcommand{\cpkg}{\href{http://cran.r-project.org/package=#1}{\pkg{#1}}}

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

