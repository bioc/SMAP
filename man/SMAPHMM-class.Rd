\name{SMAPHMM-class}
\docType{class}
\alias{SMAPHMM-class}
\alias{A}
\alias{Phi}
\alias{Pi}
\alias{noStates}

\title{Class "SMAPHMM": A class to manage HMMs for the SMAP package}
\description{Holds parameters for a Hidden Markov Model (HMM) used in
  the \pkg{SMAP} package.}
\section{Objects from the Class}{
  Objects are not intended to be created directly but via the
  constructor function \code{\link[SMAP:SMAPHMM]{SMAPHMM}}.
}
\section{Slots}{
  \describe{
    \item{\code{A}:}{Object of class \code{"matrix"} The transition
	  probability matrix between states.}
    \item{\code{Pi}:}{Object of class \code{"numeric"} The initial
	  probabilities of starting in a certain state.}
    \item{\code{Phi}:}{Object of class \code{"matrix"} A matrix that
	  specifies the paramaters of Gaussian distributions associated with
	  each hidden state. The first column specifies standard deviations,
	  the second specifies means.}
    \item{\code{noStates}:}{Object of class \code{"numeric"} The number
	  of hidden states in the HMM.}
  }
}
\section{Methods}{
  \describe{
    \item{A}{\code{signature(object = "SMAPHMM")}: Returns the
	  transition matrix.}
    \item{Pi}{\code{signature(object = "SMAPHMM")}: Returns the initial
	  probabilities.}
	\item{Phi}{\code{signature(object = "SMAPHMM")}: Returns the
	  distribution parameter matrix.}
	\item{noStates}{\code{signature(object = "SMAPHMM")}: Returns the
	  number of hidden states in the HMM.}
  }
}
\seealso{
  \code{\link[SMAP:smap]{smap}}, \code{\link[SMAP:SMAPHMM]{SMAPHMM}}
}
\author{Robin Andersson, \email{robin.andersson@lcb.uu.se}}
\references{
Andersson, R., Bruder, C. E. G., Piotrowski, A., Menzel, U., Nord, H.,
  Sandgren, J., Hvidsten, T. R., Diaz de Stahl, T., Dumanski, J. P.,
Komorowski, J., A Segmental Maximum A Posteriori Approach to Array-CGH
Copy Number Profiling, submitted
}
\keyword{classes}
