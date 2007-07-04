\name{SMAPHMM}
\alias{SMAPHMM}

\title{Constructor for "SMAPHMM" objects}
\description{
  A constructor for \code{\link[SMAP:SMAPHMM-class]{SMAPHMM-class}} objects.
}
\usage{
SMAPHMM(noStates, Phi, A=NULL,
        Pi=rep(1/noStates,noStates),
        initTrans=0.2/(noStates - 1))
}
\arguments{
  \item{noStates}{The number of hidden states in the HMM (numeric).}
  \item{Phi}{A Gaussian distribution parameter matrix (numeric).}
  \item{A}{A \code{noStates * noStates} matrix of transition
	probabilities betweeen the hidden states (numeric).}
  \item{Pi}{A vector of initial probabilities of starting in a certain
	state (numeric).}
  \item{initTrans}{Specifies the transition probability between
	non-equal states (numeric).}
}
\details{
  \code{Phi} is a \code{noStates * 2} matrix that specifies the
  paramaters of Gaussian distributions associated with each hidden
  state.  The first column specifies standard deviations, the second
  specifies means.

  If \code{A == NULL}, \code{initTrans} specifies the transition
  probability between states \code{i} and \code{j} in \code{1:noStates},
  such that \code{i != j}. Only used if \code{A ==
	NULL}. \code{initTrans * noStates} must be smaller than (or equal
  to) 1.
}
\value{
  An object of class \code{\link[SMAP:SMAPHMM-class]{SMAPHMM-class}}.
}
\seealso{
  \code{\link[SMAP:smap]{smap}},
  \code{\link[SMAP:SMAPHMM-class]{SMAPHMM-class}},
  \code{\link[SMAP:SMAPObservations-class]{SMAPObservations-class}}
}
\author{Robin Andersson, \email{robin.andersson@lcb.uu.se}}
\references{
 Andersson, R., Bruder, C. E. G., Piotrowski, A., Menzel, U., Nord, H.,
  Sandgren, J., Hvidsten, T. R., Diaz de Stahl, T., Dumanski, J. P.,
Komorowski, J., A Segmental Maximum A Posteriori Approach to Array-CGH
Copy Number Profiling, submitted
}
\keyword{classes}
