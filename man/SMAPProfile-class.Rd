\name{SMAPProfile-class}
\docType{class}
\alias{SMAPProfile-class}
\alias{HMM}
\alias{observations}
\alias{P}
\alias{Q}

\title{Class "SMAPProfile"}
\description{Holds results from running \code{\link[SMAP:smap]{smap}}.}
\section{Objects from the Class}{
  Objects are not intended to be created directly but as a result from
  running \code{\link[SMAP:smap]{smap}.}
}
\section{Slots}{
  \describe{
    \item{\code{HMM}:}{Object of class \code{"SMAPHMM"}}
    \item{\code{observations}:}{Object of class \code{"SMAPObservations"}}
    \item{\code{P}:}{Object of class \code{"numeric"} The log joint
	  posterior probability of the state sequence \code{Q} and
	  parameters of \code{HMM} given the \code{observations}.}
    \item{\code{Q}:}{Object of class \code{"numeric"} The optimal state
	  sequance (path) in the HMM.}
	\item{\code{name}:}{Object of class \code{"character"} The name of
	  the object.}
  }
}
\section{Methods}{
  \describe{
    \item{P}{\code{signature(object = "SMAPProfile")}: Returns the log
	  joint posterior probability. }
    \item{Q}{\code{signature(object = "SMAPProfile")}: Returns the
	  optimal state sequence.}
    \item{HMM}{\code{signature(object = "SMAPProfile")}: Returns the
	  optimized HMM.}
	\item{name}{\code{signature(object = "SMAPProfile")}: Returns the
	  name of the profile.}
    \item{observations}{\code{signature(object = "SMAPProfile")}:
	  Returns the observations.}
    \item{plot}{\code{signature(x = "SMAPProfile", y = "missing")}: A
	  plot method for the result profile.}
	\item{[}{\code{signature(x = "SMAPProfile")}: Creates a new object
	  of class SMAPProfile with extracted elements as specified by the
	  indices provided.}
  }
}
\seealso{
  \code{\link[SMAP:smap]{smap}},
  \code{\link[SMAP:SMAPProfiles-class]{SMAPProfiles-class}}
}
\author{Robin Andersson, \email{robin.andersson@lcb.uu.se}}
\references{
 Andersson, R., Bruder, C. E. G., Piotrowski, A., Menzel, U., Nord, H.,
  Sandgren, J., Hvidsten, T. R., Diaz de Stahl, T., Dumanski, J. P.,
Komorowski, J., A Segmental Maximum A Posteriori Approach to Array-CGH
Copy Number Profiling, submitted
}
\keyword{classes}
