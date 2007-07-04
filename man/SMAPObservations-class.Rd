\name{SMAPObservations-class}
\docType{class}
\alias{SMAPObservations-class}
\alias{chromosome}
\alias{endPosition}
\alias{name}
\alias{noObservations}
\alias{reporterId}
\alias{startPosition}
\alias{value}
\alias{plot,SMAPObservations,missing-method}

\title{Class "SMAPObservations": A class to manage microarray
  observations for the SMAP package}
\description{Holds observed microarray intensity ratios and clone
  annotations for the \pkg{SMAP} package.}
\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("SMAPObservations", value, chromosome, startPosition,
	endPosition, name, reporterId)}.  You can also use the convenience
  function \code{\link[SMAP:SMAPObservations]{SMAPObservations}}.
}
\section{Slots}{
  \describe{
    \item{\code{value}:}{Object of class \code{"numeric"} Microarray
	  intensity ratios.}
    \item{\code{chromosome}:}{Object of class \code{"character"}
	  Associated chromosomes for the observations.}
    \item{\code{startPosition}:}{Object of class \code{"numeric"}
	  Associated start positions for the observations.}
    \item{\code{endPosition}:}{Object of class \code{"numeric"}
	  Associated end positions for the observations.}
    \item{\code{reporterId}:}{Object of class \code{"character"}
	  Identifiers of the observations, e.g., probe ids.}
	\item{\code{name}:}{Object of class \code{"character"} An identifier
	  of the observation set.}
    \item{\code{noObservations}:}{Object of class \code{"numeric"} The
	  number of observations in the set.}
  }
}
\section{Methods}{
  \describe{
	\item{value}{\code{signature(object = "SMAPObservations")}: Returns
	  the values of the observations.}
	\item{chromosome}{\code{signature(object = "SMAPObservations")}:
	  Returns the chromosome annotations of the observations.}
	\item{startPosition}{\code{signature(object = "SMAPObservations")}:
	  Returns the start positions of the observations.}
    \item{endPosition}{\code{signature(object = "SMAPObservations")}:
	  Returns the end positions of the observations.}
	\item{reporterId}{\code{signature(object = "SMAPObservations")}:
	  Returns the identifiers of the observations.}
    \item{name}{\code{signature(object = "SMAPObservations")}: Returns
	  the name of the observation set.}
    \item{noObservations}{\code{signature(object = "SMAPObservations")}:
	  Returns the number of observations in the set.}
    \item{initialize}{\code{signature(.Object = "SMAPObservations")}:
	  Creates an instance.}
    \item{plot}{\code{signature(x = "SMAPObservations", y = "missing")}:
	  A plot method for the observations.}
	\item{[}{\code{signature(x = "SMAPPObservations")}: Creates a new
	  object of class SMAPObservations with extracted elements as
	  specified by the indices provided.}
  }
}
\seealso{
  \code{\link[SMAP:smap]{smap}},
  \code{\link[SMAP:SMAPObservations]{SMAPObservations}}
}
\author{Robin Andersson, \email{robin.andersson@lcb.uu.se}}
\references{
 Andersson, R., Bruder, C. E. G., Piotrowski, A., Menzel, U., Nord, H.,
  Sandgren, J., Hvidsten, T. R., Diaz de Stahl, T., Dumanski, J. P.,
Komorowski, J., A Segmental Maximum A Posteriori Approach to Array-CGH
Copy Number Profiling, submitted
}
\keyword{classes}
