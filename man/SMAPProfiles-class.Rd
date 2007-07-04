\name{SMAPProfiles-class}
\docType{class}
\alias{SMAPProfiles-class}

\title{Class "SMAPProfiles"}
\description{Holds results from running \code{\link[SMAP:smap]{smap}}.}
\section{Objects from the Class}{
  Objects are not intended to be created directly but as a result from
  running \code{\link[SMAP:smap]{smap}}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} A list of objects
	  of class \code{\link[SMAP:SMAPProfile-class]{SMAPProfile-class}}.}
	\item{\code{name}:}{Object of class \code{"character"} The name of
	  the object.}
  }
}
\section{Extends}{
  Class \code{"list"}, from data part.
  Class \code{"vector"}, by class \code{"list"}.
}
\section{Methods}{
  \describe{
	\item{Q}{\code{signature(object = "SMAPProfiles")}: Returns the
	  optimal state sequence of the list elements.}
	\item{observations}{\code{signature(object = "SMAPProfiles")}:
	  Returns the observations of the list elements.}
	\item{name}{\code{signature(object = "SMAPProfiles")}: Returns the
	  name of the profile.}
    \item{plot}{\code{signature(x = "SMAPProfiles", y = "missing")}: A
	  plot method for the result profiles.}
  }
}
\seealso{
  \code{\link[SMAP:smap]{smap}},
  \code{\link[SMAP:SMAPProfile-class]{SMAPProfile-class}}
}
\author{Robin Andersson, \email{robin.andersson@lcb.uu.se}}
\references{
 Andersson, R., Bruder, C. E. G., Piotrowski, A., Menzel, U., Nord, H.,
  Sandgren, J., Hvidsten, T. R., Diaz de Stahl, T., Dumanski, J. P.,
Komorowski, J., A Segmental Maximum A Posteriori Approach to Array-CGH
Copy Number Profiling, submitted
}
\keyword{classes}
