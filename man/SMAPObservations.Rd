\name{SMAPObservations}
\alias{SMAPObservations}

\title{Constructor for "SMAPObservations" objects}
\description{
  A constructor for
  \code{\link[SMAP:SMAPObservations-class]{SMAPObservations-class}} objects.
}
\usage{
SMAPObservations(value, chromosome, startPosition, endPosition,
                 name=character(0),
                 reporterId=as.character(1:length(value)))
}
\arguments{
  \item{value}{A vector of microarray intensity ratios (numeric).}
  \item{chromosome}{A vector of chromosome annotations (character).}
  \item{startPosition}{A vector of start positions (numeric).}
  \item{endPosition}{A vector of end positions (numeric).}
  \item{name}{The name of the observation set (character).}
  \item{reporterId}{A vector of observation identifiers, e.g., probe ids
	(character).}
}
\details{
  The vectors \code{value}, \code{chromosome}, \code{startPosition},
  \code{endPosition}, and \code{reporterId} must be of equal length.
}
\value{
  An object of class
  \code{\link[SMAP:SMAPObservations-class]{SMAPObservations-class}}.
}
\seealso{
  \code{\link[SMAP:smap]{smap}},
  \code{\link[SMAP:SMAPObservations-class]{SMAPObservations-class}},
  \code{\link[SMAP:SMAPHMM-class]{SMAPHMM-class}}
}
\author{Robin Andersson, \email{robin.andersson@lcb.uu.se}}
\references{
  Andersson, R., Bruder, C. E. G., Piotrowski, A., Menzel, U., Nord, H.,
  Sandgren, J., Hvidsten, T. R., Diaz de Stahl, T., Dumanski, J. P.,
  Komorowski, J., A Segmental Maximum A Posteriori Approach to Array-CGH
  Copy Number Profiling, submitted
}
\examples{
## Load Glioblastoma multiforme data
data(GBM)
observations <- SMAPObservations(value=as.numeric(GBM[,2]),
                                 chromosome=as.character(GBM[,3]),
                                 startPosition=as.numeric(GBM[,4]),
                                 endPosition=as.numeric(GBM[,5]),
                                 name="G24460",
                                 reporterId=as.character(GBM[,1]))
## plot observations
plot(observations, ylim=c(0,2))
## plot subset of observations (chromosome 9)
ids <- which(chromosome(observations) == "9")
plot(observations[ids])
}
\keyword{classes}
