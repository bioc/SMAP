\name{GBM}
\alias{GBM}
\docType{data}
\title{Glioblastoma multiforme array CGH data}
\description{
  Array CGH data measurements of glioblastoma multiforme sample G24460.
}
\usage{data(GBM)}
\source{
  Genome wide array CGH data from Diaz de Stahl, T., et \emph{al}. (2005).
}
\references{
  Diaz de Stahl, T., et \emph{al}. (2005) Chromosome 22 tiling-path
  array-CGH analysis identifies germ-line- and tumor-specific
  aberrations in patients with glioblastoma multiforme. \emph{Genes
	Chromosomes Cancer} \bold{44}(2), 161--169
}
\seealso{
  \code{\link[SMAP:smap]{smap}}
}
\examples{
data(GBM)
observations <- SMAPObservations(value=as.numeric(GBM[,2]),
                                 chromosome=as.character(GBM[,3]),
                                 startPosition=as.numeric(GBM[,4]),
                                 endPosition=as.numeric(GBM[,5]),
                                 name="G24460",
                                 reporterId=as.character(GBM[,1]))
plot(observations)
}
\keyword{datasets}
