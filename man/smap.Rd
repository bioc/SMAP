\name{smap}
\docType{methods}
\alias{smap}

\title{smap: A Segmental Maximum A Posteriori Approach to Array-CGH Copy
  Number Profiling}
\description{
This function fits a Hidden Markov Model (HMM) to a set of observed
microarray intensity ratios and outputs the most plausible state sequence
in the HMM through segmental a posteriori maximization.

Briefly, given
an HMM with initial parameter settings \emph{lambda} and a set of
observations \emph{O}, the method alternates maximization of the joint
posterior probability of the state sequence \emph{Q} and \emph{lambda}
given \emph{O}, \code{p(Q,lambda|O)}, over \emph{Q} (using a modified
Viterbi algorithm) and \emph{lambda} (using a gradient descent scheme
with individual learning rate adaptation).
}
\usage{
smap(x, Obs, sd.min=0.05, mean.sd=0.05,
     max.iters=Inf, gd.max.iters=Inf, tau=0.05,
     eta=0.01, e.change=0.5, e.same=1.2,
     e.min=0.0001, e.max=0.5, adaptive=TRUE,
     overlap=TRUE, distance=TRUE, chrom.wise=FALSE,
     verbose=1, L=5000000)
}
\arguments{
  \item{x}{An object of class
	\code{\link[SMAP:SMAPHMM-class]{SMAPHMM-class}}.}
  \item{Obs}{An object of class
	\code{\link[SMAP:SMAPObservations-class]{SMAPObservations-class}}.}
  \item{sd.min}{The minimum allowed standard deviation of state
	associated Gaussian distributions (numeric).}
  \item{mean.sd}{Prior standard deviation of state associated Gaussian
	means (numeric).}
  \item{max.iters}{Maximum number of iterations in the SMAP algorithm
	(numeric).}
  \item{gd.max.iters}{Maximum number of iterations in the gradient
	descent algorithm per SMAP iteration (numeric).}
  \item{tau}{Minimum log probability improvement required in the SMAP
	and gradient descent optimization (numeric).}
  \item{eta}{Initial learning rate in the gradient descent optimization
	(numeric).}
  \item{e.change}{Multiplier for individual learning rate adaptation if
	the sign of partial derivative changes (numeric). Only used if
	\code{adaptive == TRUE}.}
  \item{e.same}{Multiplier for individual learning rate adaptation if
	the sign of partial derivative stays the same (numeric). Only used
	if \code{adaptive == TRUE}.}
  \item{e.min}{Minimum allowed learning rate (numeric).}
  \item{e.max}{Maximum allowed learning rate (numeric).}
  \item{adaptive}{If \code{TRUE}, individual learning rate adaptation
	according to Algorithm 1 in Bagos et al. (2004) is used in the
	gradient descent optimization.
  }
  \item{overlap}{If \code{TRUE}, genomic overlap of clones is considered
	in the optimization.}
  \item{distance}{If \code{TRUE}, genomic distance between clones is considered
	in the optimization, in terms of distance based transition probabilities.}
  \item{chrom.wise}{If \code{TRUE}, the observations are analyzed
	chromosome-wise rather than genome-wise.}
  \item{verbose}{Specifies the amount of output produced; 0 means no
	information and 3 a lot of information (numeric).}
  \item{L}{A positive length parameter that controls the convergence of distance
	based transition probabilities towards 1 / \code{noStates(x)} (numeric).}
}
\details{
  \code{sd.min},
  \code{mean.sd}, and \code{eta} must all be greater than 0.  \code{tau}
  must be greater than 0 if \code{max.iters} or \code{gd.max.iters} are
  infinite, and can be 0 otherwise.  If \code{adaptive} is \code{TRUE},
  then \code{e.change} is required to be in the interval (0,1],
  \code{e.same} must be greater than or equal to 1, and \code{e.max}
  must be greater than 0.
}
\value{
  The method returns an object of class
  \code{\link[SMAP:SMAPProfile-class]{SMAPProfile-class}} or
  \code{\link[SMAP:SMAPProfiles-class]{SMAPProfiles-class}} if
  \code{chrom.wise} is set to \code{FALSE} or \code{TRUE}, respectively.
}
\seealso{
  \code{\link[SMAP:SMAPHMM]{SMAPHMM}},
  \code{\link[SMAP:SMAPObservations]{SMAPObservations}}
}
\author{Robin Andersson \email{robin.andersson@lcb.uu.se}}
\references{
  Andersson, R., Bruder, C. E. G., Piotrowski, A., Menzel, U., Nord, H.,
  Sandgren, J., Hvidsten, T. R., Diaz de Stahl, T., Dumanski, J. P.,
  Komorowski, J., A Segmental Maximum A Posteriori Approach to Array-CGH
  Copy Number Profiling, submitted
  
  Bagos P. G., Liakopoulos T. D., Hamodrakas, S. J. (2004) Faster
  Gradient Descent Training of Hidden Markov Models, Using Individual
  Learning Rate Adaptation. In Paliouras, G., Sakakibara, Y., editors,
  \emph{ICGI}, volume 3264 of \emph{Lecture Notes in Computer Science},
  pages 40--52.
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
plot(observations, ylim=c(0,2))
## Initiate HMM
init.means <- c(0.4, 0.7, 1, 1.3, 1.6, 3)
init.sds <- rep(0.1, 6)
phi <- cbind(init.means, init.sds)
hmm <- SMAPHMM(6, phi, initTrans=0.02)
hmm
## RUN SMAP:
profile <- smap(hmm, observations, verbose=2)
## genome profile
plot(profile, ylim=c(0,2))
## chromosome 9 profile
ids <- which(chromosome(observations) == "9")
plot(profile[ids], ylim=c(0,2), main="chromosome 9")
## output results for chromosome 9
#cbind(reporterId(observations[ids]), Q(profile[ids]))
}
\keyword{methods}
