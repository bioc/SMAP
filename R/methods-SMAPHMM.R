######################################################
##             INITIALIZE METHOD                    ##
######################################################

## Initialize method for class SMAPHMM
setMethod("initialize", "SMAPHMM",
          function(.Object,
                   noStates,
                   Phi,
                   A=NULL,
                   Pi=rep(1/noStates,noStates),
                   initTrans=0.2/(noStates-1)) {

              ## Check valid number of states
              if (noStates < 1)
                  stop("noStates must be > 0")

              ## Check validity of distribution parameter matrix
              Phi.dim <- dim(Phi)
              if (Phi.dim[1] != noStates)
                  stop("Number of rows in Phi must be of length noStates")
              if (Phi.dim[2] != 2)
                  stop(paste("Number of cols in Phi must be of length 2",
                             "(mean and SD)"))
              Phi.list <- lapply(1:Phi.dim[1],function(i) {
                  gaussparam(Phi[i,1],Phi[i,2])
              })
              .Object@Phi <- Phi.list

              .Object@noStates <- noStates
              states <- 1:noStates

              ## Check validity of A matrix
              if (is.null(A)) {
                  if ((initTrans * noStates) > 1)
                      stop(paste("initTrans * noStates must be smaller than",
                                 "(or equal to) 1"))
                  A <- matrix(initTrans, ncol=noStates, nrow=noStates,
                              dimnames=list(states,states))
                  for (i in 1:noStates)
                      A[i,i] <- 1-initTrans*(noStates-1)
                  .Object@A <- A
              } else {
                  dimnames(A) <- list(states,states)
                  .Object@A <- A
                  A.dim <- dim(A)
                  if ((A.dim[1] != A.dim[2]) || (A.dim[1] != noStates))
                      stop("A must be a matrix of size noStates * noStates")
                  if (length(which(sapply(1:noStates,function(i){
                      round(sum(A[i,]),3)}) != 1) > 0))
                      stop("Rows in A matrix must sum to 1")
              }

              ## Check validity of Pi vector
              .Object@Pi <- Pi
              Pi.len <- length(Pi)
              if (Pi.len != noStates)
                  stop("Pi must be a vector of size noStates")
              if (round(sum(Pi),3) != 1)
                  stop("Pi vector must sum to 1")

              .Object@Z <- .Object@A
              .Object@Y <- .Object@Pi
              .Object@eta <- new("eta",value=0,noStates=.Object@noStates)
              ## Return new object
              .Object
          })

SMAPHMM <- function(noStates,
                    Phi,
                    A=NULL,
                    Pi=rep(1/noStates,noStates),
                    initTrans=0.2/(noStates-1)) {

    new("SMAPHMM", noStates=noStates, Phi=Phi, A=A, Pi=Pi,
        initTrans=initTrans)
}

setMethod("show", "SMAPHMM", function(object) {

    cat("An object of class \"SMAPHMM\"\n")

    cat("\nSlot \"A\":\n")
    print(A(object))

    cat("\nSlot \"Pi\":\n")
    print(Pi(object))

    cat("\nSlot \"Phi\":\n")
    print(Phi(object))
})

######################################################
##             ACCESSOR METHODS                     ##
######################################################

setMethod("A", "SMAPHMM", function(object) object@A)
setMethod("Pi", "SMAPHMM", function(object) object@Pi)
setMethod(".Phi", "SMAPHMM", function(object) object@Phi)
setMethod("Phi", "SMAPHMM", function(object){
    phi.m <- t(sapply(1:noStates(object),function(i){
        c(gaussMean(.Phi(object)[[i]]),gaussSd(.Phi(object)[[i]]))
    }))

    dimnames(phi.m) <- list(1:noStates(object),c("mean", "SD"))
    phi.m
})
setMethod("noStates", "SMAPHMM", function(object) object@noStates)
setMethod(".Z", "SMAPHMM", function(object) object@Z)
setMethod(".Y", "SMAPHMM", function(object) object@Y)
setMethod(".eta", "SMAPHMM", function(object) object@eta)
setMethod(".grad", "SMAPHMM", function(object) object@grad)

######################################################
##          REPLACEMENT METHODS                     ##
######################################################

setReplaceMethod("A", "SMAPHMM", function(x, value) {
    ## Check validity of A matrix
    A.dim <- dim(value)
    if ((A.dim[1] != A.dim[2]) || (A.dim[1] != noStates(x)))
        stop("A must be a matrix of size noStates * noStates")
    if (length(which(sapply(1:noStates(x),function(i){
        round(sum(value[i,]),3)}) != 1) > 0))
        stop("Rows in A matrix must sum to 1")
    x@A <- value
    x
})

setReplaceMethod("Pi", "SMAPHMM", function(x, value) {
    ## Check validity of Pi vector
    len <- length(value)
    if (len != noStates(x))
        stop("Pi vector must be of size noStates")
    if (round(sum(value),3) != 1)
        stop("Pi vector must sum to 1")
    x@Pi <- value
    x
})

setReplaceMethod(".Phi", "SMAPHMM", function(x, value) {
    ## Check validity of Phi vector
    len <- length(value)
    if (len != noStates(x))
        stop("Phi must be of length noStates")
    for (i in 1:len)
        if (!is(value[[i]],"gaussparam"))
            stop("Phi must be a list of valid extensions to class gaussparam")
    x@Phi <- value
    x
})

setReplaceMethod(".Z", "SMAPHMM", function(x, value) {
    ## Check validity of Z matrix
    Z.dim <- dim(value)
    if ((Z.dim[1] != Z.dim[2]) || (Z.dim[1] != noStates(x)))
        stop("Z must be a matrix of size noStates * noStates")
    x@Z <- value
    x
})

setReplaceMethod(".Y", "SMAPHMM", function(x, value) {
    ## Check validity of Y vector
    len <- length(value)
    if (len != noStates(x))
        stop("Y vector must be of size noStates")
    x@Y <- value
    x
})

setReplaceMethod(".eta", "SMAPHMM", function(x, value) {
    x@eta <- value
    x
})

setReplaceMethod(".grad", "SMAPHMM", function(x, value) {
    x@grad <- value
    x
})

######################################################
##          VITERBI ALGORITHM                       ##
######################################################

setMethod(".viterbi", signature("SMAPHMM", "SMAPObservations"),
          function(x, Obs, mean.ref, sd.min, mean.sd, W.A, W.Pi,
                   overlap=TRUE, distance=TRUE, L=2000000)
      {

          if (noObservations(Obs) == 0)
              stop("Obs of non-zero length required")

          chroms <- chroms(Obs)
          chrom.sep <- c(chrom.start(Obs), noObservations(Obs)+1)

          no.chroms <- length(chroms)

          P <- 0
          Q <- NULL

          for (i in 1:no.chroms) {

              start <- chrom.sep[i]
              end <- chrom.sep[i+1] - 1
              no.obs <- (end - start) + 1

              prior <- (i == no.chroms)

              clones <- start:end
              olap.start <- startOverlaps(Obs)[start]
              olap.clones <-
                  olap.start:(startOverlaps(Obs)[end] +
                              noOverlaps(Obs)[end] - 1)

              res <- .C("viterbi",
                        as.integer(no.obs),                           # _T
                        as.integer(noStates(x)),                      # _N
                        as.double(A(x)),                              # _A
                        as.double(Pi(x)),                             # _Pi
                        as.double(sapply(1:noStates(x),function(i){
                            gaussMean(.Phi(x)[[i]])})),               # mu
                        as.double(sapply(1:noStates(x),function(i){
                            gaussSd(.Phi(x)[[i]])})),                 # sigma
                        as.double(value(Obs)[clones]),                # obs
                        as.integer(overlap),                          # overlap
                        as.double(overlaps(Obs)[olap.clones]),        # overlaps
                        as.integer(overlapIds(Obs)[olap.clones]-
                                   start),                          # overlap_ids
                        as.integer(noOverlaps(Obs)[clones]),        # no_overlaps
                        as.integer(startOverlaps(Obs)[clones]-
                                   olap.start),                  # start_overlaps
                        as.integer(distance),
                        as.integer(L),
                        as.integer(distance(Obs)[clones]),
                        P=double(1),                                  # P
                        Q=integer(no.obs),                            # Q
                        as.double(mean.ref),                          # mean_ref
                        as.double(sd.min),                            # sd_min
                        as.double(mean.sd),                           # mean_sd
                        as.integer(prior),
                        as.double(W.A),
                        as.double(W.Pi))

              P <- P + res$P
              Q <- c(Q, res$Q)
          }

          list(P=P, Q=Q+1)
      })

######################################################
##            GRADIENT DESCENT ALGORITHM            ##
######################################################

setMethod(".gradient.descent", signature("SMAPHMM", "SMAPObservations"),
          function(x, Obs, Q, P, mean.ref, sd.min, mean.sd,
                   W.A, W.Pi, max.iters=Inf, tau=0.05,
                   eta=0.005, e.change=0.5, e.same=1.2,
                   e.min=0.0001, e.max=0.5, adaptive=TRUE,
                   overlap=TRUE, distance=TRUE, verbose=1,
                   L=2000000)
      {

          if (length(Q) != noObservations(Obs))
              stop("Q and Obs are of different length")

          if (noObservations(Obs) == 0)
              stop("Obs of non-zero length required")

          states <- 1:noStates(x)

          if (!all(unique(Q) %in% states))
              stop("Q contains states not considered by the HMM")

          no.states <- noStates(x)
          no.obs <- noObservations(Obs)

          .Z(x) <- A(x)
          .Y(x) <- Pi(x)

          inf.iters <- 0
          if (max.iters == Inf) {
              max.iters <- 0
              inf.iters <- 1
          }

          chroms <- chroms(Obs)
          chrom.start <- chrom.start(Obs)

          res <- .C("gradient_descent",
                    as.integer(no.obs),                             # _T
                    as.integer(no.states),                          # _N
                    Pi=as.double(Pi(x)),                            # Pi
                    A=as.double(A(x)),                              # _A
                    Y=as.double(.Y(x)),                             # Y
                    Z=as.double(.Z(x)),                             # _Z
                    mu=as.double(sapply(1:no.states,function(i){
                        gaussMean(.Phi(x)[[i]])})),                 # mu
                    sigma=as.double(sapply(1:no.states,function(i){
                        gaussSd(.Phi(x)[[i]])})),                   # sigma
                    as.integer(Q-1),                                # Q
                    P = as.double(P),                               # _P
                    as.integer(sapply(1:no.states,function(i){
                        length(which(Q == i))})),                   # len
                    as.double(value(Obs)),                          # obs
                    as.integer(overlap),                            # overlap
                    as.double(overlaps(Obs)),                       # overlaps
                    as.integer(overlapIds(Obs)-1),                  # overlap_ids
                    as.integer(noOverlaps(Obs)),                    # no_overlaps
                    as.integer(startOverlaps(Obs)-1),            # start_overlaps
                    as.integer(distance),
                    as.integer(L),
                    as.integer(distance(Obs)),
                    as.integer(chrom.start-1),                     # chrom_starts
                    as.integer(length(chroms)),                     # chroms
                    as.double(mean.ref),                            # mean_ref
                    as.double(sd.min),                              # sd_min
                    as.double(mean.sd),                             # mean_sd
                    as.integer(max.iters),                          # max_iters
                    as.integer(inf.iters),                          # inf_iters
                    as.double(tau),                                 # _tau
                    as.double(eta),                                 # eta
                    as.double(e.change),                            # e_change
                    as.double(e.same),                              # e_same
                    as.double(e.min),                               # e_min
                    as.double(e.max),                               # e_max
                    as.integer(adaptive),                           # adaptive
                    as.integer(verbose),                            # verbose
                    as.double(W.A),
                    as.double(W.Pi))

          A(x) <- matrix(res$A, ncol=no.states)
          .Z(x) <- matrix(res$Z, ncol=no.states)
          Pi(x) <- res$Pi
          .Y(x) <- res$Y
          Phi <- lapply(1:no.states, function(i){
              gaussparam(res$mu[i], res$sigma[i])
          })
          .Phi(x) <- Phi

          list(hmm=x, P=res$P)
      })

######################################################
##                 S M A P                          ##
######################################################

setMethod("smap", signature("SMAPHMM", "SMAPObservations"),
          function(x, Obs, sd.min=0.05, mean.sd=0.05,
                   max.iters=Inf, gd.max.iters=Inf, tau=0.05,
                   eta=0.01, e.change=0.5, e.same=1.2,
                   e.min=0.0001, e.max=0.5, adaptive=TRUE,
                   overlap=TRUE, distance=TRUE,
                   chrom.wise=FALSE, verbose=1,
                   L=5000000) {

              ## Check parameters
              if (sd.min <= 0)
                  stop("sd.min must be > 0")

              if (mean.sd <= 0)
                  stop("mean.sd must be > 0")

              if (tau < 0)
                  stop("tau must be >= 0")
              else if (tau == 0) {
                  if ((max.iters == Inf) || (gd.max.iters == Inf))
                      stop("tau must be > 0 if infinite number of iterations")
              }

              if (eta <= 0)
                  stop("eta must be > 0")

              if (adaptive) {
                  if ((e.change <= 0) || (e.change > 1))
                      stop("e.change must be in interval (0,1]")

                  if (e.same < 1)
                      stop("e.same must be >= 1")

                  if (e.max <= 0)
                      stop ("e.max must be > 0")
              }

              if (distance && (L <= 0))
                  stop("L must be > 0")

              Obs.list <- list()
              if (chrom.wise)
                  Obs.list <- .split.on.chrom(Obs)
              else
                  Obs.list <- list(Obs)

              no.elem <- length(Obs.list)
              SMAPProfiles <- new("SMAPProfiles", name=name(Obs))

              mean.ref <- Phi(x)[,1]

              W.A <- matrix(1, ncol=noStates(x), nrow=noStates(x))
              W.Pi <- rep(1, noStates(x))

              for(i in 1:no.elem) {

                  hmm <- x

                  o <- Obs.list[[i]]

                  if (verbose > 1)
                      cat("Calculating overlaps\n")
                  o <- .calc.overlaps(o, overlap)

                  if (verbose > 0)
                      cat(paste("RUNNING SMAP ON \'", name(o), "\'\n", sep=""))

                  if (verbose > 2)
                      cat("**** running viterbi alg. ****\n")

                  res <- .viterbi(hmm, o, mean.ref, sd.min, mean.sd,
                                  W.A=W.A, W.Pi=W.Pi,
                                  overlap=overlap, distance=distance, L=L)

                  Q <- res$Q
                  P <- res$P

                  if (verbose > 0)
                      cat(paste("init P:", round(P,6), "\n"))

                  iteration <- 0

                  opt <- FALSE
                  while (!opt && (iteration < max.iters)) {

                      if (verbose > 2)
                          cat("**** running gradient descent ****\n")

                      grad.res <- .gradient.descent(hmm, o, Q, P, mean.ref,
                                                    sd.min, mean.sd,
                                                    W.A, W.Pi,
                                                    gd.max.iters,
                                                    tau, eta, e.change, e.same,
                                                    e.min, e.max, adaptive,
                                                    overlap, distance,
                                                    verbose, L)

                      new.hmm <- grad.res$hmm

                      if (verbose > 2)
                          cat("**** running viterbi alg. ****\n")

                      res <- .viterbi(new.hmm, o, mean.ref, sd.min, mean.sd,
                                      W.A=W.A, W.Pi=W.Pi,
                                      overlap=overlap, distance=distance, L=L)

                      if (res$P > P + tau) {

                          Q <- res$Q
                          P <- res$P

                          hmm <- new.hmm

                          iteration <- iteration + 1

                          if (P <= grad.res$P + tau)
                              opt <- TRUE

                          if (verbose > 1) {
                              cat(paste("Iteration ", iteration, ", P: ",
                                        round(P,6), "\n", sep=""))
                          }
                      } else {
                          opt <- TRUE
                          if (verbose > 1) {
                              cat(paste("Iteration ", iteration+1, ", P: ",
                                        round(res$P,6), "\n", sep=""))
                          }
                      }
                  }

                  if (verbose > 0)
                      cat(paste("Optimal P:", round(P,6), " found after",
                                iteration,
                                "iterations\n"))

                  SMAPProfiles[[i]] <-
                      SMAPProfile(HMM=hmm, observations=o, P=P,
                                  Q=as.numeric(Q), name(o))
              }

              if (!chrom.wise)
                  SMAPProfiles <- SMAPProfiles[[1]]

              SMAPProfiles
          })
