######################################################
##             INITIALIZE METHOD                    ##
######################################################

## Initialize method for class SMAPProfile
setMethod("initialize", "SMAPProfile",
          function(.Object,
                   HMM,
                   observations,
                   P,
                   Q,
                   name=character(0)) {

              ## Store values in object
              .Object@HMM <- HMM
              .Object@observations <- observations
              .Object@P <- P
              .Object@Q <- Q
              .Object@name <- name

              ## Check lengths of required vectors
              len <- c(noObservations(observations),
                       length(Q))
              len.un <- unique(len)
              if (length(len.un) > 1)
                  stop("observations and state sequence of unequal lengths")
              if (len.un == 0)
                  stop("observations and state sequence must be of length > 0")

              ## Check state consistency
              state.un <- unique(Q)
              nas <- is.na(state.un)
              if (any(!state.un[!nas] %in% 1:noStates(HMM)))
                  stop("Unrecognized states in state sequence")

              ## Return new object
              .Object
          })

SMAPProfile <- function(HMM, observations, P, Q, name=character(0)) {

    new("SMAPProfile", HMM=HMM, observations=observations,
        P=P, Q=Q, name=name)
}

setMethod("show", "SMAPProfile", function(object) {

    cat("An object of class \"SMAPProfile\"\n")

    cat(paste("\nP:", P(object), "\n"))

    cat(paste("Use methods HMM(object), observations(object),",
              "P(object), and Q(object)",
              "to access object slots.\n"))
})

setMethod("show", "SMAPProfiles", function(object) {

    cat("An object of class \"SMAPProfiles\"\n")

    cat(paste("\nName:", name(object), "\n"))

    show(object@.Data)
})

######################################################
##             ACCESSOR METHODS                     ##
######################################################

setMethod("HMM", "SMAPProfile", function(object) object@HMM)
setMethod("observations", "SMAPProfile", function(object) object@observations)
setMethod("P", "SMAPProfile", function(object) object@P)
setMethod("Q", "SMAPProfile", function(object) object@Q)
setMethod("name", "SMAPProfile", function(object) object@name)
setMethod("name", "SMAPProfiles", function(object) object@name)

setMethod("observations", "SMAPProfiles",
          function(object) {

              obs <- observations(object[[1]])
              obs@name <- name(object)

              no.profiles <- length(object)

              if (no.profiles > 1) {
                  for (i in 2:no.profiles) {

                      p.obs <- observations(object[[i]])

                      obs@value <- c(obs@value, p.obs@value)
                      obs@chromosome <- c(obs@chromosome, p.obs@chromosome)
                      obs@startPosition <- c(obs@startPosition,
                                             p.obs@startPosition)
                      obs@endPosition <- c(obs@endPosition, p.obs@endPosition)
                      obs@reporterId <- c(obs@reporterId, p.obs@reporterId)
                  }
              }
              obs@noObservations <- length(obs@value)
              obs@chroms <- unique(obs@chromosome)
              obs@chrom.start <- match(obs@chroms, obs@chromosome)

              obs
          })

setMethod("Q", "SMAPProfiles",
          function(object) {

              Q <- Q(object[[1]])

              no.profiles <- length(object)

              if (no.profiles > 1)
                  for (i in 2:no.profiles)
                      Q <- c(Q, Q(object[[i]]))

              Q
          })

setMethod("[", "SMAPProfile",
          function(x, i, j, ..., drop) {
              x@observations <- x@observations[i]
              x@Q <- x@Q[i]
              x
          })

######################################################
##          REPLACEMENT METHODS                     ##
######################################################

setReplaceMethod("[[", "SMAPProfiles",
                 function(x, i, j, value) {
                     if (!is(value, "SMAPProfile"))
                         stop("value must be of class SMAPProfile")
                     x@.Data[[i]] <- value
                     x
                 })

######################################################
##            VISUALIZATION METHODS                 ##
######################################################

setMethod("plot", signature(x="SMAPProfile", y="missing"),
          function(x, y, ...) profilePlot(x, ...))

setMethod("plot", signature(x="SMAPProfiles", y="missing"),
          function(x, y, ...) profilePlot(x, ...))

setMethod("profilePlot", signature("SMAPProfile"), function(profile, ...) {

    hmm <- HMM(profile)
    Obs <- observations(profile)
    Q <- Q(profile)

    middle.pos <- startPosition(Obs) +
        (endPosition(Obs) - startPosition(Obs)) / 2

    chromosome <- chromosome(Obs)
    startPosition <- startPosition(Obs)
    endPosition <- endPosition(Obs)

    chroms <- chroms(Obs)

    args <- list(...)

    .ylim <- args[["ylim"]]
    if (is.null(.ylim))
        args[["ylim"]] <- c(min(sapply(1:noStates(hmm),function(i){
            gaussMean(.Phi(hmm)[[i]]) -
                gaussSd(.Phi(hmm)[[i]])}),min(value(Obs),na.rm=TRUE)),
                            max(sapply(1:noStates(hmm),function(i){
                                gaussMean(.Phi(hmm)[[i]]) +
                                    gaussSd(.Phi(hmm)[[i]])}),
                                max(value(Obs),na.rm=TRUE)))

    ## Multiple chromosomes?
    no.chroms <- length(chroms)
    mult.chroms <- no.chroms > 1

    ## Manipulate positions of spots according to chromosome membership
    chrstart <- 0
    cmean <- NULL
    .axes <- TRUE
    .xlab <- "position"
    .main <- name(profile)
    op <- NULL

    if (mult.chroms) {
        sapply(2:no.chroms, function(c){
            chrstart[c] <<- chrstart[c-1] +
                max(middle.pos[which(chromosome == chroms[c-1])], na.rm=TRUE)
        })
        sapply(1:no.chroms, function(c){
            chrom <- chromosome == chroms[c]
            middle.pos[chrom] <<-
                middle.pos[chrom] + chrstart[c]
            cmean[c] <<- mean(middle.pos[chrom], na.rm=TRUE)
        })
        .axes <- FALSE
        .xlab <- "chromosome"

        op <- par(mar=c(5, 4, 6, 2))
    }

    title <- args[["main"]]
    if (is.null(title))
        args[["main"]] <- .main

    .ylab <- args[["ylab"]]
    if (is.null(.ylab))
        args[["ylab"]] <- "value"

    args[["xlab"]] <- .xlab
    args[["axes"]] <- .axes

    true.lengths <- args[["true.lengths"]]
    if (!mult.chroms && !is.null(true.lengths) && true.lengths) {
        args[["x"]] <- 1
        args[["y"]] <- 1
        args[["type"]] <- "n"
        .xlim <- args[["xlim"]]
        if (is.null(.xlim))
            args[["xlim"]] <- c(min(startPosition), max(startPosition))
    } else {
        true.lengths <- FALSE
        args[["x"]] <- middle.pos
        args[["y"]] <- value(Obs)
        args[["col"]] <- Q
    }

    args[["true.lengths"]] <- NULL
    draw.dists <- args[["draw.dists"]]
    args[["draw.dists"]] <- NULL

    do.call("plot", args)

    if (true.lengths) {
        value <- value(Obs)
        for (i in 1:noObservations(Obs))
            points(c(startPosition[i], endPosition[i]),
                   rep(value[i], 2), type="l", col=Q[i])
    }

    if (is.null(draw.dists))
        draw.dists <- TRUE

    if (draw.dists)
        for (i in 1:noStates(hmm))
            .draw.dist(.Phi(hmm)[[i]],col=i)

    if (mult.chroms) {
        abline(v=chrstart[2:no.chroms], col="grey", lty=2)
        box()
        odd <- seq(1, no.chroms, by=2)
        even <- seq(2, no.chroms, by=2)
        axis(1, at=cmean[odd], labels=chroms[odd])
        axis(3, at=cmean[even], labels=chroms[even])
        axis(2)
        par(op)
    }
})

setMethod("profilePlot", signature("SMAPProfiles"), function(profile, ...) {

    Obs <- observations(profile)
    Q <- Q(profile)

    middle.pos <- startPosition(Obs) +
        (endPosition(Obs) - startPosition(Obs)) / 2

    chromosome <- chromosome(Obs)

    chroms <- chroms(Obs)

    ## Multiple chromosomes?
    no.chroms <- length(chroms)
    mult.chroms <- no.chroms > 1

    ## Manipulate positions of spots according to chromosome membership
    chrstart <- 0
    cmean <- NULL
    .axes <- TRUE
    .xlab <- "position"
    .main <- name(profile)
    op <- NULL

    if (mult.chroms) {
        sapply(2:no.chroms, function(c){
            chrstart[c] <<- chrstart[c-1] +
                max(middle.pos[which(chromosome == chroms[c-1])], na.rm=TRUE)
        })
        sapply(1:no.chroms, function(c){
            chrom <- chromosome == chroms[c]
            middle.pos[chrom] <<-
                middle.pos[chrom] + chrstart[c]
            cmean[c] <<- mean(middle.pos[chrom], na.rm=TRUE)
        })
        .axes <- FALSE
        .xlab <- "chromosome"

        op <- par(mar=c(5, 4, 6, 2))
    }

    args <- list(...)

    title <- args[["main"]]
    if (is.null(title))
        args[["main"]] <- .main

    .ylab <- args[["ylab"]]
    if (is.null(.ylab))
        args[["ylab"]] <- "value"

    args[["xlab"]] <- .xlab
    args[["axes"]] <- .axes
    args[["x"]] <- middle.pos
    args[["y"]] <- value(Obs)
    args[["col"]] <- Q

    do.call("plot", args)

    if (mult.chroms) {
        abline(v=chrstart[2:no.chroms], col="grey", lty=2)
        box()
        odd <- seq(1, no.chroms, by=2)
        even <- seq(2, no.chroms, by=2)
        axis(1, at=cmean[odd], labels=chroms[odd])
        axis(3, at=cmean[even], labels=chroms[even])
        axis(2)
        par(op)
    }
})
