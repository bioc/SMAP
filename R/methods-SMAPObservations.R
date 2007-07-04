######################################################
##             INITIALIZE METHOD                    ##
######################################################

## Initialize method for class SMAPObservations
setMethod("initialize", "SMAPObservations",
          function(.Object,
                   value,
                   chromosome,
                   startPosition,
                   endPosition,
                   name=character(0),
                   reporterId=as.character(1:length(value))) {

              .Object@name <- name

              ## Check lengths of required vectors
              len <- c(length(value),
                       length(chromosome),
                       length(startPosition),
                       length(endPosition),
                       length(reporterId))
              len.un <- unique(len)
              if (length(len.un) > 1)
                  stop(paste("value, chromosome, startPosition, endPosition,",
                             "and reporterId vectors of unequal lengths"))
              if (len.un == 0)
                  stop(paste("value, chromosome, startPosition, endPosition,",
                             "and reporterId vectors must be of length > 0"))

              ## Check and remove NAs in data
              nas <- (is.na(value) | is.na(chromosome) |
                      is.na(startPosition) | is.na(endPosition))

              if (any(nas)) {
                  if (any(!nas)) {
                      value <- value[!nas]
                      chromosome <- chromosome[!nas]
                      startPosition <- startPosition[!nas]
                      endPosition <- endPosition[!nas]
                      reporterId <- reporterId[!nas]
                  } else {
                      stop(paste("value, chromosome, startPosition,",
                                 "and endPosition vectors must contain",
                                 "at least one non NA element"))
                  }
              }

              if (any(endPosition - startPosition < 0))
                  stop(paste("end positions must be greater or equal to start",
                             "positions"))

              ## Order by chromosome, start position and end position
              chr.un <- unique(chromosome)
              suppressWarnings(non.numeric <- is.na(as.numeric(chr.un)))
              chr.un <- c(sort(as.numeric(chr.un[!non.numeric])),
                          sort(chr.un[non.numeric]))
              sort.order <- NULL
              distance <- NULL
              for (c in chr.un) {
                  c.ids <- which(chromosome == c)
                  e.ids <- c.ids[order(startPosition[c.ids], endPosition[c.ids])]
                  sort.order <- c(sort.order, e.ids)

                  ## Calculate distances between clones
                  no.clones <- length(e.ids)
                  if (no.clones > 1) {
                      dist <- c(0, startPosition[e.ids][-1] - endPosition[e.ids][-no.clones])
                      dist[dist < 0] <- 0
                      distance <- c(distance, dist)
                  } else {
                      distance <- c(distance, 0)
                  }
              }

              ## Store value in object
              .Object@value <- value[sort.order]
              .Object@chromosome <- chromosome[sort.order]
              .Object@chroms <- unique(.Object@chromosome)
              .Object@chrom.start <- match(.Object@chroms, .Object@chromosome)
              .Object@startPosition <- startPosition[sort.order]
              .Object@endPosition <- endPosition[sort.order]
              .Object@reporterId <- reporterId[sort.order]
              .Object@distance <- distance

              no.obs <- length(value)
              .Object@noObservations <- no.obs

              ## Return new object
              .Object
          })

SMAPObservations <- function(value, chromosome, startPosition, endPosition,
                             name=character(0),
                             reporterId=as.character(1:length(value))) {

    new("SMAPObservations", value=value, chromosome=chromosome,
        startPosition=startPosition, endPosition=endPosition,
        name=name, reporterId=reporterId)
}

setMethod("show", "SMAPObservations", function(object) {

    cat("An object of class \"SMAPObservations\"\n")

    cat(paste("\nName:", name(object), "\n"))

    cat(paste("\nNumber of Observations:", noObservations(object), "\n"))

    cat(paste("Use methods value(object), chromosome(object),",
              "startPosition(object), endPosition(object), and",
              "reporterId(object) to access object slots.\n"))
})

######################################################
##             ACCESSOR METHODS                     ##
######################################################

setMethod("name", "SMAPObservations", function(object) object@name)
setMethod("value", "SMAPObservations", function(object) object@value)
setMethod("noObservations", "SMAPObservations",
          function(object) object@noObservations)
setMethod("chromosome", "SMAPObservations", function(object) object@chromosome)
setMethod("chroms", "SMAPObservations", function(object) object@chroms)
setMethod("chrom.start", "SMAPObservations", function(object) object@chrom.start)
setMethod("startPosition", "SMAPObservations",
          function(object) object@startPosition)
setMethod("endPosition", "SMAPObservations",
          function(object) object@endPosition)
setMethod("reporterId", "SMAPObservations", function(object) object@reporterId)
setMethod("distance", "SMAPObservations", function(object) object@distance)
setMethod("overlapIds", "SMAPObservations", function(object) object@overlapIds)
setMethod("overlaps", "SMAPObservations", function(object) object@overlaps)
setMethod("startOverlaps", "SMAPObservations", function(object)
          object@startOverlaps)
setMethod("noOverlaps", "SMAPObservations", function(object) object@noOverlaps)

setMethod("[", "SMAPObservations",
          function(x, i, j, ..., drop) {
              x@value <- x@value[i]
              x@chromosome <- x@chromosome[i]
              x@chroms <- unique(x@chromosome)
              x@chrom.start <- match(x@chroms, x@chromosome)
              x@startPosition <- x@startPosition[i]
              x@endPosition <- x@endPosition[i]
              x@reporterId <- x@reporterId[i]
              x@noObservations <- length(x@value)

              distance <- NULL
              for (c in x@chroms) {
                  c.ids <- which(x@chromosome == c)
                  ## Calculate distances between clones
                  no.clones <- length(c.ids)
                  if (no.clones > 1) {
                      dist <- c(0, x@startPosition[c.ids][-1] - x@endPosition[c.ids][-no.clones])
                      dist[dist < 0] <- 0
                      distance <- c(distance, dist)
                  } else {
                      distance <- c(distance, 0)
                  }
              }
              x@distance <- distance

              x
          })

######################################################
##          REPLACEMENT METHODS                     ##
######################################################

setReplaceMethod("value", "SMAPObservations", function(x, value) {
    x@noObservations <- length(value)
    x@value <- value
    x
})

setReplaceMethod("overlaps", "SMAPObservations", function(x, value) {
    x@overlaps <- value
    x
})

setReplaceMethod("overlapIds", "SMAPObservations", function(x, value) {
    x@overlapIds <- value
    x
})

setReplaceMethod("startOverlaps", "SMAPObservations", function(x, value) {
    x@startOverlaps <- value
    x
})

setReplaceMethod("noOverlaps", "SMAPObservations", function(x, value) {
    x@noOverlaps <- value
    x
})

######################################################
##          CHROMOSOME EXTRACTION METHODS           ##
######################################################

setMethod(".split.on.chrom", "SMAPObservations",
          function(Obs) {

              chrom.un <- unique(chromosome(Obs))

              lapply(chrom.un, function(c) {
                  getChromObs(Obs, c)
              })
          })

setMethod("getChromObs", "SMAPObservations",
          function(Obs, c) {

              ids <- chromosome(Obs) == c

              if (!any(ids))
                  stop(paste("Unknown chromosome:", c))

              SMAPObservations(value=value(Obs)[ids],
                               chromosome=chromosome(Obs)[ids],
                               startPosition=startPosition(Obs)[ids],
                               endPosition=endPosition(Obs)[ids],
                               name=paste(name(Obs),"_",c,sep=""),
                               reporterId=reporterId(Obs)[ids])
          })

######################################################
##          OVERLAP CALCULAION METHODS              ##
######################################################

setMethod(".calc.overlaps", "SMAPObservations",
          function(Obs, overlap=TRUE) {

              no.obs <- noObservations(Obs)

              ## Calculate overlaps
              overlaps <- NULL
              overlapIds <- NULL
              startOverlaps <- NULL
              noOverlaps <- NULL

              if (overlap) {

                  chroms <- chroms(Obs)
                  chromosome <- chromosome(Obs)
                  int.chrom <- NULL
                  sapply(1:length(chroms), function(c) {
                      int.chrom[which(chromosome == chroms[c])] <<- c})

                  res <- .Call("calc_overlaps",
                               no.obs,
                               startPosition(Obs),
                               endPosition(Obs),
                               int.chrom)

                  overlapIds <- res[[1]]
                  overlaps <- res[[2]]
                  startOverlaps <- res[[3]]
                  noOverlaps <- res[[4]]

              } else {
                  overlapIds <- 1:no.obs
                  overlaps <- rep(0, no.obs)
                  startOverlaps <- overlapIds
                  noOverlaps <- rep(1, no.obs)
              }
              overlapIds(Obs) <- overlapIds
              overlaps(Obs) <- overlaps
              startOverlaps(Obs) <- startOverlaps
              noOverlaps(Obs) <- noOverlaps

              Obs
          })

######################################################
##            VISUALIZATION METHODS                 ##
######################################################

setMethod("plot", signature(x="SMAPObservations", y="missing"),
          function(x, y, ...) {

              middle.pos <- startPosition(x) +
                  (endPosition(x) - startPosition(x)) / 2

              chromosome <- chromosome(x)
              startPosition <- startPosition(x)
              endPosition <- endPosition(x)

              chroms <- chroms(x)

              ## Multiple chromosomes?
              no.chroms <- length(chroms)
              mult.chroms <- no.chroms > 1

              ## Manipulate positions of spots according to chromosome
              ## membership
              chrstart <- 0
              cmean <- NULL
              .axes <- TRUE
              .xlab <- "position"
              .main <- name(x)
              op <- NULL

              if (mult.chroms) {
                  sapply(2:no.chroms, function(c){
                      chrstart[c] <<- chrstart[c-1] +
                          max(middle.pos[which(chromosome == chroms[c-1])],
                              na.rm=TRUE)
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
                  args[["y"]] <- value(x)
              }

              args[["true.lengths"]] <- NULL

              do.call("plot", args)

              if (true.lengths) {
                  value <- value(x)
                  for (i in 1:noObservations(x))
                      points(c(startPosition[i], endPosition[i]),
                             rep(value[i], 2), type="l")
              }

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
