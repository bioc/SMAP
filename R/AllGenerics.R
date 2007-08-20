## GENERIC ACCESSOR METHODS:

## SMAPObservations
if (!isGeneric("name"))
  setGeneric("name", function(object) standardGeneric("name"))
if (!isGeneric("value"))
  setGeneric("value", function(object) standardGeneric("value"))
if (!isGeneric("noObservations"))
  setGeneric("noObservations", function(object)
             standardGeneric("noObservations"))
if (!isGeneric("chromosome"))
  setGeneric("chromosome", function(object) standardGeneric("chromosome"))
if (!isGeneric("chroms"))
  setGeneric("chroms", function(object) standardGeneric("chroms"))
if (!isGeneric("chrom.start"))
  setGeneric("chrom.start", function(object) standardGeneric("chrom.start"))
if (!isGeneric("startPosition"))
  setGeneric("startPosition", function(object)
             standardGeneric("startPosition"))
if (!isGeneric("endPosition"))
  setGeneric("endPosition", function(object) standardGeneric("endPosition"))
if (!isGeneric("reporterId"))
  setGeneric("reporterId", function(object) standardGeneric("reporterId"))
if (!isGeneric("distance"))
  setGeneric("distance", function(object) standardGeneric("distance"))
if (!isGeneric("overlapIds"))
  setGeneric("overlapIds", function(object) standardGeneric("overlapIds"))
if (!isGeneric("overlaps"))
  setGeneric("overlaps", function(object) standardGeneric("overlaps"))
if (!isGeneric("startOverlaps"))
  setGeneric("startOverlaps", function(object) standardGeneric("startOverlaps"))
if (!isGeneric("noOverlaps"))
  setGeneric("noOverlaps", function(object) standardGeneric("noOverlaps"))

## SMAPHMM
if (!isGeneric("A"))
  setGeneric("A", function(object) standardGeneric("A"))
if (!isGeneric("Pi"))
  setGeneric("Pi", function(object) standardGeneric("Pi"))
if (!isGeneric("Phi"))
  setGeneric("Phi", function(object) standardGeneric("Phi"))
if (!isGeneric(".Phi"))
  setGeneric(".Phi", function(object) standardGeneric(".Phi"))
if (!isGeneric("noStates"))
  setGeneric("noStates", function(object) standardGeneric("noStates"))
if (!isGeneric(".Z"))
  setGeneric(".Z", function(object) standardGeneric(".Z"))
if (!isGeneric(".Y"))
  setGeneric(".Y", function(object) standardGeneric(".Y"))
if (!isGeneric(".eta"))
  setGeneric(".eta", function(object) standardGeneric(".eta"))
if (!isGeneric(".grad"))
  setGeneric(".grad", function(object) standardGeneric(".grad"))

## SMAPProfile
if (!isGeneric("HMM"))
    setGeneric("HMM", function(object) standardGeneric("HMM"))
if (!isGeneric("observations"))
    setGeneric("observations", function(object)
               standardGeneric("observations"))
if (!isGeneric("P"))
    setGeneric("P", function(object) standardGeneric("P"))
if (!isGeneric("Q"))
    setGeneric("Q", function(object) standardGeneric("Q"))

## gaussparam
if (!isGeneric("gaussMean"))
  setGeneric("gaussMean", function(object) standardGeneric("gaussMean"))
if (!isGeneric("gaussSd"))
  setGeneric("gaussSd", function(object) standardGeneric("gaussSd"))

## GENERIC REPLACEMENT METHODS:

## SMAPObservations
if (!isGeneric("value<-"))
  setGeneric("value<-", function(x, value) standardGeneric("value<-"))
if (!isGeneric("overlaps<-"))
  setGeneric("overlaps<-", function(x, value) standardGeneric("overlaps<-"))
if (!isGeneric("overlapIds<-"))
  setGeneric("overlapIds<-", function(x, value)
             standardGeneric("overlapIds<-"))
if (!isGeneric("startOverlaps<-"))
  setGeneric("startOverlaps<-", function(x, value)
             standardGeneric("startOverlaps<-"))
if (!isGeneric("noOverlaps<-"))
  setGeneric("noOverlaps<-", function(x, value) standardGeneric("noOverlaps<-"))

## SMAPHMM
if (!isGeneric("A<-"))
  setGeneric("A<-", function(x, value) standardGeneric("A<-"))
if (!isGeneric("Pi<-"))
  setGeneric("Pi<-", function(x, value) standardGeneric("Pi<-"))
if (!isGeneric("Phi<-"))
  setGeneric("Phi<-", function(x, value) standardGeneric("Phi<-"))
if (!isGeneric(".Phi<-"))
  setGeneric(".Phi<-", function(x, value) standardGeneric(".Phi<-"))
if (!isGeneric(".Z<-"))
  setGeneric(".Z<-", function(x, value) standardGeneric(".Z<-"))
if (!isGeneric(".Y<-"))
  setGeneric(".Y<-", function(x, value) standardGeneric(".Y<-"))
if (!isGeneric(".eta<-"))
  setGeneric(".eta<-", function(x, value) standardGeneric(".eta<-"))
if (!isGeneric(".grad<-"))
  setGeneric(".grad<-", function(x, value) standardGeneric(".grad<-"))

## gaussparam
if (!isGeneric("gaussMean<-"))
  setGeneric("gaussMean<-", function(x, value) standardGeneric("gaussMean<-"))
if (!isGeneric("gaussSd<-"))
  setGeneric("gaussSd<-", function(x, value) standardGeneric("gaussSd<-"))

## GENERIC METHODS:
if (!isGeneric(".split.on.chrom"))
    setGeneric(".split.on.chrom", function(Obs)
               standardGeneric(".split.on.chrom"))

if (!isGeneric("getChromObs"))
    setGeneric("getChromObs", function(Obs, c) standardGeneric("getChromObs"))

if (!isGeneric(".viterbi"))
    setGeneric(".viterbi", function(x, Obs, mean.ref, sd.min, mean.sd,
                                    W.A, W.Pi,
                                    overlap=TRUE, distance=TRUE, L=2000000)
               standardGeneric(".viterbi"))

if (!isGeneric(".calc.overlaps"))
  setGeneric(".calc.overlaps", function(Obs, overlap=TRUE)
             standardGeneric(".calc.overlaps"))

if (!isGeneric(".gradient.descent"))
    setGeneric(".gradient.descent", function(x, Obs, Q, P, mean.ref,
                                             sd.min, mean.sd,
                                             W.A, W.Pi,
                                             max.iters=Inf, tau=0.05,
                                             eta=0.005, e.change=0.5,
                                             e.same=1.2,
                                             e.min=0.0001, e.max=0.5,
                                             adaptive=TRUE,
                                             overlap=TRUE, distance=TRUE,
                                             verbose=1, L=2000000)
               standardGeneric(".gradient.descent"))

if (!isGeneric("smap"))
  setGeneric("smap", function(x, Obs, sd.min=0.05, mean.sd=0.05,
                              max.iters=Inf, gd.max.iters=Inf, tau=0.05,
                              eta=0.01, e.change=0.5, e.same=1.2,
                              e.min=0.0001, e.max=0.5, adaptive=TRUE,
                              overlap=TRUE, distance=TRUE,
                              chrom.wise=FALSE, verbose=1,
                              L=5000000)
             standardGeneric("smap"))

if (!isGeneric(".draw.dist"))
  setGeneric(".draw.dist", function(x,col) standardGeneric(".draw.dist"))

if (!isGeneric("profilePlot"))
    setGeneric("profilePlot", function(profile, ...)
               standardGeneric("profilePlot"))
