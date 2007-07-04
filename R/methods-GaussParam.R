######################################################
##             INITIALIZE METHOD                    ##
######################################################

## Initialize method for class gaussparam
setMethod("initialize", "gaussparam",
          function(.Object,
                   mean,
                   sd) {
              .Object@mean <- mean
              .Object@sd <- sd

              .Object
          })

gaussparam <- function(mean=0, sd=1) {
    new("gaussparam", mean=mean, sd=sd)
}

######################################################
##             ACCESSOR METHODS                     ##
######################################################

setMethod("gaussMean", "gaussparam", function(object) object@mean)
setMethod("gaussSd", "gaussparam", function(object) object@sd)

######################################################
##          REPLACEMENT METHODS                     ##
######################################################

setReplaceMethod("gaussMean", "gaussparam", function(x, value){
    x@mean <- value
    x
})

setReplaceMethod("gaussSd", "gaussparam", function(x, value){
    x@sd <- value
    x
})

######################################################
##            VISUALIZATION METHODS                 ##
######################################################

setMethod(".draw.dist", "gaussparam",
          function(x,col) {
              abline(h=gaussMean(x),lty=1,col=col)
              abline(h=gaussMean(x) + gaussSd(x),lty=2,col=col)
              abline(h=gaussMean(x) - gaussSd(x),lty=2,col=col)
          })

