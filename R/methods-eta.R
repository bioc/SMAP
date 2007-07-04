######################################################
##             INITIALIZE METHOD                    ##
######################################################

## Initialize method for class eta
setMethod("initialize", "eta",
          function(.Object, value, noStates) {
              .Object@value <- value
              .Object@A <- matrix(value,ncol=noStates,nrow=noStates)
              .Object@Pi <- rep(value,noStates)
              .Object@Phi <- lapply(1:noStates,function(i){
                  new("gaussparam", mean=value, sd=value)})
              ## Return new object
              .Object
          })

eta <- function(value=0, noStates=1) {
    new("eta", value=value, noStates=noStates)
}

######################################################
##             ACCESSOR METHODS                     ##
######################################################

setMethod("value", "eta", function(object) object@value)
setMethod("A", "eta", function(object) object@A)
setMethod("Pi", "eta", function(object) object@Pi)
setMethod("Phi", "eta", function(object) object@Phi)

######################################################
##          REPLACEMENT METHODS                     ##
######################################################

setReplaceMethod("value", "eta", function(x, value) {
    x@value <- value
    noStates <- length(x@Pi)
    x@A <- matrix(value,ncol=noStates,nrow=noStates)
    x@Pi <- array(value,noStates)
    x@Phi <- lapply(1:noStates,function(i){
        new("gaussparam", mean=value, sd=value)})
    x
})

setReplaceMethod("A", "eta", function(x, value) {
    x@A <- value
    x
})

setReplaceMethod("Pi", "eta", function(x, value) {
    x@Pi <- value
    x
})

setReplaceMethod("Phi", "eta", function(x, value) {
    x@Phi <- value
    x
})

