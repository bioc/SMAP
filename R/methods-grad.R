######################################################
##             INITIALIZE METHOD                    ##
######################################################

## Initialize method for class grad
setMethod("initialize", "grad",
          function(.Object, noStates) {
              .Object@A <- matrix(0,ncol=noStates,nrow=noStates)
              .Object@Pi <- rep(0,noStates)
              .Object@Phi <- lapply(1:noStates,function(i){
                  new("gaussparam", mean=0, sd=0)})
              ## Return new object
              .Object
          })

grad <- function(noStates=1) {
    new("grad", noStates=noStates)
}

######################################################
##             ACCESSOR METHODS                     ##
######################################################

setMethod("A", "grad", function(object) object@A)
setMethod("Pi", "grad", function(object) object@Pi)
setMethod("Phi", "grad", function(object) object@Phi)

######################################################
##          REPLACEMENT METHODS                     ##
######################################################

setReplaceMethod("A", "grad", function(x, value) {
    x@A <- value
    x
})

setReplaceMethod("Pi", "grad", function(x, value) {
    x@Pi <- value
    x
})

setReplaceMethod("Phi", "grad", function(x, value) {
    x@Phi <- value
    x
})
