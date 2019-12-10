#' PSCANseq
#' @slot expdata    
#' @slot ndata      
#' @slot fdata      
#' @slot distances  
#' @slot tsne       
#' @slot background 
#' @slot out        
#' @slot cpart      
#' @slot fcol       
#' @slot filterpar  
#' @slot clusterpar 
#' @slot outlierpar 
#' @slot kmeans     
#' @slot MBclusters 
#' @slot kordering  
#' @slot MBordering 
#' @slot MBtsne     
#' @importFrom methods new
#' @name PSCANseq
#' @rdname PSCANseq
#' @aliases PSCANseq-class, PSCANseq-class
#' @exportClass PSCANseq
PSCANseq <- setClass(
    Class = "PSCANseq",
    slots = c(
        expdata    = "data.frame",
        ndata      = "data.frame",
        fdata      = "data.frame", 
        distances  = "matrix",
        tsne       = "data.frame",
        background = "list",
        out        = "list", 
        cpart      = "vector",
        fcol       = "vector",
        filterpar  = "list",
        clusterpar = "list", 
        outlierpar = "list",
        kmeans     = "list",
        MBclusters = "vector",
        kordering  = "vector",
        MBordering = "vector",
        MBtsne     = "data.frame"
    )
)

setValidity("PSCANseq",
            function(object) {
              msg <- NULL
              if ( ! is.data.frame(object@expdata) ){
                msg <- c(msg, "input data must be data.frame")
              }else if ( nrow(object@expdata) < 2 ){
                msg <- c(msg, "input data must have more than one row")
              }else if ( ncol(object@expdata) < 2 ){
                msg <- c(msg, "input data must have more than one column")
              }else if (sum( apply( is.na(object@expdata),1,sum ) ) > 0 ){
                msg <- c(msg, "NAs are not allowed in input data")
              }else if (sum( apply( object@expdata,1,min ) ) < 0 ){
                msg <- c(msg, "negative values are not allowed in input data")
              }
              if (is.null(msg)) TRUE
              else msg
            }
            )

#' @title title
#' @description description
#' @param .Object .Object
#' @param expdata expdata
#' @importFrom methods validObject
#' @rdname initialize
setMethod("initialize",
          signature = "PSCANseq",
          definition = function(.Object, expdata ){
            .Object@expdata <- expdata
            .Object@ndata <- expdata
            .Object@fdata <- expdata
            validObject(.Object)
            return(.Object)
          }
          )