#' PSCANseq
#' @title The PSCANseq Class
#' @description The PSCANseq class is the central object storing all information generated throughout the pipeline. 
#' @slot expdata 
#' @slot expdataAll   
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
#' @slot noiseF
#' @slot FinalGeneList  
#' @importFrom methods new
#' @name PSCANseq
#' @rdname PSCANseq
#' @aliases PSCANseq-class, PSCANseq-class
#' @exportClass PSCANseq
#' @export
PSCANseq <- setClass(
    Class = "PSCANseq",
    slots = c(
        expdata    = "data.frame",
		expdataAll    = "data.frame",
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
        MBtsne     = "data.frame",
		noiseF   = "vector",
		FinalGeneList   = "vector"
    )
)

#' validity function for PSCANseq
#'
#' @param object An PSCANseq object.
#' @name PSCANseq
#' @export
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

setMethod("initialize",
          signature = "PSCANseq",
          definition = function(.Object, expdataAll ){
		    .Object@expdataAll <- expdataAll
			shortNames <- substr(rownames(expdataAll), 1, 4)
			geneTypes <- factor(c(ENSG = "ENSG", ERCC = "ERCC" )[shortNames])
			expdata <- expdataAll[which(geneTypes == "ENSG"), ]
			#countsERCC <- expdataAll[which(geneTypes == "ERCC" ), ]
            .Object@expdata <- expdata
            .Object@ndata <- expdata
            .Object@fdata <- expdata
            validObject(.Object)
            return(.Object)
          }
          )