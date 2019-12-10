#' @title title
#' @export
#' @rdname plottsneMB
#' @param K K
setGeneric("plottsneMB", function(object,K) standardGeneric("plottsneMB"))
#' @title title
#' @description description
#' @param object object
#' @importFrom graphics text
#' @export
#' @rdname plottsneMB
setMethod("plottsneMB",
          signature = "PSCANseq",
          definition = function(object){
            if ( length(object@MBtsne) == 0 ) stop("run comptsneMB before plottsneMB")
		col=c("black","blue","green","red","yellow","gray")
            part <- object@MBclusters$clusterid
            plot(object@MBtsne,xlab="Dim 1",ylab="Dim 2",pch=20,cex=1.5,col="lightgrey",las=1)
            for ( i in 1:K ){
              if ( sum(part == i) > 0 ) text(object@MBtsne[part == i,1],object@MBtsne[part == i,2],i,col=col[i],cex=.75,font=4)
            }
          }
          )    