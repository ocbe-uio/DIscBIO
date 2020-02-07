#' @title tSNE map for K-means clustering
#' @description Visualizing the K-means clusters using tSNE maps
#' @param object \code{DISCBIO} class object.
#' @importFrom graphics text
#' @examples 
#' sc <- DISCBIO(valuesG1msReduced) # changes signature of data
#' sc <- Clustexp(sc, cln=3) # data must be clustered before plotting
#' sc <- comptSNE(sc, rseed=15555, quiet=TRUE)
#' plottSNE(sc)
setGeneric("plottSNE", function(object) standardGeneric("plottSNE"))

#' @rdname plotSilhouette
#' @export
setMethod("plottSNE",
          signature = "DISCBIO",
          definition = function(object){
            if ( length(object@tsne) == 0 ) stop("run comptsne before plottsne")
			col=c("black","blue","green","red","yellow","gray")
            part <- object@kmeans$kpart
            plot(object@tsne,las=1,xlab="Dim 1",ylab="Dim 2",pch=20,cex=1.5,col="lightgrey")
            for ( i in 1:max(part) ){
              if ( sum(part == i) > 0 ) text(object@tsne[part == i,1],object@tsne[part == i,2],i,col=col[i],cex=.75,font=4)
            }
          }
          )
