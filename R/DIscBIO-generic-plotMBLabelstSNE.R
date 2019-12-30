#' title tSNE map for Model-based clustering with labels
#' @description Visualizing the Model-based clusters using tSNE maps 
#' @param object \code{PSCANseq} class object.
#' @importFrom graphics text
#' @rdname plotMBLabelstSNE
#' @export
setGeneric("plotMBLabelstSNE", function(object) standardGeneric("plotMBLabelstSNE"))

setMethod("plotMBLabelstSNE",
          signature = "PSCANseq",
          definition = function(object){
            if ( length(object@MBtsne) == 0 ) stop("run comptsneMB before plotMBLabelstSNE")
            Clusters<-object@MBclusters$clusterid
            ClustersFactor<- as.factor(Clusters)
            ClustersFactor<- gsub("1", "black", ClustersFactor)
            ClustersFactor<- gsub("2", "blue", ClustersFactor)
            ClustersFactor<- gsub("3", "green", ClustersFactor)
            ClustersFactor<- gsub("4", "red", ClustersFactor)
            ClustersFactor<- gsub("5", "yellow", ClustersFactor)
            ClustersFactor<- gsub("6", "gray", ClustersFactor)
            COL<-ClustersFactor
            labels=names(object@ndata)
            plot(object@MBtsne,xlab="Dim 1",ylab="Dim 2",pch=20,cex=.5,col="lightgrey")
            text(object@MBtsne[,1],object@MBtsne[,2],labels,cex=.7,col=COL)
          }
          ) 