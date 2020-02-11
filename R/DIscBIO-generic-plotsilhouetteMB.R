#' @title Silhouette Plot for Model-based clustering
#' @description The silhouette provides a representation of how well each point is represented by its cluster in comparison to the 
#' closest neighboring cluster. It computes for each point the difference between the average similarity 
#' to all points in the same cluster and to all points in the closest neighboring cluster. This difference it 
#' normalize such that it can take values between -1 and 1 with higher values reflecting better 
#' representation of a point by its cluster.
#' @param object \code{DISCBIO} class object.
#' @param K A numeric value of the number of clusters
#' @importFrom stats as.dist cor
#' @importFrom cluster silhouette
#' @return A silhouette plot
setGeneric("plotsilhouetteMB", function(object,K) standardGeneric("plotsilhouetteMB"))

#' @export
#' @rdname plotsilhouetteMB
setMethod("plotsilhouetteMB",
          signature = "DISCBIO",
          definition = function(object, K){
            if ( length(object@MBclusters$clusterid) == 0 ) stop("run exprmclust before plotsilhouetteMB")
            if ( length(unique(object@MBclusters$clusterid)) < 2 ) stop("only a single cluster: no silhouette plot")
      col=c("black","blue","green","red","yellow","gray")
            kpart <- object@MBclusters$clusterid
      dist.gen <- function(x,method="euclidean", ...) if ( method %in% c("spearman","pearson","kendall") ) as.dist( 1 - cor(t(x),method=method,...) ) else dist(x,method=method,...)
            distances  <- dist.gen(t(object@fdata))
            si <- silhouette(kpart,distances)
            plot(si,col=col[1:K])
          }
          )