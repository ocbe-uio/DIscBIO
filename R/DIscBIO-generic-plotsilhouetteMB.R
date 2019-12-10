#' @title title
#' @export
#' @rdname plotsilhouetteMB
#' @param K k
setGeneric("plotsilhouetteMB", function(object,K) standardGeneric("plotsilhouetteMB"))
#' @title title
#' @description description
#' @param object object
#' @importFrom cluster silhouette
#' @rdname plotsilhouetteMB
#' @export
setMethod("plotsilhouetteMB",
          signature = "PSCANseq",
          definition = function(object){
            if ( length(object@MBclusters$clusterid) == 0 ) stop("run exprmclust before plotsilhouetteMB")
            if ( length(unique(object@MBclusters$clusterid)) < 2 ) stop("only a single cluster: no silhouette plot")
		col=c("black","blue","green","red","yellow","gray")
            kpart <- object@MBclusters$clusterid
            distances  <- dist.gen(t(object@fdata))
            si <- silhouette(kpart,distances)
            plot(si,col=col[1:K])
          }
          )