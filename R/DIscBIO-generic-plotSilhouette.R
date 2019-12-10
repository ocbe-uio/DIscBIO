#' @title title
#' @export
#' @rdname plotSilhouette
#' @param K K
setGeneric("plotSilhouette", function(object,K) standardGeneric("plotSilhouette"))
#' @title title
#' @description description
#' @param object object
#' @importFrom cluster silhouette
#' @rdname plotSilhouette
#' @export
setMethod("plotSilhouette",
          signature = "PSCANseq",
          definition = function(object){
            if ( length(object@kmeans$kpart) == 0 ) stop("run clustexp before plotsilhouette")
            if ( length(unique(object@kmeans$kpart)) < 2 ) stop("only a single cluster: no silhouette plot")
            col=c("black","blue","green","red","yellow","gray")
		kpart <- object@kmeans$kpart
            distances  <- dist.gen(object@distances)
            si <- silhouette(kpart,distances)
            plot(si,col=col[1:K])
          }
          )