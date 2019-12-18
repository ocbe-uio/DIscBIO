#' @title title
#' @export
#' @rdname plotSilhouette
#' @param object object
#' @param K K
setGeneric(
  name = "plotSilhouette",
  def = function(object,K) standardGeneric("plotSilhouette")
)
#' @title title
#' @description description
#' @importFrom cluster silhouette
#' @rdname plotSilhouette
#' @export
setMethod(
  f = "plotSilhouette",
  signature = "PSCANseq",
  definition = function(object, K) {
    if (length(object@kmeans$kpart) == 0) {
      stop("run clustexp before plotsilhouette")
    }
    if (length(unique(object@kmeans$kpart)) < 2) {
      stop("only a single cluster: no silhouette plot")
    }
    col <- c("black", "blue", "green", "red", "yellow", "gray")
		kpart <- object@kmeans$kpart
    distances <- dist.gen(object@distances)
    si <- silhouette(kpart,distances)
    plot(si,col=col[1:K])
  }
)