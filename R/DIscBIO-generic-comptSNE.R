#' @title Computing tSNE for K-means clustering
#' @description This function is used to compute the t-Distributed Stochastic Neighbor Embedding (t-SNE).
#' @param object \code{PSCANseq} class object.
#' @param Integer number. Random seed to to yield exactly reproducible maps across different runs. Default is 15555. 
#' @param quiet if `TRUE`, suppresses intermediate output
#' @importFrom tsne tsne
#' @rdname comptSNE
#' @export
setGeneric(
  name = "comptSNE",
  def = function(object, rseed = 15555, quiet = FALSE) {
    standardGeneric("comptSNE")
  }
)

setMethod(
  "comptSNE",
  signature = "PSCANseq",
  definition = function(object, rseed, quiet = FALSE) {
    if (length(object@kmeans$kpart) == 0) stop("run clustexp before comptsne")
    set.seed(rseed)
	dist.gen <- function(x,method="euclidean", ...) if ( method %in% c("spearman","pearson","kendall") ) as.dist( 1 - cor(t(x),method=method,...) ) else dist(x,method=method,...)
    di <- dist.gen(as.matrix(object@distances))
    if (quiet) {
      ts <- suppressMessages(tsne(di, k = 2))
    } else {
      ts <- tsne(di, k = 2)
    }
      object@tsne <- as.data.frame(ts)
      return(object)
  }
)