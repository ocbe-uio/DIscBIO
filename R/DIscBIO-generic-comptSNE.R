#' @title Computing tSNE for K-means clustering
#' @description This function is used to compute the t-Distributed Stochastic Neighbor Embedding (t-SNE).
#' @param object \code{DISCBIO} class object.
#' @param rseed Integer number. Random seed to to yield exactly reproducible maps across different runs. Default is 15555. 
#' @param quiet if `TRUE`, suppresses intermediate output
#' @importFrom tsne tsne
#' @importFrom stats as.dist cor
#' @examples
#' sc <- DISCBIO(valuesG1ms) # changes signature of data
#' sc <- Clustexp(sc, cln=3) # data must be clustered before plottin
#' sc <- comptSNE(sc, rseed=15555, quiet=TRUE)
#' head(sc@tsne)
setGeneric(
  name = "comptSNE",
  def = function(object, rseed = 15555, quiet = FALSE) {
    standardGeneric("comptSNE")
  }
)

#' @rdname comptSNE
#' @export
setMethod(
  "comptSNE",
  signature = "DISCBIO",
  definition = function(object, rseed, quiet = FALSE) {
    if (length(object@kmeans$kpart) == 0) stop("run Clustexp before comptSNE")
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