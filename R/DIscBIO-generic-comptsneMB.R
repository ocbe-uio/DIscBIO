#' @title Computing tSNE for Model-based clustering
#' @description This function is used to compute the t-Distributed Stochastic Neighbor Embedding (t-SNE).
#' @param object \code{DISCBIO} class object.
#' @param rseed Integer number. Random seed to to yield exactly reproducible maps across different runs. Default is 15555. 
#' @param quiet if `TRUE`, suppresses intermediate output
#' @importFrom tsne tsne
#' @importFrom stats as.dist cor
#' @examples
#' sc <- DISCBIO(valuesG1ms)
#' sc <- NoiseFiltering(sc, export=FALSE)
#' sc <- Normalizedata(
#'     sc, mintotal=1000, minexpr=0, minnumber=0, maxexpr=Inf, downsample=FALSE,
#'     dsn=1, rseed=17000
#' )
#' sc <- FinalPreprocessing(sc, GeneFlitering="NoiseF", export=FALSE)
#' sc <- Exprmclust(sc, K = 3,reduce = T)
#' sc <- comptsneMB(sc, rseed=15555, quiet = TRUE)
#' print(sc@MBtsne)
setGeneric(
  name = "comptsneMB",
  def = function(object, rseed = 15555, quiet = FALSE) {
    standardGeneric("comptsneMB")
  }
)

#' @rdname comptsneMB
#' @export
setMethod(
  f = "comptsneMB",
  signature = "DISCBIO",
  definition = function(object, rseed, quiet = FALSE) {
    if (length(object@MBclusters) == 0) stop("run clustexp before comptsneMB")
    set.seed(rseed)
    dist.gen <- function(x, method = "euclidean", ...) {
      if (method %in% c("spearman","pearson","kendall")) {
        as.dist(1 - cor(t(x), method = method, ...))
      } else {
        dist(x,method=method,...)
      }
    }
    di <- dist.gen(as.matrix(t(object@fdata)))
    if (quiet) {
      ts <- suppressMessages(tsne(di,k=2,max_iter = 5000,epoch=500))
    } else {
      cat("This function takes time")
      ts <- tsne(di,k=2,max_iter = 5000,epoch=500)
    }
    object@MBtsne <- as.data.frame(ts)
    return(object)
})