#' @title title
#' @export
#' @rdname comptSNE
setGeneric(
  name = "comptSNE",
  def = function(object, rseed = 15555, quiet = FALSE) {
    standardGeneric("comptSNE")
  }
)

#' @title title
#' @description description
#' @param object object
#' @param rseed rseed
#' @param quiet if `TRUE`, suppresses intermediate output
#' @importFrom tsne tsne
#' @rdname comptSNE
#' @export
setMethod(
  "comptSNE",
  signature = "PSCANseq",
  definition = function(object, rseed, quiet = FALSE) {
    if (length(object@kmeans$kpart) == 0) stop("run clustexp before comptsne")
    set.seed(rseed)
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