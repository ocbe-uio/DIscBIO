#' @title title
#' @export
#' @rdname comptSNE
setGeneric("comptSNE", function(object,rseed=15555) standardGeneric("comptSNE"))

#' @title title
#' @description description
#' @param object object
#' @param rseed rseed
#' @importFrom tsne tsne
#' @rdname comptSNE
#' @export
setMethod("comptSNE",
          signature = "PSCANseq",
          definition = function(object,rseed){
            if ( length(object@kmeans$kpart) == 0 ) stop("run clustexp before comptsne")
            set.seed(rseed)
            di <- dist.gen(as.matrix(object@distances))
            ts <- tsne(di,k=2)
            object@tsne <- as.data.frame(ts)
            return(object)
          }
          )