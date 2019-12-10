#' @title title
#' @export
#' @rdname comptsneMB
setGeneric("comptsneMB", function(object,rseed=15555) standardGeneric("comptsneMB"))
#' @title title
#' @description description
#' @param object object
#' @param rseed rseed
#' @importFrom tsne tsne
#' @rdname comptsneMB
#' @export
setMethod("comptsneMB",
          signature = "PSCANseq",
          definition = function(object,rseed){
            if ( length(object@MBclusters) == 0 ) stop("run clustexp before comptsneMB")
            set.seed(rseed)
            dist.gen <- function(x,method="euclidean", ...) if ( method %in% c("spearman","pearson","kendall") ) as.dist( 1 - cor(t(x),method=method,...) ) else dist(x,method=method,...)
            di <- dist.gen(as.matrix(t(procdataTSCAN)))
            cat("This function takes time")
            ts <- tsne(di,k=2,max_iter = 5000,epoch=500)
            object@MBtsne <- as.data.frame(ts)
            return(object)
          }
          )