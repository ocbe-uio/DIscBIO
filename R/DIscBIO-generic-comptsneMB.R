#' @title Computing tSNE for Model-based clustering
#' @description This function is used to compute the t-Distributed Stochastic
#'   Neighbor Embedding (t-SNE).
#' @param object \code{DISCBIO} class object.
#' @param max_iter maximum number of iterations to perform.
#' @param epoch The number of iterations in between update messages.
#' @param quiet if `TRUE`, suppresses intermediate output.
#' @param rseed Random number to be set for reproducible results.
#' @param ... other parameters to be passed to `tsne::tsne`.
#' @importFrom tsne tsne
#' @importFrom stats as.dist cor
#' @return The DISCBIO-class object input with the MBtsne slot filled.
#'
setGeneric(
    name = "comptsneMB",
    def = function(object, rseed = NULL, max_iter = 5000, epoch = 500,
        quiet = FALSE, ...) {
        standardGeneric("comptsneMB")
    }
)

#' @rdname comptsneMB
#' @export
setMethod(
    f = "comptsneMB",
    signature = "DISCBIO",
    definition = function(object, rseed, max_iter, epoch, quiet, ...) {
        if (length(object@MBclusters) == 0)
            stop("run clustexp before comptsneMB")
        set.seed(rseed)
        dist.gen <- function(x, method = "euclidean") {
            if (method %in% c("spearman", "pearson", "kendall")) {
                as.dist(1 - cor(t(x), method = method))
            } else {
                dist(x, method = method)
            }
        }
        di <- dist.gen(as.matrix(t(object@fdata)))
        if (quiet) {
            ts <- suppressMessages(
                tsne(di, max_iter = max_iter, epoch = epoch, ...)
            )
        } else {
            message("This function may take time")
            ts <- tsne(di, max_iter = max_iter, epoch = epoch, ...)
        }
        object@MBtsne <- as.data.frame(ts)
        return(object)
    }
)