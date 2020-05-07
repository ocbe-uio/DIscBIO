#' @title Computing tSNE for Model-based clustering
#' @description This function is used to compute the t-Distributed Stochastic
#'   Neighbor Embedding (t-SNE).
#' @param object \code{DISCBIO} class object.
#' @param rseed Integer number. Random seed to to yield exactly reproducible
#'   maps across different runs. Default is 15555.
#' @param max_iter maximum number of iterations to perform.
#' @param epoch The number of iterations in between update messages.
#' @param quiet if `TRUE`, suppresses intermediate output
#' @param ... other parameters to be passed to `tsne::tsne`
#' @importFrom tsne tsne
#' @importFrom stats as.dist cor
#' @return The DISCBIO-class object input with the MBtsne slot filled.
#' @examples
#' sc <- DISCBIO(valuesG1msReduced)
#' sc <- NoiseFiltering(sc, percentile=0.9, CV=0.2, export=FALSE)
#' sc <- Normalizedata(
#'     sc, mintotal=1000, minexpr=0, minnumber=0, maxexpr=Inf, downsample=FALSE,
#'     dsn=1, rseed=17000
#' )
#' sc <- FinalPreprocessing(sc, GeneFlitering="NoiseF", export=FALSE)
#' sc <- Exprmclust(sc)
#' sc <- comptsneMB(sc, rseed=15555, max_iter = 1000)
#' print(sc@MBtsne)
setGeneric(
    name = "comptsneMB",
    def = function(object, rseed = 15555, max_iter = 5000, epoch = 500,
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