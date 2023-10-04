#' @title Computing tSNE
#' @description This function is used to compute the t-Distributed Stochastic
#'   Neighbor Embedding (t-SNE).
#' @param object \code{DISCBIO} class object.
#' @param rseed Random integer to to yield reproducible maps across different
#' runs
#' @param max_iter maximum number of iterations to perform.
#' @param epoch The number of iterations in between update messages.
#' @param quiet if `TRUE`, suppresses intermediate output
#' @param ... other parameters to be passed to `tsne::tsne`
#' @importFrom tsne tsne
#' @importFrom stats as.dist cor
#' @return The DISCBIO-class object input with the tsne slot filled.
#' @examples
#' sc <- DISCBIO(valuesG1msTest) # changes signature of data
#' sc <- Clustexp(sc, cln = 2) # data must be clustered before plottin
#' sc <- comptSNE(sc, max_iter = 30)
#' head(sc@tsne)
#'
setGeneric(
  name = "comptSNE",
  def = function(
      object, rseed = NULL, max_iter = 5000, epoch = 500, quiet = FALSE, ...) {
    standardGeneric("comptSNE")
  }
)

#' @rdname comptSNE
#' @export
setMethod(
  f = "comptSNE",
  signature = "DISCBIO",
  definition = function(object, rseed, max_iter, epoch, quiet, ...) {
    # ======================================================================
    # Validating
    # ======================================================================
    ran_k <- length(object@kmeans$kpart) > 0
    ran_m <- length(object@MBclusters) > 0
    if (ran_k) {
      di <- dist.gen(as.matrix(object@distances))
    } else if (ran_m) {
      di <- dist.gen(as.matrix(t(object@fdata)))
    } else {
      stop("run clustexp before comptSNE")
    }
    # ======================================================================
    # Computing
    # ======================================================================
    set.seed(rseed)
    if (quiet) {
      ts <- suppressMessages(
        tsne(di, max_iter = max_iter, epoch = epoch, ...)
      )
    } else {
      message("This function may take time")
      ts <- tsne(di, max_iter = max_iter, epoch = epoch, ...)
    }
    # ======================================================================
    # Filling output
    # ======================================================================
    if (ran_k) {
      object@tsne <- as.data.frame(ts)
    } else if (ran_m) {
      object@MBtsne <- as.data.frame(ts)
    }
    return(object)
  }
)
