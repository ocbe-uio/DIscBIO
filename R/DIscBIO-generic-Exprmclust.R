#' @title Performing Model-based clustering on expression values
#' @description this function first uses principal component analysis (PCA) to
#'   reduce dimensionality of original data. It then performs model-based
#'   clustering on the transformed expression values.
#' @param object \code{DISCBIO} class object.
#' @param K An integer vector specifying all possible cluster numbers. Default
#'   is 3.
#' @param modelNames model to be used in model-based clustering. By default
#'   "ellipsoidal, varying volume, shape, and orientation" is used.
#' @param reduce A logical vector that allows performing the PCA on the
#'   expression data. Default is TRUE.
#' @param cluster A vector showing the ID of cells in the clusters.
#' @param quiet if `TRUE`, suppresses intermediary output
#' @importFrom mclust Mclust mclustBIC
#' @importFrom stats dist prcomp lm
#' @importFrom igraph graph.adjacency minimum.spanning.tree
#' @return If `object` is of class DISCBIO, the output is the same object  with
#'   the MBclusters slot filled. If the `object` is a data frame, the function
#'   returns a named list containing the four objects that together correspond
#'   to the contents of the MBclusters slot.
#'
setGeneric(
  name = "Exprmclust",
  def = function(
      object,
      K = 3,
      modelNames = "VVV",
      reduce = TRUE,
      cluster = NULL,
      quiet = FALSE) {
    standardGeneric("Exprmclust")
  }
)

#' @export
#' @rdname Exprmclust
setMethod(
  f = "Exprmclust",
  signature = "DISCBIO",
  definition = function(object, K, modelNames, reduce, cluster, quiet) {
    pcareduceres <- calc_pcareduceres(object@fdata, reduce)
    if (is.null(cluster)) {
      K <- K[K > 1]
      res <- Mclust(
        data = pcareduceres,
        G = K,
        modelNames = modelNames,
        warn = FALSE,
        verbose = !quiet
      )
      if (is.null(res)) stop("Unable to cluster. Try a lower value for K.")
      clusterid <- apply(res$z, 1, which.max)
      clunum <- res$G
    } else {
      clunum <- length(unique(cluster))
      clusterid <- cluster
    }
    clucenter <- matrix(0, ncol = ncol(pcareduceres), nrow = clunum)
    for (cid in 1:clunum) {
      clucenter[cid, ] <- colMeans(
        pcareduceres[names(clusterid[clusterid == cid]), , drop = FALSE]
      )
    }
    dp <- as.matrix(dist(clucenter))
    gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
    dp_mst <- minimum.spanning.tree(gp)
    full_List <- list(
      pcareduceres = pcareduceres,
      MSTtree = dp_mst,
      clusterid = clusterid,
      clucenter = clucenter
    )
    object@MBclusters <- full_List
    return(object)
  }
)

#' @export
#' @rdname Exprmclust
setMethod(
  f = "Exprmclust",
  signature = "data.frame",
  definition = function(object,
                        K = 3,
                        modelNames = "VVV",
                        reduce = TRUE,
                        cluster = NULL,
                        quiet = FALSE) {
    pcareduceres <- calc_pcareduceres(object, reduce)
    if (is.null(cluster)) {
      K <- K[K > 1]
      res <- Mclust(
        data = pcareduceres,
        G = K,
        modelNames = modelNames,
        warn = FALSE,
        verbose = !quiet
      )
      if (is.null(res)) {
        stop("Unable to cluster. Try a lower value for K.")
      }
      clusterid <- apply(res$z, 1, which.max)
      clunum <- res$G
    } else {
      clunum <- length(unique(cluster))
      clusterid <- cluster
    }
    clucenter <- matrix(0, ncol = ncol(pcareduceres), nrow = clunum)
    for (cid in 1:clunum) {
      clucenter[cid, ] <- colMeans(
        pcareduceres[names(clusterid[clusterid == cid]), , drop = FALSE]
      )
    }
    dp <- as.matrix(dist(clucenter))
    gp <-
      graph.adjacency(dp, mode = "undirected", weighted = TRUE)
    dp_mst <- minimum.spanning.tree(gp)
    object <- list(
      pcareduceres = pcareduceres,
      MSTtree = dp_mst,
      clusterid = clusterid,
      clucenter = clucenter
    )
    return(object)
  }
)
