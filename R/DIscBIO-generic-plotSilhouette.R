#' @title Silhouette Plot for K-means clustering
#' @description The silhouette provides a representation of how well each point
#'   is represented by its cluster in comparison to the closest neighboring
#'   cluster. It computes for each point the difference between the average
#'   similarity to all points in the same cluster and to all points in the
#'   closest neighboring cluster. This difference it normalize such that it can
#'   take values between -1 and 1 with higher values reflecting better
#'   representation of a point by its cluster.
#' @param object \code{DISCBIO} class object.
#' @param K A numeric value of the number of clusters
#' @importFrom stats as.dist cor
#' @importFrom cluster silhouette
#' @return A silhouette plot
setGeneric(
    name = "plotSilhouette",
    def = function(object, K) standardGeneric("plotSilhouette")
)

#' @export
#' @rdname plotSilhouette
setMethod(
    f = "plotSilhouette",
    signature = "DISCBIO",
    definition = function(object, K) {
        # ======================================================================
        # Validation
        # ======================================================================
        ran_clustexp <- length(object@kmeans$kpart) > 0
        ran_exprmclust <- length(object@MBclusters$clusterid) > 0
        if (ran_clustexp) {
            kpart <- object@kmeans$kpart
        } else if (ran_exprmclust) {
            kpart <- object@MBclusters$clusterid
        } else {
            stop("run clustexp or exprmclust before plotSilhouette")
        }
        if (length(unique(kpart)) < 2) {
            stop("only a single cluster: no silhouette plot")
        }
        # ======================================================================
        # Plotting
        # ======================================================================
        col <- c("black", "blue", "green", "red", "yellow", "gray")
        distances <- dist.gen(object@distances)
        si <- silhouette(kpart, distances)
        plot(si, col = col[1:K])
    }
)