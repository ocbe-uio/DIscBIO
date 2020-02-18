#' @title tSNE map for Model-based clustering with labels
#' @description Visualizing the Model-based clusters using tSNE maps
#' @param object \code{DISCBIO} class object.
#' @importFrom graphics text
#' @return A plot of the `object@MBtsne` values
#' @examples
#' sc <- DISCBIO(valuesG1msReduced)
#' sc <- NoiseFiltering(sc, percentile=0.9, CV=0.2, export=FALSE)
#' sc <- Normalizedata(
#'     sc, mintotal=1000, minexpr=0, minnumber=0, maxexpr=Inf, downsample=FALSE,
#'     dsn=1, rseed=17000
#' )
#' sc <- FinalPreprocessing(sc, GeneFlitering="NoiseF", export=FALSE)
#' sc <- Exprmclust(sc, K=2)
#' sc <- comptsneMB(sc, rseed=15555, quiet = TRUE)
#' plotMBLabelstSNE(sc)
setGeneric("plotMBLabelstSNE", function(object)
    standardGeneric("plotMBLabelstSNE"))

#' @rdname plotMBLabelstSNE
#' @export
setMethod(
    "plotMBLabelstSNE",
    signature = "DISCBIO",
    definition = function(object) {
        if (length(object@MBtsne) == 0)
            stop("run comptsneMB before plotMBLabelstSNE")
        Clusters <- object@MBclusters$clusterid
        ClustersFactor <- as.factor(Clusters)
        ClustersFactor <- gsub("1", "black", ClustersFactor)
        ClustersFactor <- gsub("2", "blue", ClustersFactor)
        ClustersFactor <- gsub("3", "green", ClustersFactor)
        ClustersFactor <- gsub("4", "red", ClustersFactor)
        ClustersFactor <- gsub("5", "yellow", ClustersFactor)
        ClustersFactor <- gsub("6", "gray", ClustersFactor)
        COL <- ClustersFactor
        labels = names(object@ndata)
        plot(
            object@MBtsne,
            xlab = "Dim 1",
            ylab = "Dim 2",
            pch = 20,
            cex = .5,
            col = "lightgrey"
        )
        text(object@MBtsne[, 1],
             object@MBtsne[, 2],
             labels,
             cex = .7,
             col = COL)
    }
) 