#' @title tSNE map for K-means clustering with labels
#' @description Visualizing the K-means clusters using tSNE maps
#' @param object \code{DISCBIO} class object.
#' @rdname plotKmeansLabelstSNE
#' @importFrom graphics text
#' @return Plot containing the ID of the cells in each cluster
#' @examples
#' sc <- DISCBIO(valuesG1msReduced) # changes signature of data
#' sc <- Clustexp(sc, cln=3) # data must be clustered before plottin
#' sc <- comptSNE(sc, rseed=15555, quiet=TRUE)
#' plotKmeansLabelstSNE(sc) # Plots the ID of the cells in each cluster
setGeneric("plotKmeansLabelstSNE", function(object)
    standardGeneric("plotKmeansLabelstSNE"))

#' @rdname plotKmeansLabelstSNE
#' @export
setMethod(
    "plotKmeansLabelstSNE",
    signature = "DISCBIO",
    definition = function(object) {
        if (length(object@tsne) == 0)
            stop("run comptsne before plotKmeansLabelstSNE")
        Clusters <- object@kmeans$kpart
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
            object@tsne,
            xlab = "Dim 1",
            ylab = "Dim 2",
            pch = 20,
            cex = .5,
            col = "lightgrey"
        )
        text(object@tsne[, 1],
             object@tsne[, 2],
             labels,
             cex = .7,
             col = COL)
    }
)