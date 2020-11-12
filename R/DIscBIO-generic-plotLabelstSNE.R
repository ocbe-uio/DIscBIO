#' @title tSNE map with labels
#' @description Visualizing k-means or model-based clusters using tSNE maps
#' @param object \code{DISCBIO} class object.
#' @rdname plotLabelstSNE
#' @importFrom graphics text
#' @return Plot containing the ID of the cells in each cluster
setGeneric("plotLabelstSNE", function(object)
	{
	standardGeneric("plotLabelstSNE")
	}
)

#' @rdname plotLabelstSNE
#' @export
setMethod(
	"plotLabelstSNE",
	signature = "DISCBIO",
	definition = function(object)
	{
		# ======================================================================
		# Validating
		# ======================================================================
		ran_k <- length(object@tsne) > 0
		ran_m <- length(object@MBtsne) > 0
		if (ran_k) {
			Clusters <- object@kmeans$kpart
			x <- object@tsne
		} else if (ran_m) {
			Clusters <- object@MBclusters$clusterid
			x <- object@MBtsne
		} else {
			stop("run comptsne before plotLabelstSNE")
		}
		# ======================================================================
		# Plotting
		# ======================================================================
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
			x,
			xlab = "Dim 1",
			ylab = "Dim 2",
			pch = 20,
			cex = .5,
			col = "lightgrey"
		)
		text(
			x[, 1],
			x[, 2],
			labels,
			cex = .7,
			col = COL
		)
	}
)