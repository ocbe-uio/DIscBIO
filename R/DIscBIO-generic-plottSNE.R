#' @title tSNE map
#' @description Visualizing the k-means or model-based clusters using tSNE maps
#' @param object \code{DISCBIO} class object.
#' @importFrom graphics text
#' @return A plot of t-SNEs.
setGeneric("plottSNE", function(object)
	standardGeneric("plottSNE"))

#' @rdname plottSNE
#' @export
setMethod(
	"plottSNE",
	signature = "DISCBIO",
	definition = function(object) {
		# ======================================================================
		# Validating
		# ======================================================================
		ran_k <- length(object@tsne) > 0
		ran_m <- length(object@MBtsne) > 0
		if (ran_k) {
			part <- object@kmeans$kpart
			x <- object@tsne
		} else if (ran_m) {
			part <- object@MBclusters$clusterid
			x <- object@MBtsne
		} else {
			stop("run comptsne before plottSNE")
		}
		# ======================================================================
		# Plotting
		# ======================================================================
		col <- c("black", "blue", "green", "red", "yellow", "gray")
	LEN<-length(levels(factor(part)))
		plot(
			x,
			las = 1,
			xlab = "Dim 1",
			ylab = "Dim 2",
			pch = 20,
			cex = 1.5,
			col = "lightgrey"
		)
		for (i in seq_len(LEN)) {
			if (sum(part == i) > 0) {
				text(
					x[part == i, 1],
					x[part == i, 2],
					i,
					col = col[i],
					cex = .75,
					font = 4
				)
			}
		}
	}
)
