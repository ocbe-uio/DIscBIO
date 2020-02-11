#' @title tSNE map for Model-based clustering
#' @description Visualizing the Model-based clusters using tSNE maps
#' @param object \code{DISCBIO} class object.
#' @param K A numeric value of the number of clusters
#' @importFrom graphics text
#' @return A plot of t-SNEs.
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
#' plottsneMB(sc)
setGeneric(
	name = "plottsneMB",
	def = function(object, K = length(table(object@MBclusters$clusterid))) {
		standardGeneric("plottsneMB")		
	}
)

#' @export
#' @rdname plottsneMB
setMethod(
	f = "plottsneMB",
	signature = "DISCBIO",
	definition = function(object, K) {
		if (length(object@MBtsne) == 0) stop("run comptsneMB before plottsneMB")
		col = c("black", "blue", "green", "red", "yellow", "gray")
		part <- object@MBclusters$clusterid
		plot(
			object@MBtsne, xlab="Dim 1", ylab="Dim 2", pch=20, cex=1.5,
			col="lightgrey", las=1
		)
		for (i in 1:K) {
			if (sum(part == i) > 0) {
				text(
					object@MBtsne[part == i, 1], object@MBtsne[part == i, 2], i,
					col=col[i], cex=.75, font=4
				)
			}
		}
	}
)    