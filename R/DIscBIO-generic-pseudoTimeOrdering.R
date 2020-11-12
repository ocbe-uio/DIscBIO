#' @title Pseudo-time ordering
#' @description This function takes the exact output of exprmclust function and
#'   construct Pseudo-time ordering by mapping all cells onto the path that
#'   connects cluster centers.
#' @param object \code{DISCBIO} class object.
#' @param quiet if `TRUE`, suppresses intermediary output
#' @param export if `TRUE`, exports order table to csv
#' @param filename Name of the exported file (if `export=TRUE`)
#' @importFrom TSCAN TSCANorder
#' @return The DISCBIO-class object input with the kordering slot filled.
setGeneric("pseudoTimeOrdering", function(
	object,
	quiet=FALSE,
	export=FALSE,
	filename="Cellular_pseudo-time_ordering"
	)
	{
		standardGeneric("pseudoTimeOrdering")
	}
)

#' @export
#' @rdname pseudoTimeOrdering
setMethod(
	"pseudoTimeOrdering",
	signature = "DISCBIO",
	definition = function(object, quiet, export, filename) {
		# ======================================================================
		# Validating
		# ======================================================================
		ran_k <- length(object@kmeans$kpart) > 0
		ran_m <- length(object@MBclusters) > 0
		if (ran_k) {
			Obj <- object@fdata
			Names <- object@cpart
			lpsmclust <- Exprmclust(Obj, K = 4, reduce = F, cluster = Names)
			lpsorder <- TSCANorder(lpsmclust)
		} else if (ran_m) {
			Obj <- object@fdata
			Names <- names(object@MBclusters$clusterid)
			lpsmclust <- object@MBclusters
			lpsorder <- TSCANorder(lpsmclust)
		} else {
			stop("run clustexp before this pseudoTimeOrdering")
		}
		# ======================================================================
		# Ordering
		# ======================================================================
		sampleNames <- colnames(Obj)
		orderID <- lpsorder
		order <- c(1:length(lpsorder))
		orderTable <- data.frame(order, orderID)
		if (export) write.csv(orderTable, file = paste0(filename, ".csv"))
		if (!quiet) print(orderTable)
		FinalOrder <- orderTable[match(sampleNames, orderTable$orderID), ]
		out_order <- FinalOrder[, 1]
		names(out_order) <- names(Names)
		if (ran_k) {
			object@kordering <- out_order
		} else if (ran_m) {
			object@MBordering <- out_order
		}
		return(object)
	}
)
