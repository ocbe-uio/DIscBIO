#' @title Plotting Gap Statistics
#' @export
#' @rdname plotGap
#' @param object \code{DISCBIO} class object.
#' @param y_limits 2-length numeric vector with the limits for the gap plot
#' @examples
#' sc <- DISCBIO(valuesG1msReduced) # changes signature of data
#' sc <- Clustexp(sc, cln=3) # data must be clustered before plotting
#' plotGap(sc)
setGeneric(
  "plotGap", function(object, y_limits = NULL) standardGeneric("plotGap")
)

#' @rdname plotGap
setMethod(
	f = "plotGap",
	signature = "DISCBIO",
	definition = function(object, y_limits) {
		if (length(object@kmeans$kpart) == 0) {
			stop("run clustexp before plotgap")
		}
		if (is.null(y_limits)) {
			y_lo <- min(x$Tab[, "gap"]) - 1 * max(x$Tab[, "SE.sim"])
			y_up <- max(x$Tab[, "gap"]) + 1 * max(x$Tab[, "SE.sim"])
			y_limits <- c(y_lo, y_up)
		}
		plot(
			object@kmeans$gap, las=1, ylim = y_limits, main="Gap Statistics"
		)
	}
)		 