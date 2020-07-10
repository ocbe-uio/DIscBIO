#' @title Plotting Gap Statistics
#' @export
#' @rdname plotGap
#' @param object \code{DISCBIO} class object.
#' @param y_limits 2-length numeric vector with the limits for the gap plot
#' @return A plot of the gap statistics
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
        gap <- object@kmeans$gap
        if (is.null(y_limits)) {
            y_lo <- min(gap$Tab[, "gap"]) - 1 * max(gap$Tab[, "SE.sim"])
            y_up <- max(gap$Tab[, "gap"]) + 1 * max(gap$Tab[, "SE.sim"])
            y_limits <- c(y_lo, y_up)
        }
        plot(
            gap, las=1, ylim = y_limits, main="Gap Statistics"
        )
    }
)