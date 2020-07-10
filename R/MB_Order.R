#' @title Pseudo-time ordering based on Model-based clusters
#' @description This function takes the exact output of exprmclust function and
#'   construct Pseudo-time ordering by mapping all cells onto the path that
#'   connects cluster centers.
#' @export
#' @param object \code{DISCBIO} class object.
#' @param quiet if `TRUE`, intermediary output is suppressed
#' @param export if `TRUE`, exports the results as a CSV file
#' @param filename Name of the exported file (if `export=TRUE`)
#' @importFrom TSCAN TSCANorder
#' @return The DISCBIO-class object input with the MBordering slot filled.
MB_Order <- function(
    object, quiet = FALSE, export = FALSE,
    filename = "Cellular_pseudo-time_ordering_based_on_Model-based_clusters"
) {
    data = object@MBclusters
    lpsorderMB <- TSCANorder(data)
    Names <- names(object@MBclusters$clusterid)
    sampleNames <- colnames(object@fdata)
    orderID <- lpsorderMB
    order <- c(1:length(lpsorderMB))
    orderTableMB <- data.frame(order, orderID)
    if (export) {
        write.csv(orderTableMB, file = paste0(filename, ".csv"))
    }
    if (!quiet) {
        print(orderTableMB)
    }
    FinalOrder <- orderTableMB[match(sampleNames, orderTableMB$orderID), ]
    MBordering <- FinalOrder[, 1]
    names(MBordering) <- names(Names)
    object@MBordering <- MBordering
    return(object)
}