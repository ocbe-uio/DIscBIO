#' @title Pseudo-time ordering based on k-means clusters
#' @description This function takes the exact output of exprmclust function and
#'   construct Pseudo-time ordering by mapping all cells onto the path that
#'   connects cluster centers.
#' @param object \code{DISCBIO} class object.
#' @param quiet if `TRUE`, suppresses intermediary output
#' @param export if `TRUE`, exports order table to csv
#' @importFrom TSCAN TSCANorder
#' @return The DISCBIO-class object input with the kordering slot filled.
#' @examples
#' sc <- DISCBIO(valuesG1msReduced)
#' sc <- Clustexp(sc, cln=3) # K-means clustering
#' Order <- KmeanOrder(sc, export = FALSE)
#' Order@kordering
setGeneric("KmeanOrder", function(object, quiet = FALSE, export = TRUE)
    standardGeneric("KmeanOrder")
)

#' @export
#' @rdname KmeanOrder
setMethod(
    "KmeanOrder",
    signature = "DISCBIO",
    definition = function(object, quiet = FALSE, export = TRUE) {
        # Validation
        if (length(object@kmeans$kpart) == 0) {
            stop("run Clustexp before KmeanOrder")
        }

        Obj <- object@fdata
        Clusters <- object@cpart
        sampleNames <- colnames(object@fdata)
        lpsmclust <- Exprmclust(Obj, K = 4, reduce = F, cluster = Clusters)
        lpsorder <- TSCANorder(lpsmclust)
        orderID <- lpsorder
        order <- c(1:length(lpsorder))
        orderTable <- data.frame(order, orderID)
        if (export) {
            nm <- "Cellular_pseudo-time_ordering_based_on_k-meansc-lusters.csv"
            write.csv(orderTable, file = nm)
        }
        if (!quiet) {
            print(orderTable)
        }
        FinalOrder <-
            orderTable[match(sampleNames, orderTable$orderID), ]
        out_order <- FinalOrder[, 1]
        names(out_order) <- names(Clusters)
        object@kordering <- out_order
        return(object)
    }
)