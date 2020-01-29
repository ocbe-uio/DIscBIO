#' @title Pseudo-time ordering based on k-means clusters
#' @description This function takes the exact output of exprmclust function and construct Pseudo-time ordering by mapping all cells onto the path that connects cluster centers. 
#' @param object \code{DISCBIO} class object.
#' @param quiet if `TRUE`, suppresses intermediary output
#' @param export if `TRUE`, exports order table to csv
#' @importFrom TSCAN TSCANorder
setGeneric("KmeanOrder", function(object, quiet = FALSE, export = TRUE) standardGeneric("KmeanOrder"))

#' @export
#' @rdname KmeanOrder
setMethod("KmeanOrder",
          signature = "DISCBIO",
          definition = function(object, quiet = FALSE, export = TRUE) {
			Obj <- object@fdata
			Clusters<-object@cpart
			sampleNames<-colnames(object@fdata)
			lpsmclust <- Exprmclust(Obj, clusternum = 4, reduce = F, cluster = Clusters)
			lpsorder <- TSCANorder(lpsmclust)
			orderID <- lpsorder
			order <- c(1:length(lpsorder))
			orderTable <- data.frame(order, orderID)
			if (export) {
				write.csv(orderTable,file = "Cellular_pseudo-time_ordering_based_on_k-meansc-lusters.csv")
			}
			if (!quiet) {
				print(orderTable)
			}
			FinalOrder <- orderTable[match(sampleNames, orderTable$orderID), ]
			out_order <- FinalOrder[, 1]
			names(out_order) <- names(Clusters)
			object@kordering<-out_order
			return(object)
})