#' @title title
#' @description description 
#' @export
#' @param object object
#' @param sampleNames sampleNames
#' @param Names Names
#' @param quiet if `TRUE`, intermediary output is suppressed
#' @param export if `TRUE`, exports the results as a CSV file
#' @importFrom TSCAN TSCANorder
MB_Order<-function(object, sampleNames,Names, quiet = FALSE, export = TRUE){
	lpsorderMB <- TSCANorder(object)
	orderID<-lpsorderMB
	order<-c(1:length(lpsorderMB))
	orderTableMB<-data.frame(order,orderID)
	if (export) {
		write.csv(
			orderTableMB,
			file = "Cellular_pseudo-time_ordering_based_on_Model-based_clusters.csv"
		)
	}
	if (!quiet) {
		print(orderTableMB)
	}
	FinalOrder<-orderTableMB[match(sampleNames, orderTableMB$orderID),]
	MBordering<-FinalOrder[,1]
	names(MBordering)<-names(Names)
	return(MBordering)
}