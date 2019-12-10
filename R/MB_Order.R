#' @title title
#' @description description 
#' @export
#' @param object object
#' @param sampleNames sampleNames
#' @param Names Names
#' @importFrom TSCAN TSCANorder
MB_Order<-function(object,sampleNames,Names){
	lpsorderMB <- TSCANorder(object)
	orderID<-lpsorderMB
	order<-c(1:length(lpsorderMB))
	orderTableMB<-data.frame(order,orderID)
	write.csv(orderTableMB, file = "Cellular_pseudo-time_ordering_based_on_Model-based_clusters.csv")
	print(orderTableMB)
	FinalOrder<-orderTableMB[match(sampleNames, orderTableMB$orderID),]
	sc@MBordering<-FinalOrder[,1]
	names(sc@MBordering)<-names(Names)
	return(sc@MBordering)
}