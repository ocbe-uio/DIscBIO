#' @title title
#' @description description 
#' @export
#' @param object object
#' @param Clusters Clusters
#' @param sampleNames sampleNames
#' @param quiet if `TRUE`, suppresses intermediary output
#' @param export if `TRUE`, exports order table to csv
#' @importFrom TSCAN TSCANorder
KmeanOrder<-function(object,Clusters,sampleNames, quiet = FALSE, export = TRUE){
	lpsmclust <- Exprmclust(object,clusternum =4,reduce = F, cluster = Clusters)
	lpsorder <- TSCANorder(lpsmclust)
	orderID<-lpsorder
	order<-c(1:length(lpsorder))
	orderTable<-data.frame(order,orderID)
	if (export) {
		write.csv(orderTable, file = "Cellular_pseudo-time_ordering_based_on_k-meansc-lusters.csv")
	}
	if (!quiet) {
		print(orderTable)
	}
	FinalOrder<-orderTable[match(sampleNames, orderTable$orderID),]
	sc@kordering<-FinalOrder[,1]
	names(sc@kordering)<-names(Clusters)
	return(sc@kordering)
}