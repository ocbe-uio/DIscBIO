#' @title title
#' @description description 
#' @export
#' @param object object
#' @param Clusters Clusters
#' @param sampleNames sampleNames
#' @importFrom TSCAN TSCANorder
KmeanOrder<-function(object,Clusters,sampleNames){
	lpsmclust <- Exprmclust(object,clusternum =4,reduce = F, cluster = Clusters)
	lpsorder <- TSCANorder(lpsmclust)
	orderID<-lpsorder
	order<-c(1:length(lpsorder))
	orderTable<-data.frame(order,orderID)
	write.csv(orderTable, file = "Cellular_pseudo-time_ordering_based_on_k-meansc-lusters.csv")
	print(orderTable)
	FinalOrder<-orderTable[match(sampleNames, orderTable$orderID),]
	sc@kordering<-FinalOrder[,1]
	names(sc@kordering)<-names(Clusters)
	return(sc@kordering)
}