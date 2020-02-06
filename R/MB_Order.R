#' @title Pseudo-time ordering based on Model-based clusters
#' @description This function takes the exact output of exprmclust function and construct Pseudo-time ordering by mapping all cells onto the path that connects cluster centers. 
#' @export
#' @param object \code{DISCBIO} class object.
#' @param quiet if `TRUE`, intermediary output is suppressed
#' @param export if `TRUE`, exports the results as a CSV file
#' @importFrom TSCAN TSCANorder
#' @examples
#' sc <- DISCBIO(valuesG1ms)
#' sc <- NoiseFiltering(sc, export=FALSE)
#' sc <- Normalizedata(
#'     sc, mintotal=1000, minexpr=0, minnumber=0, maxexpr=Inf, downsample=FALSE,
#'     dsn=1, rseed=17000
#' )
#' sc <- FinalPreprocessing(sc, GeneFlitering="NoiseF", export=FALSE)
#' sc <- Exprmclust(sc)
#' sc <- comptsneMB(sc, rseed=15555)
#' sc <- Clustexp(sc, cln=3)
#' sc <- MB_Order(sc, export = FALSE)
#' sc@MBordering
MB_Order<-function(object, quiet = FALSE, export = TRUE){
	data=object@MBclusters
	lpsorderMB <- TSCANorder(data)
	Names<-names(object@MBclusters$clusterid)
	sampleNames<-colnames(object@fdata)
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
	object@MBordering<-MBordering
	return(object)
}