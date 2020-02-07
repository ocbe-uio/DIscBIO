#' @title RPART Decision Tree
#' @description The decision tree analysis is implemented over a training dataset, which consisted of the DEGs obtained by either SAMseq or the binomial differential expression.  
#' @param data The exact output of the exprmclust function.
#' @param quiet If `TRUE`, suppresses intermediary output
#' @param plot If `FALSE`, suppresses plot output
#' @export
#' @importFrom rpart rpart
#' @importFrom rpart.plot rpart.plot
#' @examples
#' \dontrun{
#' sc <- DISCBIO(valuesG1msReduced)
#' sc <- NoiseFiltering(sc, export=FALSE)
#' sc <- Normalizedata(
#'     sc, mintotal=1000, minexpr=0, minnumber=0, maxexpr=Inf, downsample=FALSE,
#'     dsn=1, rseed=17000
#' )
#' sc <- FinalPreprocessing(sc, GeneFlitering="NoiseF", export=FALSE)
#' sc <- Clustexp(sc, cln=3) # K-means clustering
#' sc <- comptSNE(sc, rseed=15555)
#' cdiff <- DEGanalysis2clust(
#'     sc, Clustering="K-means", K=3, fdr=.2, name="Name", First="CL1",
#'     Second="CL2", export=FALSE
#' )
#' sigDEG <- cdiff[[1]]
#' DATAforDT <- ClassVectoringDT(
#'     sc, Clustering="K-means", K=3, First="CL1", Second="CL2", sigDEG,
#' )
#' RpartDT(DATAforDT)
#' }
RpartDT<-function(data, quiet = FALSE, plot = TRUE){
	exp.df<-as.data.frame(t(data))
	classVector<- factor(colnames(data))
	model<-rpart(classVector~.,exp.df,method="class",minsplit = 1, cp=-1)
	if (!quiet) print(model)
	if (plot) rpart.plot(model,type=4,extra=101)
	return(model)
}


