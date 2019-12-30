#' @title J48 Decision Tree
#' @description The decision tree analysis is implemented over a training dataset, which consisted of the DEGs obtained by either SAMseq or the binomial differential expression.  
#' @export
#' @param data A data frame resulted from running the function ClassVectoringDT.
#' @importFrom RWeka J48
#' @importFrom graphics plot
#' @importFrom partykit as.party
#' @importFrom grid gpar
J48DT<-function(object){
	exp.df<-as.data.frame(t(object))
	classVector<- factor(colnames(object))
	j48.model<-J48(classVector~.,exp.df)
	print(j48.model)
      plot(as.party(j48.model),gp = gpar(cex=0.65,col="black", lty = "solid", lwd = 1.5, fontsize = 12))
	return(j48.model)
}