#' @title RPART Decision Tree
#' @description The decision tree analysis is implemented over a training dataset, which consisted of the DEGs obtained by either SAMseq or the binomial differential expression.  
#' @export
#' @importFrom rpart rpart
#' @importFrom rpart.plot rpart.plot
RpartDT<-function(data){
	exp.df<-as.data.frame(t(data))
	classVector<- factor(colnames(data))
	model<-rpart(classVector~.,exp.df,method="class",minsplit = 1, cp=-1)
	print(model)
	rpart.plot(model,type=4,extra=101)
	return(model)
}


