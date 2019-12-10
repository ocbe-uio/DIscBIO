#' @title title
#' @description description 
#' @export
#' @param object object
#' @importFrom rpart rpart
#' @importFrom rpart.plot rpart.plot
RpartDT<-function(object){
	exp.df<-as.data.frame(t(object))
	classVector<- factor(colnames(object))
	model<-rpart(classVector~.,exp.df,method="class",minsplit = 1, cp=-1)
	print(model)
	rpart.plot(model,type=4,extra=101)
	return(model)
}
