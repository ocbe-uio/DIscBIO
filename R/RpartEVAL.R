#' @title Evaluating the performance of the RPART Decision Tree.
#' @description This function evaluates the performance of the generated trees for error estimation by ten-fold cross validation assessment.
#' @export
#' @param data The resulted data from running the function J48DT.
#' @param num.folds A numeric value of the number of folds for the cross validation assessment. Default is 10.
#' @param First A string vector showing the first target cluster.  Default is "CL1"
#' @param Second A string vector showing the second target cluster.  Default is "CL2"

	num.instances<-nrow(exp.imput.df)
	indices<-1:num.instances
	classVector<- factor(colnames(object))

	cv.segments<-split(sample(indices),rep(1:num.folds,length=num.instances))
	Rpart.performance<-c("TP"=0,"FN"=0,"FP"=0,"TN"=0)
	Rpart.performance<-cross.val(exp.imput.df,classVector,cv.segments,Rpart.performance,"rpart")
	print(Rpart.performance)
	Rpart.confusion.matrix<-matrix(Rpart.performance,nrow=2)
	rownames(Rpart.confusion.matrix)<-c(paste0("Predicted",First), paste0("Predicted",Second))
	colnames(Rpart.confusion.matrix)<-c(First,Second)
	print(Rpart.confusion.matrix)

	Rpart.sn<-SN(Rpart.confusion.matrix)
	Rpart.sp<-SP(Rpart.confusion.matrix)
	Rpart.acc<-ACC(Rpart.confusion.matrix)
	Rpart.mcc<-MCC(Rpart.confusion.matrix)

	cat("Rpart SN: ", Rpart.sn, "\n",
	"Rpart SP: ", Rpart.sp, "\n",
	"Rpart ACC: ", Rpart.acc, "\n",
	"Rpart MCC: ", Rpart.mcc, "\n",sep="")

	return(Rpart.performance)
}