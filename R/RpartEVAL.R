RpartEVAL<- function(object,num.folds=10,First,Second){
	exp.imput.df<-as.data.frame(t(object))
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