J48DTeval<- function(object,num.folds=10,First,Second){
	exp.imput.df<-as.data.frame(t(object))
	num.instances<-nrow(exp.imput.df)
	indices<-1:num.instances
	classVector<- factor(colnames(object))

	cv.segments<-split(sample(indices),rep(1:num.folds,length=num.instances))
	j48.performance<-c("TP"=0,"FN"=0,"FP"=0,"TN"=0)
	j48.performance<-cross.val(exp.imput.df,classVector,cv.segments,j48.performance,"J48")
	print(j48.performance)

	j48.confusion.matrix<-matrix(j48.performance,nrow=2)
	rownames(j48.confusion.matrix)<-c(paste0("Predicted",First), paste0("Predicted",Second))
	colnames(j48.confusion.matrix)<-c(First,Second)
	print(j48.confusion.matrix)

	j48.sn<-SN(j48.confusion.matrix)
	j48.sp<-SP(j48.confusion.matrix)
	j48.acc<-ACC(j48.confusion.matrix)
	j48.mcc<-MCC(j48.confusion.matrix)

	cat("J48 SN: ", j48.sn, "\n",
		"J48 SP: ", j48.sp, "\n",
		"J48 ACC: ", j48.acc, "\n",
		"J48 MCC: ", j48.mcc, "\n",sep="")
	return(j48.performance)
}