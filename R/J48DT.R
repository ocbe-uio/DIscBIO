#' @title title
#' @description description 
#' @export
#' @param object object
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