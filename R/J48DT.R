#' @title J48 Decision Tree
#' @description The decision tree analysis is implemented over a training dataset, which consisted of the DEGs obtained by either SAMseq or the binomial differential expression.  
#' @export
#' @param data A data frame resulted from running the function ClassVectoringDT.
#' @importFrom RWeka J48
#' @importFrom graphics plot
#' @importFrom partykit as.party
#' @importFrom grid gpar
J48DT<-function(data){
		msg <- NULL
		if ( ! is.data.frame(data) ){
            msg <- c(msg, "input data must be data.frame")
        }else if ( nrow(data) < 2 ){
            msg <- c(msg, "input data must have more than one row")
        }else if ( ncol(data) < 2 ){
            msg <- c(msg, "input data must have more than one column")
        }else if (sum( apply( is.na(data),1,sum ) ) > 0 ){
            msg <- c(msg, "NAs are not allowed in input data")
        }else if (sum( apply( data,1,min ) ) < 0 ){
            msg <- c(msg, "negative values are not allowed in input data")
        }
        if (is.null(msg)) TRUE
        else msg
		
		exp.df<-as.data.frame(t(data))
		classVector<- factor(colnames(data))
		j48.model<-J48(classVector~.,exp.df)
		print(j48.model)
		plot(as.party(j48.model),gp = gpar(cex=0.65,col="black", lty = "solid", lwd = 1.5, fontsize = 12))
		return(j48.model)
}