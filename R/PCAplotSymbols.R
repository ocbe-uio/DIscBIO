#' @title title
#' @description description 
#' @export
#' @param object1 object1
#' @param object2 object2
#' @param types types
#' @importFrom grDevices rainbow
#' @importFrom graphics legend
PCAplotSymbols= function(object1,object2,types=NULL){
	types <- names(object2)
	types<- gsub("_[0-9]+","",types)
	coloc <- rainbow(length(unique(types)))
	syms<-c()
	plot(object1[,1],object1[,2],xlab="PC1",ylab="PC2",pch=20,cex=0,col="grey",las=1)
    for ( i in 1:length(unique(types)) ){
		f <- types == sort(unique(types))[i]
        syms <- append( syms, ( (i-1) %% 25 ) + 1 )
		points(object1[f,1],object1[f,2],col=coloc[i],pch=( (i-1) %% 25 ) + 1,cex=1)
        }
    legend("topright", legend=sort(unique(types)), col=coloc, pch=syms)
}       
