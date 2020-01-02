#' @title title
#' @description description 
#' @param object \code{PSCANseq} class object.
#' @param types If types=NULL then the names of the cells will be grouped automatically. Default is NULL
#' @importFrom grDevices rainbow
#' @importFrom graphics legend
setGeneric("PCAplotSymbols", function(object,types=NULL) standardGeneric("PCAplotSymbols"))

#' @export
#' @rdname PCAplotSymbols
setMethod("PCAplotSymbols",
          signature = "PSCANseq",
          definition = function(object,types=NULL){
            if ( length(object@MBclusters) == 0 ) stop("run ExprmclustMB before PCAplotSymbols")
			total<-object@MBclusters
			types <- names(object@fdata)
			types<- gsub("_[0-9]+","",types)
			coloc <- rainbow(length(unique(types)))
			total<-object@MBclusters$pcareduceres
			syms<-c()
			plot(total[,1],total[,2],xlab="PC1",ylab="PC2",pch=20,cex=0,col="grey",las=1)
			for ( i in 1:length(unique(types)) ){
				f <- types == sort(unique(types))[i]
				syms <- append( syms, ( (i-1) %% 25 ) + 1 )
				points(total[f,1],total[f,2],col=coloc[i],pch=( (i-1) %% 25 ) + 1,cex=1)
			}
			legend("topright", legend=sort(unique(types)), col=coloc, pch=syms)
	} )      
