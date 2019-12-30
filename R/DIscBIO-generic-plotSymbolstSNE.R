#' @title tSNE map for K-means clustering with symbols
#' @description Visualizing the K-means clusters using tSNE maps 
#' @param object \code{PSCANseq} class object.
#' @param types If types=NULL then the names of the cells will be grouped automatically. Default is NULL
#' @export
#' @rdname plotSymbolstSNE
setGeneric("plotSymbolstSNE", function(object,types=NULL) standardGeneric("plotSymbolstSNE"))

setMethod("plotSymbolstSNE",
          signature = "PSCANseq",
          definition = function(object,types){
            if ( is.null(types) ) types <- names(object@fdata)
            if ( length(object@tsne) == 0 ) stop("run comptsne before plotSymbolstSNE")
            if ( length(types) != ncol(object@fdata) ) stop("types argument has wrong length. Length has to equal to the column number of object@ndata")
			
            coloc <- rainbow(length(unique(types)))
            syms <- c()
            plot(object@tsne,xlab="Dim 1",ylab="Dim 2",pch=20,col="grey")
            for ( i in 1:length(unique(types)) ){
              f <- types == sort(unique(types))[i]
              syms <- append( syms, ( (i-1) %% 25 ) + 1 )
              points(object@tsne[f,1],object@tsne[f,2],col=coloc[i],pch=( (i-1) %% 25 ) + 1,cex=1)
            }
            legend("topright", legend=sort(unique(types)), col=coloc, pch=syms)
          }
          )