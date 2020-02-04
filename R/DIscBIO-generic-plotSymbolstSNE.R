#' @title tSNE map for K-means clustering with symbols
#' @description Visualizing the K-means clusters using tSNE maps 
#' @param object \code{DISCBIO} class object.
#' @param types If types=NULL then the names of the cells will be grouped automatically. Default is NULL
#' @param legloc A keyword from the list "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center". Default is "bottomright"
#' @examples 
#' sc <- DISCBIO(valuesG1ms) # changes signature of data
#' sc <- Clustexp(sc, cln=3) # data must be clustered before plottin
#' sc <- comptSNE(sc, rseed=15555, quiet=TRUE)
#' plotSymbolstSNE(sc,types=sub("(\\_\\d+)$","", names(sc@ndata)))    
setGeneric("plotSymbolstSNE", function(object,types=NULL,legloc="bottomright") standardGeneric("plotSymbolstSNE"))

#' @export
#' @rdname plotSymbolstSNE
setMethod("plotSymbolstSNE",
          signature = "DISCBIO",
          definition = function(object,types,legloc="bottomright"){
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
            legend(legloc, legend=sort(unique(types)), col=coloc, pch=syms)
          }
          )