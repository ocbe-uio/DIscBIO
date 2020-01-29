#' @title Highlighting gene expression in K-means clustering in the t-SNE map
#' @description The t-SNE map representation can also be used to analyze expression of a gene or a group of genes, 
#' to investigate cluster specific gene expression patterns
#' @param object \code{DISCBIO} class object.
#' @param g  Individual gene name or vector with a group of gene names corresponding to a subset of valid row names of the \code{ndata} slot
#' of the \code{DISCBIO} object.
#' @param n String of characters representing the title of the plot. Default is NULL and the first element of \code{g} is chosen.
#' @examples 
#' sc <- DISCBIO(valuesG1ms)
#' sc <- Clustexp(sc, cln=3, quiet=TRUE) # K-means clustering
#' sc <- comptSNE(sc, rseed=15555, quiet=TRUE)
#' g <- 'ENSG00000001460'
#' plotExptSNE(sc, g)

setGeneric("plotExptSNE", function(object,g,n=NULL) standardGeneric("plotExptSNE"))

#' @export
#' @rdname plotExptSNE
setMethod("plotExptSNE",
          signature = "DISCBIO",
          definition = function(object,g,n=NULL){
            if ( length(object@tsne) == 0 ) stop("run comptSNE before plotExptSNE")
            if ( length(intersect(g,rownames(object@ndata))) < length(unique(g)) ) stop("second argument does not correspond to set of rownames slot ndata of SCseq object")
            if (is.null(n))  n <- g[1]
            l <- apply(object@ndata[g,] - .1,2,sum) + .1
            mi <- min(l,na.rm=TRUE)
            ma <- max(l,na.rm=TRUE)
            ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7,name = "RdYlBu")))(100)
            ColorLevels <- seq(mi, ma, length=length(ColorRamp))
            v <- round((l - mi)/(ma - mi)*99 + 1,0)
            layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
            par(mar = c(3,5,2.5,2))
            plot(object@tsne,xlab="Dim 1",ylab="Dim 2",main=n,pch=20,cex=0,col="grey",las=1)
            for ( k in 1:length(v) ){
              points(object@tsne[k,1],object@tsne[k,2],col=ColorRamp[v[k]],pch=20,cex=1.5)
            }
            par(mar = c(3,2.5,2.5,2))
            image(1, ColorLevels,
                  matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
                  col=ColorRamp,
                  xlab="",ylab="",las=1,
                  xaxt="n")
            layout(1)
          }
          )