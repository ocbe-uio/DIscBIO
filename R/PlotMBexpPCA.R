#' @title Plotting gene expression in Model-based clustering in PCA.
#' @description The PCA representation can also be used to show the gene expression of a particular gene.
#' @param object \code{DISCBIO} class object.
#' @param g  Individual gene name or vector with a group of gene names corresponding to a subset of valid row names of the \code{ndata} slot
#' of the \code{DISCBIO} object.
#' @param n String of characters representing the title of the plot. Default is NULL and the first element of \code{g} is chosen.
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics layout par image
#' @export

PlotMBexpPCA<- function(object,g,n=NULL) {
    if ( length(intersect(g,rownames(object@ndata))) < length(unique(g)) ) stop("second argument does not correspond to set of rownames slot ndata of SCseq object")
   	if (is.null(n))  n <- g[1]
    data=object@MBclusters
	#Expression<-cbind(data$pcareduceres,object@ndata[g,])
    l <- apply(object@ndata[g,] - .1,2,sum) + .1
    #l <- Expression[,3]
    mi <- min(l,na.rm=TRUE)
    ma <- max(l,na.rm=TRUE)
    ColorRamp <- colorRampPalette(rev(brewer.pal(n = 11,name = "RdYlBu")))(100)
    ColorLevels <- seq(mi, ma, length=length(ColorRamp))
    v <- round((l - mi)/(ma - mi)*99 + 1,0)
    layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
    par(mar = c(5,5,2.5,2))
    plot(data$pcareduceres[,1],data$pcareduceres[,2],xlab="PC1",ylab="PC2",pch=20,cex=0,col="grey",las=1,main=n)
    for ( k in 1:length(v) ){
        points(data$pcareduceres[k,1],data$pcareduceres[k,2],col=ColorRamp[v[k]],pch=20,cex=2)
    }
    par(mar = c(3,2.5,2.5,2))
    image(1, ColorLevels,matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),col=ColorRamp,xlab="",ylab="",las=2,xaxt="n")
    layout(1)
}