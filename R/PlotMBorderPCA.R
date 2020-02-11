#' @title Plotting pseudo-time ordering in Model-based clustering in PCA.
#' @description The PCA representation can also be used to show the pseudo-time ordering.
#' @param object \code{DISCBIO} class object.
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics layout par image
#' @export
#' @return A plot of the PCA.
#' @examples
#' \dontrun{
#' sc <- DISCBIO(valuesG1msReduced)
#' sc <- NoiseFiltering(sc, percentile=0.9, CV=0.2, export=FALSE)
#' sc <- Normalizedata(
#'     sc, mintotal=1000, minexpr=0, minnumber=0, maxexpr=Inf, downsample=FALSE,
#'     dsn=1, rseed=17000
#' )
#' sc <- FinalPreprocessing(sc, GeneFlitering="NoiseF", export=FALSE)
#' sc <- Exprmclust(sc, K=3)
#' sc <- comptsneMB(sc, rseed=15555)
#' sc <- Clustexp(sc, cln=3)
#' sc <- MB_Order(sc, export = FALSE)
#' PlotMBorderPCA(sc)
#' }

PlotMBorderPCA<- function(object) {
	data=object@MBclusters
	MBordertable<-cbind(data$pcareduceres,object@MBordering)
    l <- MBordertable[,3]
    mi <- min(l,na.rm=TRUE)
    ma <- max(l,na.rm=TRUE)
    ColorRamp <- colorRampPalette(rev(brewer.pal(n = 11,name = "RdYlBu")))(100)
    ColorLevels <- seq(mi, ma, length=length(ColorRamp))
    v <- round((l - mi)/(ma - mi)*99 + 1,0)
    layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
    par(mar = c(5,5,2.5,2))
    plot(MBordertable[,1],MBordertable[,2],xlab="PC1",ylab="PC2",pch=20,cex=0,col="grey",las=1)
    for ( k in 1:length(v) ){
        points(MBordertable[k,1],MBordertable[k,2],col=ColorRamp[v[k]],pch=20,cex=2)
    }
    par(mar = c(3,2.5,2.5,2))
    image(1, ColorLevels,matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),col=ColorRamp,xlab="",ylab="",las=2,xaxt="n")
    layout(1)
}