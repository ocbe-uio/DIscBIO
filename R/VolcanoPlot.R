#' @title Volcano Plot
#' @description Plotting differentially expressed genes (DEGs) in a particular cluster. 
#' Volcano plots are used to readily show the DEGs by plotting significance versus fold-change on the y and x axes, respectively.
#' @param object A data frame showing the differentially expressed genes (DEGs) in a particular cluster
#' @param Value A numeric value of the false discovery rate. Default is 0.05.. Default is 0.05
#' @param fc A numeric value of the fold change. Default is 0.5.
#' @param FS A numeric value of the font size. Default is 0.4.
#' @param name A string vector showing the name to be used to save the resulted tables.
#' @param adj A logical vector that allows adjusting the y axis by adding a very smal value (0.00000000001) to logFDR. Default is TRUE. 
#' @param export A logical vector that allows writing the final gene list in excel file. Default is TRUE. 
#' @importFrom samr samr samr.compute.delta.table samr.plot samr.compute.siggenes.table
#' @importFrom graphics title
#' @importFrom utils write.csv
#' @export
#' @rdname DEGanalysis

    if (length(object[1,])>8) {object<-object[,-1]}
    object[,8] <- if ( adj ) object[,8]+0.00000000001 else object[,8]
    with(object, plot(object[,7], -log10(object[,8]), pch=20,cex=2, las=1,xlab="log2 Fold Change",ylab="-log10 FDR",sub=paste0("Volcano plot ",DEGs),font.sub=4,col.sub="black"))
    FC<-subset(object, abs(object[,7])>fc)    # Fold Change
    sigFC<-subset(object, object[,8]<Value & abs(object[,7])>fc) # Significant genes
    with(FC, points(FC[,7], -log10(FC[,8]), pch=20,cex=2, col="red"))
    with(sigFC, points(sigFC[,7], -log10(sigFC[,8]), pch=20,cex=2, col="blue"))
    with(sigFC, textxy(sigFC[,7], -log10(sigFC[,8]), labs=sigFC[,2], cex=FS,col="blue"))
    add_legend("topleft", legend=c(paste0("DEGs (FC < ",fc," - FDR> ",Value,")   "), paste0("DEGs (FC > ",fc," - FDR> ",Value,")"),paste0("DEGs (FC > ",fc," - FDR< ",Value,")   ")), pch=20, col=c("black", "red","blue"),horiz=TRUE, bty='n', cex=0.7)
}