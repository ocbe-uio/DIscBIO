VolcanoPlot<-function(object,Value,DEGs,fc,adj=FALSE,FS=.4){
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