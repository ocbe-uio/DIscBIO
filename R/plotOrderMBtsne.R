#' @title Plotting the pseudo-time ordering based on Model-based clusters in the t-SNE map
#' @description The tSNE representation can also be used to show the pseudo-time ordering.
#' @param object \code{PSCANseq} class object.
#' of the \code{PSCANseq} object.
#' @export
#' @rdname plotOrderMBtsne
setGeneric("plotOrderMBtsne", function(object) standardGeneric("plotOrderMBtsne"))

setMethod("plotOrderMBtsne",
          signature = "PSCANseq",
          definition = function(object){
            if ( length(object@MBtsne) == 0 ) stop("run comptsneMB before plotOrderMBtsne")
			total<-rbind(object@ndata,object@MBordering)
			rownames(total)[nrow(total)]<-"Pseudo-time ordering of MBclustering" 
			g<-rownames(total)[nrow(total)]
            n <- g[1]
            l <- apply(total[g,] - .1,2,sum) + .1
            
            mi <- min(l,na.rm=TRUE)
            ma <- max(l,na.rm=TRUE)
            ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7,name = "RdYlBu")))(100)
            ColorLevels <- seq(mi, ma, length=length(ColorRamp))
            v <- round((l - mi)/(ma - mi)*99 + 1,0)
            layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
            par(mar = c(3,5,2.5,2))
            plot(object@MBtsne,xlab="Dim 1",ylab="Dim 2",main=n,pch=20,cex=0,col="grey",las=1)
            for ( k in 1:length(v) ){
              points(object@MBtsne[k,1],object@MBtsne[k,2],col=ColorRamp[v[k]],pch=20,cex=1.5)
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