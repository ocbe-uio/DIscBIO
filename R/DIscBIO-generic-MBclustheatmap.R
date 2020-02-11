#' @title Plotting the Model-based clusters in a heatmap representation of the cell-to-cell distances
#' @description  This functions plots a heatmap of the distance matrix grouped by clusters. Individual clusters are highlighted with rainbow colors along the x and y-axes.
#' @param object \code{DISCBIO} class object.
#' @param hmethod  Agglomeration method used for determining the cluster order from hierarchical clustering of the cluster medoids. 
#' This should be one of "ward.D", "ward.D2", "single", "complete", "average". Default is "single".
#' @param plot if `TRUE`, plots the heatmap; otherwise, just prints cclmo
#' @param quiet if `TRUE`, intermediary output is suppressed
#' @importFrom stats hclust as.dist cor kmeans
#' @importFrom cluster clusGap maxSE
#' @importFrom fpc clusterboot kmeansCBI
#' @return Unless otherwise specified, a heatmap and a vector of the underlying cluster order.
#' @examples
#' sc <- DISCBIO(valuesG1msReduced)
#' sc <- NoiseFiltering(sc, export=FALSE)
#' sc <- Normalizedata(
#'     sc, mintotal=1000, minexpr=0, minnumber=0, maxexpr=Inf, downsample=FALSE,
#'     dsn=1, rseed=17000
#' )
#' sc <- FinalPreprocessing(sc, GeneFlitering="NoiseF", export=FALSE)
#' sc <- Exprmclust(sc,K = 2)
#' sc <- comptsneMB(sc, rseed=15555)
#' sc <- Clustexp(sc, cln=3)
#' sc <- MB_Order(sc, export = FALSE)
#' MBclustheatmap(sc, hmethod="single")
setGeneric("MBclustheatmap", function(object,hmethod="single", plot = TRUE, quiet = FALSE) standardGeneric("MBclustheatmap"))

#' @export
#' @rdname MBclustheatmap
setMethod("MBclustheatmap",
          signature = "DISCBIO",
          definition = function(object,hmethod, plot = TRUE, quiet = FALSE){
            x <- object@fdata  
    object@clusterpar$metric <- "pearson"
    dist.gen <- function(x,method="euclidean", ...) if ( method %in% c("spearman","pearson","kendall") ) as.dist( 1 - cor(t(x),method=method,...) ) else dist(x,method=method,...)
    dist.gen.pairs <- function(x,y,...) dist.gen(t(cbind(x,y)),...)
    clustfun <- function(x,clustnr=20,bootnr=50,metric="pearson",do.gap=TRUE,SE.method="Tibs2001SEmax",SE.factor=.25,B.gap=50,cln=0,rseed=17000,quiet = FALSE){
      if ( clustnr < 2) stop("Choose clustnr > 1")
      di <- dist.gen(t(x),method=metric)
      if ( do.gap | cln > 0 ){
        gpr <- NULL
        if ( do.gap ){
          set.seed(rseed)
          gpr <- clusGap(
            as.matrix(di), FUNcluster = kmeans, K.max = clustnr, B = B.gap,
            verbose = !quiet
          ) 
          if ( cln == 0 ) cln <- maxSE(gpr$Tab[,3],gpr$Tab[,4],method=SE.method,SE.factor)
        }    
        if ( cln <= 1 ) {
          clb <- list(result=list(partition=rep(1,dim(x)[2])),bootmean=1)
          names(clb$result$partition) <- names(x)
          return(list(x=x,clb=clb,gpr=gpr,di=di))
        }
        # FUN <- match.fun(clustermethod)
        clb <- clusterboot(
          di, B=bootnr, distances=FALSE, bootmethod="boot",
          clustermethod=fpc::kmeansCBI, krange=cln, scaling=FALSE, 
          multipleboot=FALSE, bscompare=TRUE, seed=rseed, count = !quiet
        )
        return(list(x=x,clb=clb,gpr=gpr,di=di))
      }
    }
     y <- clustfun(object@fdata,clustnr=20,bootnr=50,metric="pearson",do.gap=TRUE,SE.method="Tibs2001SEmax",SE.factor=.25,B.gap=50,cln=0,rseed=17000, quiet = quiet) 
    object@distances <- as.matrix( y$di )
            part <- object@MBclusters$clusterid
            na <- c()
            j <- 0
            for ( i in 1:max(part) ){
              if ( sum(part == i) == 0 ) next
              j <- j + 1
              na <- append(na,i)
              d <- x[,part == i]
              if ( sum(part == i) == 1 ) cent <- d else cent <- apply(d,1,mean)
              if ( j == 1 ) tmp <- data.frame(cent) else tmp <- cbind(tmp,cent)
            }
            names(tmp) <- paste("cl",na,sep=".")
            if ( max(part) > 1 )  cclmo <- hclust(dist.gen(as.matrix(dist.gen(t(tmp),method=object@clusterpar$metric))),method=hmethod)$order else cclmo <- 1
            q <- part
            for ( i in 1:max(part) ){
              q[part == na[cclmo[i]]] <- i
            }
            part <- q
            di <- as.data.frame( as.matrix( dist.gen(t(object@distances)) ) )
            pto <- part[order(part,decreasing=FALSE)]
            ptn <- c()
            for ( i in 1:max(pto) ){ pt <- names(pto)[pto == i]; z <- if ( length(pt) == 1 ) pt else pt[hclust(as.dist(t(di[pt,pt])),method=hmethod)$order]; ptn <- append(ptn,z) }
            col=c("black","blue","green","red","yellow","gray")
            mi  <- min(di,na.rm=TRUE)
            ma  <- max(di,na.rm=TRUE)

            if (plot) {
              layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
              ColorRamp   <- colorRampPalette(brewer.pal(n = 7,name = "RdYlBu"))(100)
              ColorLevels <- seq(mi, ma, length=length(ColorRamp))
              if ( mi == ma ){
                ColorLevels <- seq(0.99*mi, 1.01*ma, length=length(ColorRamp))
              }
              par(mar = c(3,5,2.5,2))
              image(as.matrix(di[ptn,ptn]),col=ColorRamp,axes=FALSE)
              abline(0,1)
              box()
              
              tmp <- c()
              for ( u in 1:max(part) ){
                ol <- (0:(length(part) - 1)/(length(part) - 1))[ptn %in% names(x)[part == u]]
                points(rep(0,length(ol)),ol,col=col[cclmo[u]],pch=15,cex=.75)
                points(ol,rep(0,length(ol)),col=col[cclmo[u]],pch=15,cex=.75)
                tmp <- append(tmp,mean(ol))
              }
              axis(1,at=tmp,labels=cclmo)
              axis(2,at=tmp,labels=cclmo)
              par(mar = c(3,2.5,2.5,2))
              image(1, ColorLevels,
                    matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
                    col=ColorRamp,
                    xlab="",ylab="",las=2,
                    xaxt="n")
              layout(1)
            }
            return(cclmo)
          }
          )