#' PSCANseq
#' @slot expdata    
#' @slot ndata      
#' @slot fdata      
#' @slot distances  
#' @slot tsne       
#' @slot background 
#' @slot out        
#' @slot cpart      
#' @slot fcol       
#' @slot filterpar  
#' @slot clusterpar 
#' @slot outlierpar 
#' @slot kmeans     
#' @slot MBclusters 
#' @slot kordering  
#' @slot MBordering 
#' @slot MBtsne     
#' @importFrom methods new
#' @name PSCANseq
#' @rdname PSCANseq
#' @aliases PSCANseq-class, PSCANseq-class
#' @exportClass PSCANseq
PSCANseq <- setClass(
    Class = "PSCANseq",
    slots = c(
        expdata    = "data.frame",
        ndata      = "data.frame",
        fdata      = "data.frame", 
        distances  = "matrix",
        tsne       = "data.frame",
        background = "list",
        out        = "list", 
        cpart      = "vector",
        fcol       = "vector",
        filterpar  = "list",
        clusterpar = "list", 
        outlierpar = "list",
        kmeans     = "list",
        MBclusters = "vector",
        kordering  = "vector",
        MBordering = "vector",
        MBtsne     = "data.frame"
    )
)

setValidity("PSCANseq",
            function(object) {
              msg <- NULL
              if ( ! is.data.frame(object@expdata) ){
                msg <- c(msg, "input data must be data.frame")
              }else if ( nrow(object@expdata) < 2 ){
                msg <- c(msg, "input data must have more than one row")
              }else if ( ncol(object@expdata) < 2 ){
                msg <- c(msg, "input data must have more than one column")
              }else if (sum( apply( is.na(object@expdata),1,sum ) ) > 0 ){
                msg <- c(msg, "NAs are not allowed in input data")
              }else if (sum( apply( object@expdata,1,min ) ) < 0 ){
                msg <- c(msg, "negative values are not allowed in input data")
              }
              if (is.null(msg)) TRUE
              else msg
            }
            )

#' @title title
#' @description description
#' @param .Object .Object
#' @param expdata expdata
#' @importFrom methods validObject
#' @rdname initialize
setMethod("initialize",
          signature = "PSCANseq",
          definition = function(.Object, expdata ){
            .Object@expdata <- expdata
            .Object@ndata <- expdata
            .Object@fdata <- expdata
            validObject(.Object)
            return(.Object)
          }
          )

#' @title title
#' @export
#' @rdname Normalizedata
#' @param object object
#' @param mintotal mintotal
#' @param minexpr minexpr
#' @param minnumber minnumber
#' @param maxexpr maxexpr
#' @param downsample downsample
#' @param dsn dsn
#' @param rseed rseed
setGeneric("Normalizedata", function(object, mintotal=1000, minexpr=0, minnumber=0, maxexpr=Inf, downsample=FALSE, dsn=1, rseed=17000) standardGeneric("Normalizedata"))

#' @rdname Normalizedata
setMethod("Normalizedata",
          signature = "PSCANseq",
          definition = function(object,mintotal,minexpr,minnumber,maxexpr,downsample,dsn,rseed) {
            if ( ! is.numeric(mintotal) ) stop( "mintotal has to be a positive number" ) else if ( mintotal <= 0 ) stop( "mintotal has to be a positive number" )
            if ( ! is.numeric(minexpr) ) stop( "minexpr has to be a non-negative number" ) else if ( minexpr < 0 ) stop( "minexpr has to be a non-negative number" )
            if ( ! is.numeric(minnumber) ) stop( "minnumber has to be a non-negative integer number" ) else if ( round(minnumber) != minnumber | minnumber < 0 ) stop( "minnumber has to be a non-negative integer number" )
            if ( ! ( is.numeric(downsample) | is.logical(downsample) ) ) stop( "downsample has to be logical (TRUE/FALSE)" )
            if ( ! is.numeric(dsn) ) stop( "dsn has to be a positive integer number" ) else if ( round(dsn) != dsn | dsn <= 0 ) stop( "dsn has to be a positive integer number" )
            object@filterpar <- list(mintotal=mintotal, minexpr=minexpr, minnumber=minnumber, maxexpr=maxexpr, downsample=downsample, dsn=dsn)
            object@ndata <- object@expdata[,apply(object@expdata,2,sum,na.rm=TRUE) >= mintotal]
            if ( downsample ){
              set.seed(rseed)
              object@ndata <- downsample(object@expdata,n=mintotal,dsn=dsn)
            }else{
              x <- object@ndata
              object@ndata <- as.data.frame( t(t(x)/apply(x,2,sum))*median(apply(x,2,sum,na.rm=TRUE)) + .1 )
            }
            x <- object@ndata
            object@fdata <- x[apply(x>=minexpr,1,sum,na.rm=TRUE) >= minnumber,]
            x <- object@fdata
            object@fdata <- x[apply(x,1,max,na.rm=TRUE) < maxexpr,]
            return(object)
          }
          )

#' @title title
#' @export
#' @rdname Clustexp
#' @docType methods
#' @param object object
#' @param clustnr clustnr
#' @param bootnr bootnr
#' @param metric metric
#' @param do.gap do.gap
#' @param SE.method SE.method
#' @param SE.factor SE.factor
#' @param B.gap B.gap
#' @param cln cln
#' @param rseed rseed
setGeneric("Clustexp", function(object, clustnr = 20, bootnr = 50,
                                metric = "pearson", do.gap = TRUE,
                                SE.method = "Tibs2001SEmax", SE.factor = .25,
                                B.gap = 50, cln = 0, rseed = 17000) {
        standardGeneric("Clustexp")
})

#' @rdname Clustexp
setMethod(
    f = "Clustexp",
    signature = "PSCANseq",
    definition = function(object, clustnr, bootnr, metric, do.gap, SE.method,
                          SE.factor, B.gap, cln, rseed) {
        # Validation
        if ( ! is.numeric(clustnr) ) stop("clustnr has to be a positive integer") else if ( round(clustnr) != clustnr | clustnr <= 0 ) stop("clustnr has to be a positive integer")
        if ( ! is.numeric(bootnr) ) stop("bootnr has to be a positive integer") else if ( round(bootnr) != bootnr | bootnr <= 0 ) stop("bootnr has to be a positive integer")
        if ( ! ( metric %in% c( "spearman","pearson","kendall","euclidean","maximum","manhattan","canberra","binary","minkowski") ) ) stop("metric has to be one of the following: spearman, pearson, kendall, euclidean, maximum, manhattan, canberra, binary, minkowski")
        if ( ! ( SE.method %in% c( "firstSEmax","Tibs2001SEmax","globalSEmax","firstmax","globalmax") ) ) stop("SE.method has to be one of the following: firstSEmax, Tibs2001SEmax, globalSEmax, firstmax, globalmax")
        if ( ! is.numeric(SE.factor) ) stop("SE.factor has to be a non-negative integer") else if  ( SE.factor < 0 )  stop("SE.factor has to be a non-negative integer")
        if ( ! ( is.numeric(do.gap) | is.logical(do.gap) ) ) stop( "do.gap has to be logical (TRUE/FALSE)" )
        if ( ! is.numeric(B.gap) ) stop("B.gap has to be a positive integer") else if ( round(B.gap) != B.gap | B.gap <= 0 ) stop("B.gap has to be a positive integer")
        if ( ! is.numeric(cln) ) stop("cln has to be a non-negative integer") else if ( round(cln) != cln | cln < 0 ) stop("cln has to be a non-negative integer")          
        if ( ! is.numeric(rseed) ) stop("rseed has to be numeric")
        if ( !do.gap & cln == 0 ) stop("cln has to be a positive integer or do.gap has to be TRUE")
        
        # Operations
        object@clusterpar <- list(clustnr=clustnr,bootnr=bootnr,metric=metric,do.gap=do.gap,SE.method=SE.method,SE.factor=SE.factor,B.gap=B.gap,cln=cln,rseed=rseed)
        y <- clustfun(object@fdata,clustnr,bootnr,metric,do.gap,SE.method,SE.factor,B.gap,cln,rseed)
        object@kmeans    <- list(kpart=y$clb$result$partition, jaccard=y$clb$bootmean, gap=y$gpr)
        object@distances <- as.matrix( y$di )
        set.seed(111111)
        object@fcol <- sample(rainbow(max(y$clb$result$partition)))
        return(object)
})

#' @title title
#' @export
#' @rdname plotGap
#' @param object object
setGeneric("plotGap", function(object) standardGeneric("plotGap"))

#' @rdname plotGap
setMethod("plotGap",
          signature = "PSCANseq",
          definition = function(object){
            if ( length(object@kmeans$kpart) == 0 ) stop("run clustexp before plotgap")
            plot(object@kmeans$gap,ylim=c(0.1,0.5),las=1,main="Gap Statistics")
          }
          )

#' @title title
#' @export
#' @rdname comptSNE
setGeneric("comptSNE", function(object,rseed=15555) standardGeneric("comptSNE"))

#' @title title
#' @description description
#' @param object object
#' @param rseed rseed
#' @importFrom tsne tsne
#' @rdname comptSNE
#' @export
setMethod("comptSNE",
          signature = "PSCANseq",
          definition = function(object,rseed){
            if ( length(object@kmeans$kpart) == 0 ) stop("run clustexp before comptsne")
            set.seed(rseed)
            di <- dist.gen(as.matrix(object@distances))
            ts <- tsne(di,k=2)
            object@tsne <- as.data.frame(ts)
            return(object)
          }
          )






            

#' @title title
#' @export
#' @rdname plotSilhouette
#' @param K K
setGeneric("plotSilhouette", function(object,K) standardGeneric("plotSilhouette"))
#' @title title
#' @description description
#' @param object object
#' @importFrom cluster silhouette
#' @rdname plotSilhouette
#' @export
setMethod("plotSilhouette",
          signature = "PSCANseq",
          definition = function(object){
            if ( length(object@kmeans$kpart) == 0 ) stop("run clustexp before plotsilhouette")
            if ( length(unique(object@kmeans$kpart)) < 2 ) stop("only a single cluster: no silhouette plot")
            col=c("black","blue","green","red","yellow","gray")
		kpart <- object@kmeans$kpart
            distances  <- dist.gen(object@distances)
            si <- silhouette(kpart,distances)
            plot(si,col=col[1:K])
          }
          )

#' @export
#' @title title
#' @rdname plottSNE
#' @param object object
setGeneric("plottSNE", function(object) standardGeneric("plottSNE"))
#' @title title
#' @description description
#' @importFrom graphics text
#' @rdname plotSilhouette
#' @export
setMethod("plottSNE",
          signature = "PSCANseq",
          definition = function(object){
            if ( length(object@tsne) == 0 ) stop("run comptsne before plottsne")
		col=c("black","blue","green","red","yellow","gray")
            part <- object@kmeans$kpart
            plot(object@tsne,las=1,xlab="Dim 1",ylab="Dim 2",pch=20,cex=1.5,col="lightgrey")
            for ( i in 1:max(part) ){
              if ( sum(part == i) > 0 ) text(object@tsne[part == i,1],object@tsne[part == i,2],i,col=col[i],cex=.75,font=4)
            }
          }
          )


#' @title title
#' @rdname plotKmeansLabelstSNE
setGeneric("plotKmeansLabelstSNE", function(object) standardGeneric("plotKmeansLabelstSNE"))
#' @title title
#' @description description
#' @param object object
#' @rdname plotKmeansLabelstSNE
#' @importFrom graphics text
setMethod("plotKmeansLabelstSNE",
          signature = "PSCANseq",
          definition = function(object){
            if ( length(object@tsne) == 0 ) stop("run comptsne before plotKmeansLabelstSNE")
            Clusters<-object@kmeans$kpart
            ClustersFactor<- as.factor(Clusters)
            ClustersFactor<- gsub("1", "black", ClustersFactor)
            ClustersFactor<- gsub("2", "blue", ClustersFactor)
            ClustersFactor<- gsub("3", "green", ClustersFactor)
            ClustersFactor<- gsub("4", "red", ClustersFactor)
            ClustersFactor<- gsub("5", "yellow", ClustersFactor)
            ClustersFactor<- gsub("6", "gray", ClustersFactor)
            COL<-ClustersFactor
            labels=names(object@ndata)
            plot(object@tsne,xlab="Dim 1",ylab="Dim 2",pch=20,cex=.5,col="lightgrey")
            text(object@tsne[,1],object@tsne[,2],labels,cex=.7,col=COL)
          }
          )  

#' @title title
#' @export
#' @rdname plotSymbolstSNE
#' @param object object
#' @param types types
setGeneric("plotSymbolstSNE", function(object,types=NULL) standardGeneric("plotSymbolstSNE"))

#' @export
#' @rdname plotSymbolstSNE
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

#' @title title
#' @export
#' @rdname KMclustheatmap
setGeneric("KMclustheatmap", function(object,hmethod="single") standardGeneric("KMclustheatmap"))

#' @title title
#' @description description
#' @param object object
#' @param hmethod hmethod
#' @importFrom stats hclust
#' @export
#' @rdname KMclustheatmap
setMethod("KMclustheatmap",
          signature = "PSCANseq",
          definition = function(object,hmethod){
            x <- object@fdata  
            part <- object@kmeans$kpart
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
            return(cclmo)
          }
          )

#' @title title
#' @export
#' @rdname MBclustheatmap
#' @param object object
#' @param hmethod hmethod
setGeneric("MBclustheatmap", function(object,hmethod="single") standardGeneric("MBclustheatmap"))

#' @export
#' @rdname MBclustheatmap
setMethod("MBclustheatmap",
          signature = "PSCANseq",
          definition = function(object,hmethod){
            x <- object@fdata  
		object@clusterpar$metric <- "pearson"
 		y <- clustfun(object@fdata,clustnr=20,bootnr=50,metric="pearson",do.gap=TRUE,SE.method="Tibs2001SEmax",SE.factor=.25,B.gap=50,cln=0,rseed=17000) 
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
            return(cclmo)
          }
          )

#' @title title
#' @export
#' @rdname plotExptSNE
#' @param object object
#' @param g g
#' @param n n
setGeneric("plotExptSNE", function(object,g,n="") standardGeneric("plotExptSNE"))

#' @export
#' @rdname plotExptSNE
setMethod("plotExptSNE",
          signature = "PSCANseq",
          definition = function(object,g,n=""){
            if ( length(object@tsne) == 0 ) stop("run comptsne before plotExptSNE")
            if ( length(intersect(g,rownames(object@ndata))) < length(unique(g)) ) stop("second argument does not correspond to set of rownames slot ndata of SCseq object")
            if ( n == "" ) n <- g[1]
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

#' @title title
#' @export
#' @rdname comptsneMB
setGeneric("comptsneMB", function(object,rseed=15555) standardGeneric("comptsneMB"))
#' @title title
#' @description description
#' @param object object
#' @param rseed rseed
#' @importFrom tsne tsne
#' @rdname comptsneMB
#' @export
setMethod("comptsneMB",
          signature = "PSCANseq",
          definition = function(object,rseed){
            if ( length(object@MBclusters) == 0 ) stop("run clustexp before comptsneMB")
            set.seed(rseed)
            dist.gen <- function(x,method="euclidean", ...) if ( method %in% c("spearman","pearson","kendall") ) as.dist( 1 - cor(t(x),method=method,...) ) else dist(x,method=method,...)
            di <- dist.gen(as.matrix(t(procdataTSCAN)))
            cat("This function takes time")
            ts <- tsne(di,k=2,max_iter = 5000,epoch=500)
            object@MBtsne <- as.data.frame(ts)
            return(object)
          }
          )
           
#' @title title
#' @export
#' @rdname plottsneMB
#' @param K K
setGeneric("plottsneMB", function(object,K) standardGeneric("plottsneMB"))
#' @title title
#' @description description
#' @param object object
#' @importFrom graphics text
#' @export
#' @rdname plottsneMB
setMethod("plottsneMB",
          signature = "PSCANseq",
          definition = function(object){
            if ( length(object@MBtsne) == 0 ) stop("run comptsneMB before plottsneMB")
		col=c("black","blue","green","red","yellow","gray")
            part <- object@MBclusters$clusterid
            plot(object@MBtsne,xlab="Dim 1",ylab="Dim 2",pch=20,cex=1.5,col="lightgrey",las=1)
            for ( i in 1:K ){
              if ( sum(part == i) > 0 ) text(object@MBtsne[part == i,1],object@MBtsne[part == i,2],i,col=col[i],cex=.75,font=4)
            }
          }
          )    

#' @title title
#' @export
#' @rdname plotMBLabelstSNE
setGeneric("plotMBLabelstSNE", function(object) standardGeneric("plotMBLabelstSNE"))
#' @title title
#' @description description
#' @param object object
#' @importFrom graphics text
#' @rdname plotMBLabelstSNE
#' @export
setMethod("plotMBLabelstSNE",
          signature = "PSCANseq",
          definition = function(object){
            if ( length(object@MBtsne) == 0 ) stop("run comptsneMB before plotMBLabelstSNE")
            Clusters<-object@MBclusters$clusterid
            ClustersFactor<- as.factor(Clusters)
            ClustersFactor<- gsub("1", "black", ClustersFactor)
            ClustersFactor<- gsub("2", "blue", ClustersFactor)
            ClustersFactor<- gsub("3", "green", ClustersFactor)
            ClustersFactor<- gsub("4", "red", ClustersFactor)
            ClustersFactor<- gsub("5", "yellow", ClustersFactor)
            ClustersFactor<- gsub("6", "gray", ClustersFactor)
            COL<-ClustersFactor
            labels=names(object@ndata)
            plot(object@MBtsne,xlab="Dim 1",ylab="Dim 2",pch=20,cex=.5,col="lightgrey")
            text(object@MBtsne[,1],object@MBtsne[,2],labels,cex=.7,col=COL)
          }
          ) 

#' @title title
#' @export
#' @rdname plotsilhouetteMB
#' @param K k
setGeneric("plotsilhouetteMB", function(object,K) standardGeneric("plotsilhouetteMB"))
#' @title title
#' @description description
#' @param object object
#' @importFrom cluster silhouette
#' @rdname plotsilhouetteMB
#' @export
setMethod("plotsilhouetteMB",
          signature = "PSCANseq",
          definition = function(object){
            if ( length(object@MBclusters$clusterid) == 0 ) stop("run exprmclust before plotsilhouetteMB")
            if ( length(unique(object@MBclusters$clusterid)) < 2 ) stop("only a single cluster: no silhouette plot")
		col=c("black","blue","green","red","yellow","gray")
            kpart <- object@MBclusters$clusterid
            distances  <- dist.gen(t(object@fdata))
            si <- silhouette(kpart,distances)
            plot(si,col=col[1:K])
          }
          )

#' @title title
#' @export
#' @rdname plotexptsneMB
#' @param object object
#' @param g g
#' @param n n
setGeneric("plotexptsneMB", function(object,g,n="") standardGeneric("plotexptsneMB"))

#' @export
#' @rdname plotexptsneMB
setMethod("plotexptsneMB",
          signature = "PSCANseq",
          definition = function(object,g,n=""){
            if ( length(object@MBtsne) == 0 ) stop("run comptsneMB before plotexptsneMB")
            if ( length(intersect(g,rownames(object@ndata))) < length(unique(g)) ) stop("second argument does not correspond to set of rownames slot ndata of SCseq object")
            if ( n == "" ) n <- g[1]
            l <- apply(object@ndata[g,] - .1,2,sum) + .1
            
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

#' @title title
#' @export
#' @rdname FindOutliersKM
setGeneric("FindOutliersKM", function(object,outminc=5,outlg=2,probthr=1e-3,thr=2**-(1:40),outdistquant=.75) standardGeneric("FindOutliersKM"))

#' @title title
#' @description description
#' @param object object
#' @param outminc outminc
#' @param outlg outlg
#' @param probthr probthr
#' @param thr thr
#' @param outdistquant outdistquant
#' @importFrom stats coef pnbinom
#' @importFrom amap K
#' @rdname FindOutliersKM
#' @export
setMethod("FindOutliersKM",
          signature = "PSCANseq",
          definition = function(object,outminc,outlg,probthr,thr,outdistquant) {
            if ( length(object@kmeans$kpart) == 0 ) stop("run clustexp before FindOutliersKM")
            if ( ! is.numeric(outminc) ) stop("outminc has to be a non-negative integer") else if ( round(outminc) != outminc | outminc < 0 ) stop("outminc has to be a non-negative integer")
            if ( ! is.numeric(outlg) ) stop("outlg has to be a non-negative integer") else if ( round(outlg) != outlg | outlg < 0 ) stop("outlg has to be a non-negative integer")
            if ( ! is.numeric(probthr) ) stop("probthr has to be a number between 0 and 1") else if (  probthr < 0 | probthr > 1 ) stop("probthr has to be a number between 0 and 1")
            if ( ! is.numeric(thr) ) stop("thr hast to be a vector of numbers between 0 and 1") else if ( min(thr) < 0 | max(thr) > 1 ) stop("thr hast to be a vector of numbers between 0 and 1")
            if ( ! is.numeric(outdistquant) ) stop("outdistquant has to be a number between 0 and 1") else if (  outdistquant < 0 | outdistquant > 1 ) stop("outdistquant has to be a number between 0 and 1")
                      
            object@outlierpar <- list( outminc=outminc,outlg=outlg,probthr=probthr,thr=thr,outdistquant=outdistquant )
            ### calibrate background model
            m <- log2(apply(object@fdata,1,mean))
            v <- log2(apply(object@fdata,1,var))
            f <- m > -Inf & v > -Inf
            m <- m[f]
            v <- v[f]
            mm <- -8
            repeat{
              fit <- lm(v ~ m + I(m^2)) 
              if( coef(fit)[3] >= 0 | mm >= 3){
                break
              }
              mm <- mm + .5
              f <- m > mm
              m <- m[f]
              v <- v[f]
            }
            object@background <- list()
            object@background$vfit <- fit
            object@background$lvar <- function(x,object) 2**(coef(object@background$vfit)[1] + log2(x)*coef(object@background$vfit)[2] + coef(object@background$vfit)[3] * log2(x)**2)
            object@background$lsize <- function(x,object) x**2/(max(x + 1e-6,object@background$lvar(x,object)) - x)

            ### identify outliers
            out   <- c()
            stest <- rep(0,length(thr))
            cprobs <- c()
            for ( n in 1:max(object@kmeans$kpart) ){     
              if ( sum(object@kmeans$kpart == n) == 1 ){ cprobs <- append(cprobs,.5); names(cprobs)[length(cprobs)] <- names(object@kmeans$kpart)[object@kmeans$kpart == n]; next }
              x <- object@fdata[,object@kmeans$kpart == n]
              x <- x[apply(x,1,max) > outminc,]
              z <- t( apply(x,1,function(x){ apply( cbind( pnbinom(round(x,0),mu=mean(x),size=object@background$lsize(mean(x),object)) , 1 - pnbinom(round(x,0),mu=mean(x),size=object@background$lsize(mean(x),object)) ),1, min) } ) )
              cp <- apply(z,2,function(x){ y <- p.adjust(x,method="BH"); y <- y[order(y,decreasing=FALSE)]; return(y[outlg]);})
              f <- cp < probthr
              cprobs <- append(cprobs,cp)
              if ( sum(f) > 0 ) out <- append(out,names(x)[f])
              for ( j in 1:length(thr) )  stest[j] <-  stest[j] + sum( cp < thr[j] )  
            }
            object@out <-list(out=out,stest=stest,thr=thr,cprobs=cprobs)

            ### cluster outliers
            clp2p.cl <- c()
            cols <- names(object@fdata)
            di <- as.data.frame(object@distances)
            for ( i in 1:max(object@kmeans$kpart) ) {
              tcol <- cols[object@kmeans$kpart == i]
              if ( sum(!(tcol %in% out)) > 1 ) clp2p.cl <- append(clp2p.cl,as.vector(t(di[tcol[!(tcol %in% out)],tcol[!(tcol %in% out)]])))
            }
            clp2p.cl <- clp2p.cl[clp2p.cl>0]
  
            cpart <- object@kmeans$kpart
            cadd  <- list()
            if ( length(out) > 0 ){
              if (length(out) == 1){
                cadd <- list(out)
              }else{
                n <- out
                m <- as.data.frame(di[out,out])
                
                for ( i in 1:length(out) ){
                  if ( length(n) > 1 ){
                    o   <- order(apply(cbind(m,1:dim(m)[1]),1,function(x)  min(x[1:(length(x)-1)][-x[length(x)]])),decreasing=FALSE)
                    m <- m[o,o]
                    n <- n[o]          
                    f <- m[,1] < quantile(clp2p.cl,outdistquant) | m[,1] == min(clp2p.cl)
                    ind <- 1  
                    if ( sum(f) > 1 ) for ( j in 2:sum(f) ) if ( apply(m[f,f][j,c(ind,j)] > quantile(clp2p.cl,outdistquant) ,1,sum) == 0 ) ind <- append(ind,j)
                    cadd[[i]] <- n[f][ind]
                    g <- ! n %in% n[f][ind]
                    n <- n[g]
                    m <- m[g,g]
                    if ( sum(g) == 0 ) break
          
                  }else if (length(n) == 1){
                    cadd[[i]] <- n
                    break
                  }
                }
              }
    
              for ( i in 1:length(cadd) ){
                cpart[cols %in% cadd[[i]]] <- max(cpart) + 1
              }
            }

            ### determine final clusters
            
            object@cpart <- cpart

            set.seed(111111)
            object@fcol <- sample(rainbow(max(cpart)))
                           p <- object@kmeans$kpart[ order(object@kmeans$kpart,decreasing=FALSE)]
            x <- object@out$cprobs[names(p)]
            fcol <- c("black","blue","green","red","yellow","gray")

            for ( i in 1:max(p) ){
              y <- -log10(x + 2.2e-16)
              y[p != i] <- 0
              if ( i == 1 ) b <- barplot(y,ylim=c(0,max(-log10(x + 2.2e-16))*1.1),col=fcol[i],border=fcol[i],names.arg=FALSE,ylab="-log10prob") else barplot(y,add=TRUE,col=fcol[i],border=fcol[i],names.arg=FALSE,axes=FALSE)
  }
            abline(-log10(object@outlierpar$probthr),0,col="black",lty=2)
            d <- b[2,1] - b[1,1]
            y <- 0
            for ( i in 1:max(p) ) y <- append(y,b[sum(p <=i),1] + d/2)
            axis(1,at=(y[1:(length(y)-1)] + y[-1])/2,labels=1:max(p))
            box()
            cat("The following cells are considered as outlier cells:",which(object@cpart>K),"\n")
            print(which(object@cpart>K))
            LL= which(object@cpart>K)  
            return(LL)
          }
        )

#' @title title
#' @export
#' @rdname FindOutliersMB
setGeneric("FindOutliersMB", function(object,outminc=5,outlg=2,probthr=1e-3,thr=2**-(1:40),outdistquant=.75) standardGeneric("FindOutliersMB"))

#' @title title
#' @description description
#' @param object object
#' @param outminc outminc
#' @param outlg outlg
#' @param probthr probthr
#' @param thr thr
#' @param outdistquant outdistquant
#' @importFrom stats pnbinom
#' @rdname FindOutliersMB
#' @export
setMethod("FindOutliersMB",
          signature = "PSCANseq",
          definition = function(object,outminc,outlg,probthr,thr,outdistquant) {
            if ( length(object@MBclusters$clusterid) == 0 ) stop("run exprmclust before FindOutliersMB")
            if ( ! is.numeric(outminc) ) stop("outminc has to be a non-negative integer") else if ( round(outminc) != outminc | outminc < 0 ) stop("outminc has to be a non-negative integer")
            if ( ! is.numeric(outlg) ) stop("outlg has to be a non-negative integer") else if ( round(outlg) != outlg | outlg < 0 ) stop("outlg has to be a non-negative integer")
            if ( ! is.numeric(probthr) ) stop("probthr has to be a number between 0 and 1") else if (  probthr < 0 | probthr > 1 ) stop("probthr has to be a number between 0 and 1")
            if ( ! is.numeric(thr) ) stop("thr hast to be a vector of numbers between 0 and 1") else if ( min(thr) < 0 | max(thr) > 1 ) stop("thr hast to be a vector of numbers between 0 and 1")
            if ( ! is.numeric(outdistquant) ) stop("outdistquant has to be a number between 0 and 1") else if (  outdistquant < 0 | outdistquant > 1 ) stop("outdistquant has to be a number between 0 and 1")
            object<- Clustexp(object, clustnr=20,bootnr=50,metric="pearson",do.gap=T,SE.method="Tibs2001SEmax",SE.factor=.25,B.gap=50,cln=K,rseed=17000)
			object@outlierpar <- list( outminc=outminc,outlg=outlg,probthr=probthr,thr=thr,outdistquant=outdistquant )
            
			### calibrate background model
            m <- log2(apply(object@fdata,1,mean))
            v <- log2(apply(object@fdata,1,var))
            f <- m > -Inf & v > -Inf
            m <- m[f]
            v <- v[f]
            mm <- -8
            repeat{
              fit <- lm(v ~ m + I(m^2)) 
              if( coef(fit)[3] >= 0 | mm >= 3){
                break
              }
              mm <- mm + .5
              f <- m > mm
              m <- m[f]
              v <- v[f]
            }
            object@background <- list()
            object@background$vfit <- fit
            object@background$lvar <- function(x,object) 2**(coef(object@background$vfit)[1] + log2(x)*coef(object@background$vfit)[2] + coef(object@background$vfit)[3] * log2(x)**2)
            object@background$lsize <- function(x,object) x**2/(max(x + 1e-6,object@background$lvar(x,object)) - x)

            ### identify outliers
            out   <- c()
            stest <- rep(0,length(thr))
            cprobs <- c()
            for ( n in 1:max(object@MBclusters$clusterid) ){     
              if ( sum(object@MBclusters$clusterid == n) == 1 ){ cprobs <- append(cprobs,.5); names(cprobs)[length(cprobs)] <- names(object@MBclusters$clusterid)[object@MBclusters$clusterid == n]; next }
              x <- object@fdata[,object@MBclusters$clusterid == n]
              x <- x[apply(x,1,max) > outminc,]
              z <- t( apply(x,1,function(x){ apply( cbind( pnbinom(round(x,0),mu=mean(x),size=object@background$lsize(mean(x),object)) , 1 - pnbinom(round(x,0),mu=mean(x),size=object@background$lsize(mean(x),object)) ),1, min) } ) )
              cp <- apply(z,2,function(x){ y <- p.adjust(x,method="BH"); y <- y[order(y,decreasing=FALSE)]; return(y[outlg]);})
              f <- cp < probthr
              cprobs <- append(cprobs,cp)
              if ( sum(f) > 0 ) out <- append(out,names(x)[f])
              for ( j in 1:length(thr) )  stest[j] <-  stest[j] + sum( cp < thr[j] )  
            }
            object@out <-list(out=out,stest=stest,thr=thr,cprobs=cprobs)

            ### cluster outliers
            clp2p.cl <- c()
            cols <- names(object@fdata)
            di <- as.data.frame(object@distances)
            for ( i in 1:max(object@MBclusters$clusterid) ) {
              tcol <- cols[object@MBclusters$clusterid == i]
              if ( sum(!(tcol %in% out)) > 1 ) clp2p.cl <- append(clp2p.cl,as.vector(t(di[tcol[!(tcol %in% out)],tcol[!(tcol %in% out)]])))
            }
            clp2p.cl <- clp2p.cl[clp2p.cl>0]
  
            cpart <- object@MBclusters$clusterid
            cadd  <- list()
            if ( length(out) > 0 ){
              if (length(out) == 1){
                cadd <- list(out)
              }else{
                n <- out
                m <- as.data.frame(di[out,out])
                
                for ( i in 1:length(out) ){
                  if ( length(n) > 1 ){
                    o   <- order(apply(cbind(m,1:dim(m)[1]),1,function(x)  min(x[1:(length(x)-1)][-x[length(x)]])),decreasing=FALSE)
                    m <- m[o,o]
                    n <- n[o]          
                    f <- m[,1] < quantile(clp2p.cl,outdistquant) | m[,1] == min(clp2p.cl)
                    ind <- 1  
                    if ( sum(f) > 1 ) for ( j in 2:sum(f) ) if ( apply(m[f,f][j,c(ind,j)] > quantile(clp2p.cl,outdistquant) ,1,sum) == 0 ) ind <- append(ind,j)
                    cadd[[i]] <- n[f][ind]
                    g <- ! n %in% n[f][ind]
                    n <- n[g]
                    m <- m[g,g]
                    if ( sum(g) == 0 ) break
          
                  }else if (length(n) == 1){
                    cadd[[i]] <- n
                    break
                  }
                }
              }
    
              for ( i in 1:length(cadd) ){
                cpart[cols %in% cadd[[i]]] <- max(cpart) + 1
              }
            }

            ### determine final clusters
            
            object@cpart <- cpart

            set.seed(111111)
            object@fcol <- sample(rainbow(max(cpart)))
                           p <- object@MBclusters$clusterid[ order(object@MBclusters$clusterid,decreasing=FALSE)]
            x <- object@out$cprobs[names(p)]
            fcol <- c("black","blue","green","red","yellow","gray")

            for ( i in 1:max(p) ){
              y <- -log10(x + 2.2e-16)
              y[p != i] <- 0
              if ( i == 1 ) b <- barplot(y,ylim=c(0,max(-log10(x + 2.2e-16))*1.1),col=fcol[i],border=fcol[i],names.arg=FALSE,ylab="-log10prob") else barplot(y,add=TRUE,col=fcol[i],border=fcol[i],names.arg=FALSE,axes=FALSE)
  }
            abline(-log10(object@outlierpar$probthr),0,col="black",lty=2)
            d <- b[2,1] - b[1,1]
            y <- 0
            for ( i in 1:max(p) ) y <- append(y,b[sum(p <=i),1] + d/2)
            axis(1,at=(y[1:(length(y)-1)] + y[-1])/2,labels=1:max(p))
            box()
            cat("The following cells are considered as outlier cells:",which(object@cpart>K),"\n")
            print(which(object@cpart>K))
            LL= which(object@cpart>K)  
            return(LL)
          }
        )

#' @export
#' @title title
#' @rdname MBClustDiffGenes
setGeneric("MBClustDiffGenes", function(object,fdr=.01) standardGeneric("MBClustDiffGenes"))
#' @export
#' @rdname MBClustDiffGenes
#' @param object object
#' @param fdr fdr
setMethod("MBClustDiffGenes",
          signature = "PSCANseq",
          definition = function(object,fdr){
            if ( ! is.numeric(fdr) ) stop("pvalue has to be a number between 0 and 1") else if (  fdr < 0 | fdr > 1 ) stop("fdr has to be a number between 0 and 1")
            cdiff <- list()
            x     <- object@ndata
            y     <- object@expdata[,names(object@ndata)]
            part  <- object@MBclusters$clusterid
            for ( i in 1:max(part) ){
              if ( sum(part == i) == 0 ) next
              m <- apply(x,1,mean)
              n <- if ( sum(part == i) > 1 ) apply(x[,part == i],1,mean) else x[,part == i]
              no <- if ( sum(part == i) > 1 ) median(apply(y[,part == i],2,sum))/median(apply(x[,part == i],2,sum)) else sum(y[,part == i])/sum(x[,part == i])
              m <- m*no
              n <- n*no
              pv <- binompval(m/sum(m),sum(n),n)
              d <- data.frame(mean.all=m,mean.cl=n,fc=n/m,pv=pv)[order(pv,decreasing=FALSE),]
              #cdiff[[paste("cl",i,sep=".")]] <- d[d$pv < pvalue,]
              cdiff[[i]] <- d[d$pv < fdr,]
            }
            DEGsE<-c()
            DEGsS<-c()
            DEGsTable<-data.frame()

            for (n in 1:K){
                p.adj<- p.adjust(cdiff[[n]][,4],method="bonferroni")
                out<-cbind(cdiff[[n]], p.adj)
                out<-subset(out,out[,5]<fdr)
                Regulation<-c()
                for (i in 1:length(out[,1])){
                    if(out[i,1]>out[i,2]){
                        Regulation[i]="Down"
                    }else{
                        Regulation[i]="Up"
                    }
                }
                out<-cbind(out,Regulation)
                mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
                genes <- rownames(out)
                G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
                Final<-cbind(genes,out)
                Final<-merge(Final,G_list,by.x="genes",by.y="ensembl_gene_id")
                Final<-Final[!duplicated(Final[,8]), ]
                rownames(Final)<-Final[,1]
                Final[,1]<-Final[,8]
                Final<-Final[,-8]
    
                DEGsS<-c(DEGsS,Final[,1])
                DEGsE<-c(DEGsE,as.character(rownames(Final)))
                
                Up<-subset(Final,Final[,7]=="Up")
                Up<-select(Up, "Regulation","genes","pv","mean.all", "mean.cl","fc","p.adj")
                Up[,3]<-rownames(Up)
                Up[,6]<-log2(Up[,6])
                Up[,1]<-Up[,2]
                colnames(Up)<-c("Genes","genes","E.genes","mean.all", "mean.cl","log2.fc","p.adj")
                write.csv(Up, file = paste0("Up-DEG-cluster",n,".csv"))

                Down<-subset(Final,Final[,7]=="Down")
                Down<-select(Down, "Regulation","genes","pv","mean.all", "mean.cl","fc","p.adj")
                Down[,3]<-rownames(Down)
                Down[,6]<-log2(Down[,6])
                Down[,1]<- Down[,2]
                colnames(Down)<-c("Genes","genes","E.genes","mean.all", "mean.cl","log2.fc","p.adj")
                write.csv(Down, file = paste0("Down-DEG-cluster",n,".csv"))
    
    
                sigDEG<-cbind(DEGsE,DEGsS)
                write.csv(sigDEG, file = "binomial-sigDEG.csv")
                
                DEGsTable[n,1]<-paste0("Cluster ",n)
                DEGsTable[n,2]<-"Remaining Clusters"
                DEGsTable[n,3]<-length(Up[,1])
                DEGsTable[n,4]<-paste0("Up-DEG-cluster",n,".csv")
                DEGsTable[n,5]<-length(Down[,1])
                DEGsTable[n,6]<- paste0("Down-DEG-cluster",n,".csv")
            }
            colnames(DEGsTable)<-c("Target Cluster","VS","Up-regulated genes","File name","Low-regulated genes","File name")
            write.csv(DEGsTable, file = "binomial-DEGsTable.csv")

            return(list(sigDEG,DEGsTable))
                          
          }
          )

#' @title title
#' @description description
#' @param object object
#' @param fdr fdr
#' @importFrom dplyr select
#' @rdname KMClustDiffGenes
#' @export
setGeneric("KMClustDiffGenes", function(object,fdr=.01) standardGeneric("KMClustDiffGenes"))
#' @export
#' @rdname KMClustDiffGenes
setMethod("KMClustDiffGenes",
          signature = "PSCANseq",
          definition = function(object,fdr){
            if ( ! is.numeric(fdr) ) stop("pvalue has to be a number between 0 and 1") else if (  fdr < 0 | fdr > 1 ) stop("fdr has to be a number between 0 and 1")
            cdiff <- list()
            x     <- object@ndata
            y     <- object@expdata[,names(object@ndata)]
            part  <- object@kmeans$kpart
            for ( i in 1:max(part) ){
              if ( sum(part == i) == 0 ) next
              m <- apply(x,1,mean)
              n <- if ( sum(part == i) > 1 ) apply(x[,part == i],1,mean) else x[,part == i]
              no <- if ( sum(part == i) > 1 ) median(apply(y[,part == i],2,sum))/median(apply(x[,part == i],2,sum)) else sum(y[,part == i])/sum(x[,part == i])
              m <- m*no
              n <- n*no
              pv <- binompval(m/sum(m),sum(n),n)
              d <- data.frame(mean.all=m,mean.cl=n,fc=n/m,pv=pv)[order(pv,decreasing=FALSE),]
              #cdiff[[paste("cl",i,sep=".")]] <- d[d$pv < pvalue,]
              cdiff[[i]] <- d[d$pv < fdr,]
            }
            DEGsE<-c()
            DEGsS<-c()
            DEGsTable<-data.frame()

            for (n in 1:K){
                p.adj<- p.adjust(cdiff[[n]][,4],method="bonferroni")
                out<-cbind(cdiff[[n]], p.adj)
                out<-subset(out,out[,5]<fdr)
                Regulation<-c()
                for (i in 1:length(out[,1])){
                    if(out[i,1]>out[i,2]){
                        Regulation[i]="Down"
                    }else{
                        Regulation[i]="Up"
                    }
                }
                out<-cbind(out,Regulation)
                mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
                genes <- rownames(out)
                G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
                Final<-cbind(genes,out)
                Final<-merge(Final,G_list,by.x="genes",by.y="ensembl_gene_id")
                Final<-Final[!duplicated(Final[,8]), ]
                rownames(Final)<-Final[,1]
                Final[,1]<-Final[,8]
                Final<-Final[,-8]
    
                DEGsS<-c(DEGsS,Final[,1])
                DEGsE<-c(DEGsE,as.character(rownames(Final)))
                
                Up<-subset(Final,Final[,7]=="Up")
                Up<-select(Up, "Regulation","genes","pv","mean.all", "mean.cl","fc","p.adj")
                Up[,3]<-rownames(Up)
                Up[,6]<-log2(Up[,6])
                Up[,1]<-Up[,2]
                colnames(Up)<-c("Genes","genes","E.genes","mean.all", "mean.cl","log2.fc","p.adj")
                write.csv(Up, file = paste0("Up-DEG-cluster",n,".csv"))

                Down<-subset(Final,Final[,7]=="Down")
                Down<-select(Down, "Regulation","genes","pv","mean.all", "mean.cl","fc","p.adj")
                Down[,3]<-rownames(Down)
                Down[,6]<-log2(Down[,6])
                Down[,1]<- Down[,2]
                colnames(Down)<-c("Genes","genes","E.genes","mean.all", "mean.cl","log2.fc","p.adj")
                write.csv(Down, file = paste0("Down-DEG-cluster",n,".csv"))
    
                sigDEG<-cbind(DEGsE,DEGsS)
                write.csv(sigDEG, file = "binomial-sigDEG.csv")
                
                DEGsTable[n,1]<-paste0("Cluster ",n)
                DEGsTable[n,2]<-"Remaining Clusters"
                DEGsTable[n,3]<-length(Up[,1])
                DEGsTable[n,4]<-paste0("Up-DEG-cluster",n,".csv")
                DEGsTable[n,5]<-length(Down[,1])
                DEGsTable[n,6]<- paste0("Down-DEG-cluster",n,".csv")
            }
            colnames(DEGsTable)<-c("Target Cluster","VS","Up-regulated genes","File name","Low-regulated genes","File name")
            write.csv(DEGsTable, file = "binomial-DEGsTable.csv")

            return(list(sigDEG,DEGsTable))
                          
          }
          )