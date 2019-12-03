PSCANseq <- setClass("PSCANseq ", slots = c(expdata   = "data.frame", ndata = "data.frame", fdata = "data.frame", 
                                    distances  = "matrix", tsne = "data.frame", background = "list", out = "list", 
                                    cpart      = "vector", fcol = "vector", filterpar = "list", clusterpar = "list", 
                                    outlierpar = "list", kmeans = "list",MBclusters = "vector",kordering = "vector",
                                    MBordering = "vector" , MBtsne = "data.frame"))

setValidity("PSCANseq ",
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

setMethod("initialize",
          signature = "PSCANseq ",
          definition = function(.Object, expdata ){
            .Object@expdata <- expdata
            .Object@ndata <- expdata
            .Object@fdata <- expdata
            validObject(.Object)
            return(.Object)
          }
          )




setGeneric("Normalizedata", function(object, mintotal=1000, minexpr=0, minnumber=0, maxexpr=Inf, downsample=FALSE, dsn=1, rseed=17000) standardGeneric("Normalizedata"))
setMethod("Normalizedata",
          signature = "PSCANseq ",
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


downsample <- function(x,n,dsn){
  x <- round( x[,apply(x,2,sum,na.rm=TRUE) >= n], 0)
  nn <- min( apply(x,2,sum) )
  for ( j in 1:dsn ){
    z  <- data.frame(GENEID=rownames(x))
    rownames(z) <- rownames(x)
    initv <- rep(0,nrow(z))
    for ( i in 1:dim(x)[2] ){
      y <- aggregate(rep(1,nn),list(sample(rep(rownames(x),x[,i]),nn)),sum)
      na <- names(x)[i]
      names(y) <- c("GENEID",na)
      rownames(y) <- y$GENEID
      z[,na] <- initv
      k <- intersect(rownames(z),y$GENEID)
      z[k,na] <- y[k,na]
      z[is.na(z[,na]),na] <- 0
    }
    rownames(z) <- as.vector(z$GENEID)
    ds <- if ( j == 1 ) z[,-1] else ds + z[,-1]
  }
  ds <- ds/dsn + .1
  return(ds)
}

dist.gen <- function(x,method="euclidean", ...) if ( method %in% c("spearman","pearson","kendall") ) as.dist( 1 - cor(t(x),method=method,...) ) else dist(x,method=method,...)

dist.gen.pairs <- function(x,y,...) dist.gen(t(cbind(x,y)),...)

clustfun <- function(x,clustnr=20,bootnr=50,metric="pearson",do.gap=TRUE,SE.method="Tibs2001SEmax",SE.factor=.25,B.gap=50,cln=0,rseed=17000)
{
  if ( clustnr < 2) stop("Choose clustnr > 1")
  di <- dist.gen(t(x),method=metric)
  if ( do.gap | cln > 0 ){
    gpr <- NULL
    if ( do.gap ){
      set.seed(rseed)
      gpr <- clusGap(as.matrix(di), FUN = kmeans, K.max = clustnr, B = B.gap) 
      if ( cln == 0 ) cln <- maxSE(gpr$Tab[,3],gpr$Tab[,4],method=SE.method,SE.factor)
    }    
    if ( cln <= 1 ) {
      clb <- list(result=list(partition=rep(1,dim(x)[2])),bootmean=1)
      names(clb$result$partition) <- names(x)
      return(list(x=x,clb=clb,gpr=gpr,di=di))
    }
    clb <- clusterboot(di,B=bootnr,distances=FALSE,bootmethod="boot",clustermethod=KmeansCBI,krange=cln,scaling=FALSE,multipleboot=FALSE,bscompare=TRUE,seed=rseed)
    return(list(x=x,clb=clb,gpr=gpr,di=di))
  }
}


setGeneric("Clustexp", function(object,clustnr=20,bootnr=50,metric="pearson",do.gap=TRUE,SE.method="Tibs2001SEmax",SE.factor=.25,B.gap=50,cln=0,rseed=17000) standardGeneric("Clustexp"))
setMethod("Clustexp",
          signature = "PSCANseq ",
          definition = function(object,clustnr,bootnr,metric,do.gap,SE.method,SE.factor,B.gap,cln,rseed) {
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
            object@clusterpar <- list(clustnr=clustnr,bootnr=bootnr,metric=metric,do.gap=do.gap,SE.method=SE.method,SE.factor=SE.factor,B.gap=B.gap,cln=cln,rseed=rseed)
            y <- clustfun(object@fdata,clustnr,bootnr,metric,do.gap,SE.method,SE.factor,B.gap,cln,rseed)
            object@kmeans    <- list(kpart=y$clb$result$partition, jaccard=y$clb$bootmean, gap=y$gpr)
            object@distances <- as.matrix( y$di )
            set.seed(111111)
            object@fcol <- sample(rainbow(max(y$clb$result$partition)))
            return(object)
          }
          )



setGeneric("plotGap", function(object) standardGeneric("plotGap"))
setMethod("plotGap",
          signature = "PSCANseq ",
          definition = function(object){
            if ( length(object@kmeans$kpart) == 0 ) stop("run clustexp before plotgap")
            plot(object@kmeans$gap,ylim=c(0.1,0.5),las=1,main="Gap Statistics")
          }
          )




setGeneric("comptSNE", function(object,rseed=15555) standardGeneric("comptSNE"))
setMethod("comptSNE",
          signature = "PSCANseq ",
          definition = function(object,rseed){
            if ( length(object@kmeans$kpart) == 0 ) stop("run clustexp before comptsne")
            set.seed(rseed)
            di <- dist.gen(as.matrix(object@distances))
            ts <- tsne(di,k=2)
            object@tsne <- as.data.frame(ts)
            return(object)
          }
          )






            

#Silhouette 
setGeneric("plotSilhouette", function(object,K) standardGeneric("plotSilhouette"))
setMethod("plotSilhouette",
          signature = "PSCANseq ",
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




#Jaccard 
Jaccard<- function(object,Clustering) {
    JACCARD<-c()
    if ( ! ( Clustering %in% c( "K-means","MB") ) ) stop("Clustering has to be either K-means or MB")
    JS <- function(data, indices) {
      d <- data[indices,]                                    # allows boot to select sample 
      jac<-distance(t(d), method = "jaccard")
      jac1<-1- jac
      JSmean<- mean(jac1)
      return(JSmean)
    }
    for (i in 1:K){
        if (Clustering=="K-means"){
            results <- boot(data=object@fdata[,which(object@kmeans$kpart==i)], statistic=JS, R=100,stype = "f")
            JACCARD[i]<-round(mean(results$t),digits=3)            # to get the mean of all bootstrappings (mean of mean Jaccard values)
    
        }
        if (Clustering=="MB"){
            results <- boot(data=object@fdata[,which(object@MBclusters$clusterid==i)], statistic=JS, R=100,stype = "f")
            JACCARD[i]<-round(mean(results$t),digits=3)            # to get the mean of all bootstrappings (mean of mean Jaccard values)
    
        }
    }
    barplot(JACCARD,names.arg=1:length(JACCARD),ylab="Mean Jaccard's similarity values",xlab="Clusters",
    las=1,ylim=c(0,1),col=c("black","blue","green","red","yellow","gray"))
    box()
    return(JACCARD)
}


setGeneric("plottSNE", function(object) standardGeneric("plottSNE"))
setMethod("plottSNE",
          signature = "PSCANseq ",
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




setGeneric("plotKmeansLabelstSNE", function(object) standardGeneric("plotKmeansLabelstSNE"))
setMethod("plotKmeansLabelstSNE",
          signature = "PSCANseq ",
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



setGeneric("plotSymbolstSNE", function(object,types=NULL) standardGeneric("plotSymbolstSNE"))
setMethod("plotSymbolstSNE",
          signature = "PSCANseq ",
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



PCAplotSymbols= function(object1,object2,types=NULL){
	types <- names(object2)
	types<- gsub("_[0-9]+","",types)
	coloc <- rainbow(length(unique(types)))
	syms<-c()
	plot(object1[,1],object1[,2],xlab="PC1",ylab="PC2",pch=20,cex=0,col="grey",las=1)
    for ( i in 1:length(unique(types)) ){
		f <- types == sort(unique(types))[i]
        syms <- append( syms, ( (i-1) %% 25 ) + 1 )
		points(object1[f,1],object1[f,2],col=coloc[i],pch=( (i-1) %% 25 ) + 1,cex=1)
        }
    legend("topright", legend=sort(unique(types)), col=coloc, pch=syms)
}       




setGeneric("KMclustheatmap", function(object,hmethod="single") standardGeneric("KMclustheatmap"))

setMethod("KMclustheatmap",
          signature = "PSCANseq ",
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
            axis(1,at=tmp,lab=cclmo)
            axis(2,at=tmp,lab=cclmo)
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



Exprmclust <- function (data, clusternum = 2:9, modelNames = "VVV", reduce = T, cluster = NULL) {
      set.seed(12345)
      if (reduce) {
            sdev <- prcomp(t(data), scale = T)$sdev[1:20]
            x <- 1:20
            optpoint <- which.min(sapply(2:10, function(i) {
                  x2 <- pmax(0, x - i)
                  sum(lm(sdev ~ x + x2)$residuals^2)
            }))
            pcadim = optpoint + 1
            tmpdata <- t(apply(data, 1, scale))
            colnames(tmpdata) <- colnames(data)
            tmppc <- prcomp(t(tmpdata), scale = T)
            pcareduceres <- t(tmpdata) %*% tmppc$rotation[, 1:pcadim]
      }
      else {
            pcareduceres <- t(data)
      }
      if (is.null(cluster)) {   
            clusternum <- clusternum[clusternum > 1]
            res <- suppressWarnings(Mclust(pcareduceres, G = clusternum, modelNames = modelNames))
            clusterid <- apply(res$z, 1, which.max)
            clunum <- res$G
      } else {
            clunum <- length(unique(cluster))
            clusterid <- cluster
      }
      clucenter <- matrix(0, ncol = ncol(pcareduceres), nrow = clunum)
      for (cid in 1:clunum) {
            clucenter[cid, ] <- colMeans(pcareduceres[names(clusterid[clusterid == cid]), , drop = F])
      }
      dp <- as.matrix(dist(clucenter))
      gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
      dp_mst <- minimum.spanning.tree(gp)
      list(pcareduceres = pcareduceres, MSTtree = dp_mst, clusterid = clusterid, clucenter = clucenter)
}




PlotMBorderPCA<- function(object) {
    l <- object[,3]
    mi <- min(l,na.rm=TRUE)
    ma <- max(l,na.rm=TRUE)
    ColorRamp <- colorRampPalette(rev(brewer.pal(n = 11,name = "RdYlBu")))(100)
    ColorLevels <- seq(mi, ma, length=length(ColorRamp))
    v <- round((l - mi)/(ma - mi)*99 + 1,0)
    layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
    par(mar = c(5,5,2.5,2))
    plot(object[,1],object[,2],xlab="PC1",ylab="PC2",pch=20,cex=0,col="grey",las=1)
    for ( k in 1:length(v) ){
        points(object[k,1],object[k,2],col=ColorRamp[v[k]],pch=20,cex=2)
    }
    par(mar = c(3,2.5,2.5,2))
    image(1, ColorLevels,matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),col=ColorRamp,xlab="",ylab="",las=2,xaxt="n")
    layout(1)
}





setGeneric("MBclustheatmap", function(object,hmethod="single") standardGeneric("MBclustheatmap"))

setMethod("MBclustheatmap",
          signature = "PSCANseq ",
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
            axis(1,at=tmp,lab=cclmo)
            axis(2,at=tmp,lab=cclmo)
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


setGeneric("plotExptSNE", function(object,g,n="") standardGeneric("plotExptSNE"))
setMethod("plotExptSNE",
          signature = "PSCANseq ",
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




Plotmclust <- function(mclustobj, x = 1, y = 2, MSTorder = NULL, show_tree = T, show_full_tree = F, show_cell_names = F, cell_name_size = 3, markerexpr = NULL, showcluster = T) {
      color_by = "State"
      lib_info_with_pseudo <- data.frame(State=mclustobj$clusterid,sample_name=names(mclustobj$clusterid))
      lib_info_with_pseudo$State <- factor(lib_info_with_pseudo$State)
      S_matrix <- mclustobj$pcareduceres
      pca_space_df <- data.frame(S_matrix[,c(x, y)])
      colnames(pca_space_df) <- c("pca_dim_1","pca_dim_2")
      pca_space_df$sample_name <- row.names(pca_space_df)
      edge_df <- merge(pca_space_df, lib_info_with_pseudo, by.x = "sample_name", by.y = "sample_name")     
      edge_df$Marker <- markerexpr[edge_df$sample_name]
      if (!is.null(markerexpr)) {
            g <- ggplot(data = edge_df, aes(x = pca_dim_1, y = pca_dim_2, size = Marker))
            if (showcluster) {
                  g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE) 
                    g <- g +scale_colour_manual(values = c("1"="black","2"="blue","3"="green","4"="red","5"="yellow","6"="gray"))

            } else {
                  g <- g + geom_point(na.rm = TRUE,color="green")
            }
      } else {
            g <- ggplot(data = edge_df, aes(x = pca_dim_1, y = pca_dim_2))
            if (showcluster) {
                  g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE, size = 3)
                    g <- g +scale_colour_manual(values = c("1"="black","2"="blue","3"="green","4"="red","5"="yellow","6"="gray"))

      
            } else {
                  g <- g + geom_point(na.rm = TRUE, size = 3)
            }
      }
      if (show_cell_names) {
            g <- g + geom_text(aes(label = sample_name), size = cell_name_size)
      }
      
      if (show_tree) {
            clucenter <- mclustobj$clucenter[,c(x,y)]
            clulines <- NULL
            if (show_full_tree) {
                  alledges <- as.data.frame(get.edgelist(mclustobj$MSTtree),stringsAsFactors=F)
                  alledges[,1] <- as.numeric(alledges[,1])
                  alledges[,2] <- as.numeric(alledges[,2])
                  for (i in 1:nrow(alledges)) {
                        clulines <- rbind(clulines, c(clucenter[alledges[i,1],],clucenter[alledges[i,2],]))
                  }      
            } else {
                  if (is.null(MSTorder)) {
                        clutable <- table(mclustobj$clusterid)
                        alldeg <- degree(mclustobj$MSTtree)
                        allcomb <- expand.grid(as.numeric(names(alldeg)[alldeg == 
                                                                              1]), as.numeric(names(alldeg)[alldeg == 1]))
                        allcomb <- allcomb[allcomb[, 1] < allcomb[, 2], ]
                        numres <- t(apply(allcomb, 1, function(i) {
                              tmp <- as.vector(get.shortest.paths(mclustobj$MSTtree, 
                                                                  i[1], i[2])$vpath[[1]])
                              c(length(tmp), sum(clutable[tmp]))
                        }))
                        optcomb <- allcomb[order(numres[, 1], numres[, 2], decreasing = T)[1], ]
                        MSTorder <- get.shortest.paths(mclustobj$MSTtree, optcomb[1], 
                                                       optcomb[2])$vpath[[1]]
                  }
                  for (i in 1:(length(MSTorder)-1)) {
                        clulines <- rbind(clulines, c(clucenter[MSTorder[i],],clucenter[MSTorder[i+1],]))
                  }      
            }
            clulines <- data.frame(x=clulines[,1],xend=clulines[,3],y=clulines[,2],yend=clulines[,4])
            g <- g + geom_segment(aes_string(x="x",xend="xend",y="y",yend="yend",size=NULL),data=clulines,size=1,color ="orange")
            
            clucenter <- data.frame(x=clucenter[,1],y=clucenter[,2],id=1:nrow(clucenter))
            g <- g + geom_text(aes_string(label="id",x="x",y="y",size=NULL),data=clucenter,size=10,color ="orange")
            
      }            
      g <- g + guides(colour = guide_legend(override.aes = list(size=5))) + 
            xlab(paste0("PCA_dimension_",x)) + ylab(paste0("PCA_dimension_",y)) +
            theme(panel.border = element_blank(), axis.line = element_line()) + 
            theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
            theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
            theme(legend.position = "top", legend.key.size = unit(0.3, "in"),legend.text = element_text(size = 20),legend.title=element_text(size = 20),legend.box = "vertical") + theme(legend.key = element_blank()) + 
            theme(panel.background = element_rect(fill = "white")) +
            theme(axis.text.x = element_text(size=17,color="black"),
                  axis.text.y = element_text(size=17,color='black'),
                  axis.title.x = element_text(size=20,vjust=-1),
                  axis.title.y = element_text(size=20,vjust=1),
                  plot.margin=unit(c(1,1,1,1),"cm"))       
      g       
}
   


setGeneric("comptsneMB", function(object,rseed=15555) standardGeneric("comptsneMB"))
setMethod("comptsneMB",
          signature = "PSCANseq ",
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
           

setGeneric("plottsneMB", function(object,K) standardGeneric("plottsneMB"))
setMethod("plottsneMB",
          signature = "PSCANseq ",
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
               

setGeneric("plotMBLabelstSNE", function(object) standardGeneric("plotMBLabelstSNE"))
setMethod("plotMBLabelstSNE",
          signature = "PSCANseq ",
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




setGeneric("plotsilhouetteMB", function(object,K) standardGeneric("plotsilhouetteMB"))

setMethod("plotsilhouetteMB",
          signature = "PSCANseq ",
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

setGeneric("plotexptsneMB", function(object,g,n="") standardGeneric("plotexptsneMB"))

setMethod("plotexptsneMB",
          signature = "PSCANseq ",
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




add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}


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

cross.val <- function(exp.df, class.vec, segments, performance, class.algo){
	
	#Start cross validation loop
	class1 <- levels(class.vec)[1]
	for(fold in 1:length(segments)){
		cat("Fold", fold, "of", length(segments), "\n")

		#Define training and test set
		test.ind <- segments[[fold]]
		training.set <- exp.df[-test.ind,]
		training.class <- class.vec[-test.ind]
		test.set <- exp.df[test.ind,, drop=FALSE]
		test.class <- class.vec[test.ind]
		
		#Train J48 on training set
		if(class.algo == "J48"){
			cv.model <- J48(training.class ~ ., training.set)
            pred.class <- predict(cv.model, test.set)
		} else if(class.algo == "rpart"){
			cv.model <- rpart(training.class ~ ., training.set,method="class")
		pred.class <- predict(cv.model, test.set,type="class")
		} else{
			stop("Unknown classification algorithm")
		}

		#Evaluate model on test set
#		pred.class <- predict(cv.model, test.set)
#		pred.class <- predict(cv.model, test.set,type="class")
		performance <- eval.pred(pred.class, test.class, class1, performance)
	}
	return(performance)

}

#Function for counting TPs, FNs, FPs and TNs
eval.pred <- function(pred.class, true.class, class1, performance){
	for(index in 1:length(pred.class)){
		pred <- pred.class[index]
		true <- true.class[index]
		if(pred == true && true == class1){
			performance["TP"] <- performance["TP"] + 1
		} else if(pred != true && true == class1){
			performance["FN"] <- performance["FN"] + 1
		} else if(pred != true && true != class1){
			performance["FP"] <- performance["FP"] + 1
		} else if(pred == true && true != class1){
			performance["TN"] <- performance["TN"] + 1
		}
	}
	return(performance)
}
    
    
SN <- function(con.mat){
    TP <- con.mat[1,1]
    FN <- con.mat[2,1]
    return(TP/(TP+FN))
}
SP <- function(con.mat){
    TN <- con.mat[2,2]
    FP <- con.mat[1,2]
    return(TN/(TN+FP))
}
ACC <- function(con.mat){
    TP <- con.mat[1,1]
    FN <- con.mat[2,1]
    TN <- con.mat[2,2]
    FP <- con.mat[1,2]
    return((TP+TN)/(TP+FN+TN+FP))
}
MCC <- function(con.mat){
    TP <- con.mat[1,1]
    FN <- con.mat[2,1]
    TN <- con.mat[2,2]
    FP <- con.mat[1,2]
    denom <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    denom <- ifelse(denom==0, NA, denom)
    return((TP*TN-FP*FN)/denom)
}









J48DT<-function(object){
	exp.df<-as.data.frame(t(object))
	classVector<- factor(colnames(object))
	j48.model<-J48(classVector~.,exp.df)
	print(j48.model)
      plot(as.party(j48.model),gp = gpar(cex=0.65,col="black", lty = "solid", lwd = 1.5, fontsize = 12))
	return(j48.model)
}

J48DTeval<- function(object,num.folds=10,First,Second){
	exp.imput.df<-as.data.frame(t(object))
	num.instances<-nrow(exp.imput.df)
	indices<-1:num.instances
	classVector<- factor(colnames(object))

	cv.segments<-split(sample(indices),rep(1:num.folds,length=num.instances))
	j48.performance<-c("TP"=0,"FN"=0,"FP"=0,"TN"=0)
	j48.performance<-cross.val(exp.imput.df,classVector,cv.segments,j48.performance,"J48")
	print(j48.performance)

	j48.confusion.matrix<-matrix(j48.performance,nrow=2)
	rownames(j48.confusion.matrix)<-c(paste0("Predicted",First), paste0("Predicted",Second))
	colnames(j48.confusion.matrix)<-c(First,Second)
	print(j48.confusion.matrix)

	j48.sn<-SN(j48.confusion.matrix)
	j48.sp<-SP(j48.confusion.matrix)
	j48.acc<-ACC(j48.confusion.matrix)
	j48.mcc<-MCC(j48.confusion.matrix)

	cat("J48 SN: ", j48.sn, "\n",
		"J48 SP: ", j48.sp, "\n",
		"J48 ACC: ", j48.acc, "\n",
		"J48 MCC: ", j48.mcc, "\n",sep="")
	return(j48.performance)
}







RpartDT<-function(object){
	exp.df<-as.data.frame(t(object))
	classVector<- factor(colnames(object))
	model<-rpart(classVector~.,exp.df,method="class",minsplit = 1, cp=-1)
	print(model)
	rpart.plot(model,type=4,extra=101)
	return(model)
}


RpartEVAL<- function(object,num.folds=10,First,Second){
	exp.imput.df<-as.data.frame(t(object))
	num.instances<-nrow(exp.imput.df)
	indices<-1:num.instances
	classVector<- factor(colnames(object))

	cv.segments<-split(sample(indices),rep(1:num.folds,length=num.instances))
	Rpart.performance<-c("TP"=0,"FN"=0,"FP"=0,"TN"=0)
	Rpart.performance<-cross.val(exp.imput.df,classVector,cv.segments,Rpart.performance,"rpart")
	print(Rpart.performance)
	Rpart.confusion.matrix<-matrix(Rpart.performance,nrow=2)
	rownames(Rpart.confusion.matrix)<-c(paste0("Predicted",First), paste0("Predicted",Second))
	colnames(Rpart.confusion.matrix)<-c(First,Second)
	print(Rpart.confusion.matrix)

	Rpart.sn<-SN(Rpart.confusion.matrix)
	Rpart.sp<-SP(Rpart.confusion.matrix)
	Rpart.acc<-ACC(Rpart.confusion.matrix)
	Rpart.mcc<-MCC(Rpart.confusion.matrix)

	cat("Rpart SN: ", Rpart.sn, "\n",
	"Rpart SP: ", Rpart.sp, "\n",
	"Rpart ACC: ", Rpart.acc, "\n",
	"Rpart MCC: ", Rpart.mcc, "\n",sep="")

	return(Rpart.performance)
}








PPI<-function(data,FileName){
	# Save base enpoint as variable
	string_api_url <- "https://string-db.org/api/"
	output_format <- "tsv" #"json", "tsv-no-header", "tsv", "xml"
	method <- "network"
	species <- "9606"
	your_identifiers <- ""
	optional_parameters <- ""
	# Construct API request
	genes <- data
	repos <- GET(url = paste0(string_api_url,output_format,'/',method,'?identifiers=',
                          paste(as.character(data), collapse = "%0d"),"&", "species=",species))
    cat("Examine response components =",status_code(repos),"\t","200 means successful","\n")
	# Process API request content 
	repo_content <- content(repos)
	results  <- read_tsv(repo_content)
	write.csv(results, file = paste0("PPI-",FileName,".csv"))
	return(results)
}



Networking<-function(data,FileName){
	string_api_url <- "https://string-db.org/api/"
	output_format <- "image"
	method <- "network"
	species <- "9606"
	your_identifiers <- ""
	optional_parameters <- ""

	# Construct API request
	genes <- data
	repos <- GET(url = paste0(string_api_url,output_format,'/',method,'?identifiers=',
                          paste(as.character(data), collapse = "%0d"),"&", "species=",species))
      cat("Examine response components =",status_code(repos),"\t","200 means successful","\n")

	y = repos$request$ url
	download.file(y,paste0("network",FileName,".png"), mode = 'wb')
	Network<- readPNG(paste0("network",FileName,".png"), native = TRUE)
	plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
	rasterImage(Network,0,0,1,1)	
}


NetAnalysis<-function(object){
    if ( length(object[,1])<1 ) stop( "No Protein-Protein Interactions" )
    df<-object[,-c(1,2)]
    gg <- graph.data.frame(df)
    betweenness<-betweenness(gg)
    betweenness.table<- data.frame(betweenness)
    names <- rownames(betweenness.table)
    rownames(betweenness.table) <- NULL
    degree<-degree(gg)
    degree.table<- data.frame(degree)
    names <- rownames(degree.table)
    rownames(degree.table) <- NULL
    AnalysisTable <- cbind(names,degree.table,betweenness.table)
    write.csv(AnalysisTable, file = paste0("NetworkAnalysisTable-",FileName,".csv"))

    test.graph.adj<-get.adjacency(gg,sparse=F)
    test.graph.properties<-GenInd(test.graph.adj)
    cat("Number of nodes: ",test.graph.properties$N,"\n")
    V(gg)
    cat("Number of links: ",test.graph.properties$Ltot,"\n")
    E(gg)
    cat("Link Density: ",test.graph.properties$LD,"\n")
    cat("The connectance of the graph: ",test.graph.properties$C,"\n")
    cat("Mean Distences",mean_distance(gg),"\n")
    cat("Average Path Length",average.path.length(gg),"\n","\n")
    return(AnalysisTable)    
}





NoiseFiltering<-function(object,percentile,CV,GeneList,geneCol,FgeneCol,erccCol,Val=TRUE){
    geneTypes <- factor( c( ENSG="ENSG", ERCC="ERCC" )[ substr( rownames(object), 1, 4 ) ] )    # Split data into sub tables based on the factor object geneTypes
    countsG1ms   <- valuesG1ms[ which( geneTypes=="ENSG" ), ]                                   # calculate normalisation for counts\n",
    countsERCC  <- valuesG1ms[ which( geneTypes=="ERCC" ), ]
    sfERCC <- estimateSizeFactorsForMatrix( countsERCC )
    # Given a matrix or data frame of count data, this function estimates the size factors as follows: Each column is divided by the geometric means of the rows. The median (or, 
    # If requested, another location estimator) of these ratios (skipping the genes with a # geometric mean of zero) is used as the size factor for this column. Source: DESeq package.
    sfERCC <- estimateSizeFactorsForMatrix( countsERCC )
    sfG1ms <-  estimateSizeFactorsForMatrix( countsG1ms )
    nCountsERCC <- t( t(countsERCC) / sfERCC )               # Divide columns by size factors to normalize counts
    nCountsG1ms <- t( t(countsG1ms) / sfG1ms )
    # perform fit, define sample moments per gene
    meansG1ms <- rowMeans(nCountsG1ms)
    varsG1ms <- rowVars(nCountsG1ms)
    cv2G1ms <- varsG1ms / meansG1ms^2
    meansERCC <- rowMeans(nCountsERCC)
    varsERCC <- rowVars(nCountsERCC)
    cv2ERCC <- varsERCC / meansERCC^2
    minMeanForFit <- unname( quantile( meansERCC[ which( cv2ERCC > .3 ) ], percentile ) )
    cat("Cut-off value for the ERCCs= ",round(minMeanForFit,digits=2),"\n","\n")
    #Perform the fit of technical noise strength on average count. We regress cv2HeLa on 1/meansForHeLa. We use the
    #glmgam.fit function from the statmod package to perform the regression as a GLM fit of the gamma family with log link.
    #The 'cbind' construct serves to produce a model matrix with an intercept.
    useForFit <- meansERCC >= minMeanForFit
    fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansERCC[useForFit] ), cv2ERCC[useForFit] )
    cat("Coefficients of the fit:","\n")
    print(fit$coefficients)
    table( useForFit )
    #To get the actual noise coefficients, we need to subtract Xi
    xi <- mean( 1 / sfERCC )
    a0 <- unname( fit$coefficients["a0"] )
    a1 <- unname( fit$coefficients["a1tilde"] - xi )
    #cat("\n","The actual noise coefficients: ",c( a0, a1 ),"\n")
    #how much variance does the fit explain?
    residual <- var( log( fitted.values(fit) ) - log( cv2ERCC[useForFit] ) )
    total <- var( log( cv2ERCC[useForFit] ) )
    cat("Explained variances of log CV^2 values= ",c(round(1 - residual / total,digits=2)),"\n","\n")
    ## Pick out genes above noise line
    # test which entries are above the line
    idx_test <- cv2G1ms>(xi+a1)/meansG1ms + a0
    # pick out genes that fulfil statement
    genes_test <- gene_names2[idx_test] #pick out genes
    genes_test <- genes_test[!is.na(genes_test)] #remove na entries
    meansG1ms_test <- meansG1ms[idx_test] #take out mean values for fulfilled genes
    meansG1ms_test <- meansG1ms_test[!is.na(meansG1ms_test)] #remove na entries
    cv2G1ms_test <- cv2G1ms[idx_test] #take out cv2 values for fulfilled genes
    cv2G1ms_test <- cv2G1ms_test[!is.na(cv2G1ms_test)] #remove na entries
    genes_test = genes_test[-which(sapply(genes_test, is.null))]
    genes_test <- sapply( genes_test, paste0, collapse="")
    cat("Number of genes that passt the filtering= ",length(genes_test),"\n","\n")
    write.csv(genes_test, file = "Noise_filtering_genes_test.csv")

    cat("The filtered gene list was saved as: Noise_filtering_genes_test","\n")
    ## test genes for variance, the following is the term Psi + a0 * Theta, that appears in the formula for Omega.
    psia1theta <- mean( 1 / sfG1ms,na.rm=TRUE ) + a1 * mean( sfERCC / sfG1ms,na.rm=TRUE )
    # Now, we perform a one-sided test against the null hypothesis that the true variance is at most the technical variation plus biological variation with a CV of at most 50% (minBiolDisp = .52).
    minBiolDisp <- 0.5^2
    #Calculate Omega, then perform the test, using the formula given in the Online methods and in Supplementary Note 6.
    m <- ncol(countsG1ms)
    cv2th <- a0 + minBiolDisp + a0 * minBiolDisp
    testDenom <- ( meansG1ms * psia1theta + meansG1ms^2 * cv2th ) / ( 1 + cv2th/m )
    p <- 1 - pchisq( varsG1ms * (m-1) / testDenom, m-1 )
    p <- subset(p, !is.nan(p))

    padj <- p.adjust( p, "BH" )       #Adjust for multiple testing with the Benjamini-Hochberg method, cut at 10%
    plot( NULL, xaxt="n", yaxt="n",log="xy", xlim = c( 1e-1, 3e5 ), ylim = c( .005, 100 ),main="Gene filtration by accounting for technical noise",
    xlab = "Average normalized read count", ylab = "Squared coefficient of variation (CV^2)" )
    axis( 1, 10^(-1:5), c("0.1", "1", "10", "100", "1000", expression(10^4), expression(10^5) ) )
    axis( 2, 10^(-2:2), c("0.01", "0.1", "1", "10", "100" ), las=2 )
    abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0", lwd=2 )
    # Plot the genes, use a different color if they are highly variable
    points( meansG1ms, cv2G1ms, pch=20, cex=.2,col = geneCol )

    #highlight gene list from test
	 points( meansG1ms_test, cv2G1ms_test, pch=20, cex=.22,col= FgeneCol)

    # Add the technical noise fit, as before
    xg <- 10^seq( -2, 6, length.out=1000 )
    lines( xg, (xi+a1)/xg + a0, col="red", lwd=5 )
    
    # Add the normalised ERCC points
    if (Val){
    	points( meansERCC[useForFit], cv2ERCC[useForFit], pch=20, cex=1.5, col=erccCol ) # Showing only the valied ERCCs
    }else{
    	points( meansERCC, cv2ERCC, pch=20, cex=2, col=erccCol) # Showing all the valied ERCCs
    }
    add_legend("topleft", legend=c("Noise curve","ERCC spike-ins","Genes above the noise line"), pch=c(15,20,20), 
     col=c("red",erccCol,FgeneCol),horiz=TRUE, bty='n', cex=0.85)

    return (genes_test)
}



KmeanOrder<-function(object,Clusters,sampleNames){
	lpsmclust <- Exprmclust(object,clusternum =4,reduce = F, cluster = Clusters)
	lpsorder <- TSCANorder(lpsmclust)
	orderID<-lpsorder
	order<-c(1:length(lpsorder))
	orderTable<-data.frame(order,orderID)
	write.csv(orderTable, file = "Cellular_pseudo-time_ordering_based_on_k-meansc-lusters.csv")
	print(orderTable)
	FinalOrder<-orderTable[match(sampleNames, orderTable$orderID),]
	sc@kordering<-FinalOrder[,1]
	names(sc@kordering)<-names(Clusters)
	return(sc@kordering)
}




MB_Order<-function(object,sampleNames,Names){
	lpsorderMB <- TSCANorder(object)
	orderID<-lpsorderMB
	order<-c(1:length(lpsorderMB))
	orderTableMB<-data.frame(order,orderID)
	write.csv(orderTableMB, file = "Cellular_pseudo-time_ordering_based_on_Model-based_clusters.csv")
	print(orderTableMB)
	FinalOrder<-orderTableMB[match(sampleNames, orderTableMB$orderID),]
	sc@MBordering<-FinalOrder[,1]
	names(sc@MBordering)<-names(Names)
	return(sc@MBordering)
}

ClassVectoringDT<- function(object,Cluster_ID,K,First,Second,sigDEG){
    SC<-PSCANseq(object)
    SC<- Normalizedata(SC)
    DatasetForDT<- SC@fdata
    Nam <-colnames(DatasetForDT)
    num<-c(1:K)
    num1<- paste("CL", num, sep="")
    for (n in num){
        Nam<-ifelse((Cluster_ID==n),num1[n],Nam)
    }
    colnames(DatasetForDT)<-Nam

    #cat("The dataset is ready for decision tree analysis","\n")
    #dim(DatasetForDT)
    sg1 <- DatasetForDT[,which(colnames(DatasetForDT)==First | colnames(DatasetForDT)==Second)]
    dim(sg1)
    # Creating a dataset that includes only the DEGs
    gene_list<-sigDEG[,1]
    gene_names<-rownames(DatasetForDT)
    idx_genes <- is.element(gene_names, gene_list)
    gene_names2 <- gene_names[idx_genes]
    DEGsfilteredDataset<- sg1[gene_names2,]
    cat("The DEGs filtered normalized dataset contains:","\n","Genes:",length(DEGsfilteredDataset[,1]),"\n","cells:",length(DEGsfilteredDataset[1,]))
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    
    genes <- rownames(DEGsfilteredDataset)
    G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
    
    DATAforDT<-cbind(genes,DEGsfilteredDataset)
    DATAforDT<-merge(DATAforDT,G_list,by.x="genes",by.y="ensembl_gene_id")
    DATAforDT
    DATAforDT[,1]<-DATAforDT[,length(DATAforDT[1,])]
    DATAforDT<-DATAforDT[!duplicated(DATAforDT[,1]), ]

    rownames(DATAforDT)<-DATAforDT[,1]
    DATAforDT<-DATAforDT[,c(-1,-length(DATAforDT[1,]))]
    sg <- factor(gsub(paste0("(",First,"|",Second,").*"), "\\1", colnames(DATAforDT)), levels = c(paste0(First), paste0(Second)))
    sg<-sg[!is.na(sg)]
    colnames(DATAforDT)<-sg
    return(DATAforDT)
}





DEGanalysis<- function(object,Cluster_ID,K,fdr,name){
    Nam <-colnames(object)
    num<-c(1:K)
    num1<- paste("CL", num, sep="")
    for (n in num){
        Nam<-ifelse((Cluster_ID==n),num1[n],Nam)
    }
    colnames(object)<-Nam
    cat("The dataset is ready for differential expression analysis")
    
    num1<- paste("CL", num, sep="")
    clustName<- paste("Cl", num, sep="")
    ClusterLength<-c()
    for (i in num){
        d1<-clustName[i]
        d2<-object[,which(names(object)==num1[i])]
        ClusterLength[i]<-length(d2)
        assign(d1,d2)
    }
    ccdf<-data.frame(clustName,ClusterLength)
    ccdff<-ccdf[order(ClusterLength),] 
    clustName<-ccdff[,1]
    print(clustName)
    first<-c()
    second<-c()
    if ( K<2 ){
        stop( "K has to be at least 2" )
    }else if (K==2){
        first<-c(paste0(clustName[1]))
        second<-c(paste0(clustName[2]))

    }else if ( K==3 ){
        first<-c(rep(paste0(clustName[1]),2),rep(paste0(clustName[2]),1))
        second<-c(paste0(clustName[2]),paste0(clustName[3]),paste0(clustName[3]))
    }else if (K==4){
        first<-c(rep(paste0(clustName[1]),3),rep(paste0(clustName[2]),1),rep(paste0(clustName[4]),1),rep(paste0(clustName[3]),1))
        second<-c(paste0(clustName[2]),paste0(clustName[3]),paste0(clustName[4]),paste0(clustName[3]),paste0(clustName[2]),paste0(clustName[4]))
    }else if (K==5){
        first<-c(rep(paste0(clustName[1]),4),rep(paste0(clustName[2]),3),rep(paste0(clustName[3]),2),rep(paste0(clustName[5]),1))
        second<-c(paste0(clustName[2]),paste0(clustName[3]),paste0(clustName[4]),paste0(clustName[5]),paste0(clustName[3]),paste0(clustName[4]),paste0(clustName[5]),paste0(clustName[4]),paste0(clustName[5]),paste0(clustName[4]))
    }else if (K==6){
        first<-c(rep(paste0(clustName[1]),3),rep(paste0(clustName[5]),1),rep(paste0(clustName[1]),1),rep(paste0(clustName[2]),2),rep(paste0(clustName[5]),1),rep(paste0(clustName[2]),1),rep(paste0(clustName[3]),1),rep(paste0(clustName[5]),1),rep(paste0(clustName[3]),1),rep(paste0(clustName[5]),1),rep(paste0(clustName[4]),1),rep(paste0(clustName[5]),1))
        second<-c(paste0(clustName[2]),paste0(clustName[3]),paste0(clustName[4]),paste0(clustName[1]),paste0(clustName[6]),paste0(clustName[3]),paste0(clustName[4]),paste0(clustName[2]),paste0(clustName[6]),paste0(clustName[4]),paste0(clustName[3]),paste0(clustName[6]),paste0(clustName[4]),paste0(clustName[6]),paste0(clustName[6]))
     }
    o<-c(1:K)
    oo<-o[-length(o)]
    com<-sum(oo)
    cat("Number of comparisons: ", com *2,"\n")
    comNUM<-paste("comp", c(1:com), sep="")
    DEGsTable<-data.frame()
    DEGsE<-c()
    DEGsS<-c()
    for (i in 1:com){
        FinalDEGsL<-data.frame()
        FinalDEGsU<-data.frame()
        FDRl<-c()
	    FDRu<-c()

        d1<-comNUM[i]
        d2<-cbind(get(first[i]), get(second[i]))
        assign(d1,d2)
        len<-c(length(get(first[i])[1,]),length(get(second[i])[1,]))
        y <- c(rep(1:2,len))
        L<-as.matrix(get(comNUM[i]))
        gname<-rownames(get(comNUM[i]))
        x<-L
        data=list(x=x,y=y, geneid=gname)
        samr.obj<-samr(data, resp.type="Two class unpaired", assay.type="seq",nperms=100,nresamp=20,testStatistic="wilcoxon",random.seed=15)
        delta.table <- samr.compute.delta.table(samr.obj)
        
        wm<-which.min(delta.table[,5])
        if (delta.table[wm,5]<=fdr){
            w<-which(delta.table[,5]<= fdr)
            delta<-delta.table[w[1],1]-0.001
    
            samr.plot(samr.obj, delta)
            title(paste0("DEGs in the ",second[i]," in ",first[i]," VS ",second[i]))
    
            siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, data, delta.table)
        
            FDRl<-as.numeric(siggenes.table$genes.lo[,8])/100
            FDRu<-as.numeric(siggenes.table$genes.up[,8])/100

            siggenes.table$genes.lo[,8]<- FDRl
            siggenes.table$genes.up[,8]<- FDRu
    
            DEGsTable[i,1]<-paste0(first[i]," VS ",second[i])
            DEGsTable[i,2]<-second[i]
            DEGsTable[i,3]<-length(FDRu)
            DEGsTable[i,4]<-paste0("Up-regulated-",name,second[i],"in",first[i],"VS",second[i],".csv")
            DEGsTable[i,5]<-length(FDRl)
            DEGsTable[i,6]<-paste0("Low-regulated-",name,second[i],"in",first[i],"VS",second[i],".csv")

            DEGsTable[i+com,1]<-paste0(first[i]," VS ",second[i])
            DEGsTable[i+com,2]<-first[i]
            DEGsTable[i+com,3]<-length(FDRu)
            DEGsTable[i+com,4]<-paste0("Low-regulated-",name,first[i],"in",first[i],"VS",second[i],".csv")
            DEGsTable[i+com,5]<-length(FDRl)
            DEGsTable[i+com,6]<-paste0("Up-regulated-",name,first[i],"in",first[i],"VS",second[i],".csv")
    
            if (length(FDRl)>0){
                mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
                genes <- siggenes.table$genes.lo[,3]
                G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
                FinalDEGsL<-cbind(genes,siggenes.table$genes.lo)
                FinalDEGsL<-merge(FinalDEGsL,G_list,by.x="genes",by.y="ensembl_gene_id")
                FinalDEGsL[,3]<-FinalDEGsL[,10]
                FinalDEGsL<-FinalDEGsL[,c(-1,-10)]
                FinalDEGsL<-FinalDEGsL[order(FinalDEGsL[,8]),]
                cat(paste0("Low-regulated genes in the ",second[i]," in ",first[i]," VS ",second[i],"\n"))
        
                write.csv(FinalDEGsL, file = paste0("Low-regulated-",name,second[i],"in",first[i],"VS",second[i],".csv"))
                write.csv(FinalDEGsL, file = paste0("Up-regulated-",name,first[i],"in",first[i],"VS",second[i],".csv"))
                DEGsS<-c(DEGsS,FinalDEGsL[,2])
                DEGsE<-c(DEGsE,as.character(FinalDEGsL[,3]))
            }
    
            if (length(FDRu)>0){
                mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
                genes <- siggenes.table$genes.up[,3]
                G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
                FinalDEGsU<-cbind(genes,siggenes.table$genes.up)
                FinalDEGsU<-merge(FinalDEGsU,G_list,by.x="genes",by.y="ensembl_gene_id")
                FinalDEGsU[,3]<-FinalDEGsU[,10]
                FinalDEGsU<-FinalDEGsU[,c(-1,-10)]
                FinalDEGsU<-FinalDEGsU[order(FinalDEGsU[,8]),]
                cat(paste0("Up-regulated genes in the ",second[i]," in ",first[i]," VS ",second[i],"\n"))
        
                write.csv(FinalDEGsU, file = paste0("Up-regulated-",name,second[i],"in",first[i],"VS",second[i],".csv"))
                write.csv(FinalDEGsU, file = paste0("Low-regulated-",name,first[i],"in",first[i],"VS",second[i],".csv"))
                DEGsS<-c(DEGsS,FinalDEGsU[,2])
                DEGsE<-c(DEGsE,as.character(FinalDEGsU[,3]))
            }

        }else{
            DEGsTable[i,1]<-paste0(first[i]," VS ",second[i])
            DEGsTable[i,2]<-second[i]
            DEGsTable[i,3]<-NA
            DEGsTable[i,4]<-NA
            DEGsTable[i,5]<-NA
            DEGsTable[i,6]<-NA

            DEGsTable[i+com,1]<-paste0(first[i]," VS ",second[i])
            DEGsTable[i+com,2]<-first[i]
            DEGsTable[i+com,3]<-NA
            DEGsTable[i+com,4]<-NA
            DEGsTable[i+com,5]<-NA
            DEGsTable[i+com,6]<-NA
        }
    }
    cat("The results of DEGs are saved in your directory","\n")
    colnames(DEGsTable)<-c("Comparisons","Target cluster","Up-regulated genes","File name","Low-regulated genes","File name")
    write.csv(DEGsTable, file = "DEGsTable.csv")
    print(DEGsTable)
    sigDEG<-cbind(DEGsE,DEGsS)
    write.csv(sigDEG, file = "sigDEG.csv")
    return(list(sigDEG,DEGsTable))
}









DEGanalysisM<- function(object,Cluster_ID,fdr,name,First,Second){
    Nam <-colnames(object)
    num<-c(1:K)
    num1<- paste("CL", num, sep="")
    for (n in num){
        Nam<-ifelse((Cluster_ID==n),num1[n],Nam)
    }
    colnames(object)<-Nam
    sg1 <- object[,which(colnames(object)==First )]
    sg2 <- object[,which(colnames(object)==Second)]
    sg<-cbind(sg1,sg2)
    
    sg3 <- factor(gsub(paste0("(",First,"|",Second,").*"), "\\1", colnames(sg)), levels = c(paste0(First), paste0(Second)))
    sg3<-sg3[!is.na(sg3)]
    
    colnames(sg)<-sg3
    len<-c(length(sg[,which(colnames(sg)==First)]),length(sg[,which(colnames(sg)==Second)]))
    print(len)
        y <- c(rep(1:2,len))
        L<-as.matrix(sg)
        gname<-rownames(sg)
        x<-L
        data=list(x=x,y=y, geneid=gname)
        samr.obj<-samr(data, resp.type="Two class unpaired", assay.type="seq",nperms=100,nresamp=20,testStatistic="wilcoxon",random.seed=15)
        delta.table <- samr.compute.delta.table(samr.obj)
        DEGsTable<-data.frame()
        DEGsE<-c()
        DEGsS<-c()
        wm<-which.min(delta.table[,5])
        if (delta.table[wm,5]<=fdr){
            w<-which(delta.table[,5]<= fdr)
            delta<-delta.table[w[1],1]-0.001
    
            samr.plot(samr.obj, delta)
            title(paste0("DEGs in the ",Second," in ",First," VS ",Second))
    
            siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, data, delta.table)
        
            FDRl<-as.numeric(siggenes.table$genes.lo[,8])/100
            FDRu<-as.numeric(siggenes.table$genes.up[,8])/100

            siggenes.table$genes.lo[,8]<- FDRl
            siggenes.table$genes.up[,8]<- FDRu
    
            DEGsTable[1,1]<-paste0(First," VS ",Second)
            DEGsTable[1,2]<-Second
            DEGsTable[1,3]<-length(FDRu)
            DEGsTable[1,4]<-paste0("Up-regulated-",name,Second,"in",First,"VS",Second,".csv")
            DEGsTable[1,5]<-length(FDRl)
            DEGsTable[1,6]<-paste0("Low-regulated-",name,Second,"in",First,"VS",Second,".csv")

            DEGsTable[2,1]<-paste0(First," VS ",Second)
            DEGsTable[2,2]<-First
            DEGsTable[2,3]<-length(FDRu)
            DEGsTable[2,4]<-paste0("Low-regulated-",name,First,"in",First,"VS",Second,".csv")
            DEGsTable[2,5]<-length(FDRl)
            DEGsTable[2,6]<-paste0("Up-regulated-",name,First,"in",First,"VS",Second,".csv")
    
            if (length(FDRl)>0){
                mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
                genes <- siggenes.table$genes.lo[,3]
                G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
                FinalDEGsL<-cbind(genes,siggenes.table$genes.lo)
                FinalDEGsL<-merge(FinalDEGsL,G_list,by.x="genes",by.y="ensembl_gene_id")
                FinalDEGsL[,3]<-FinalDEGsL[,10]
                FinalDEGsL<-FinalDEGsL[,c(-1,-10)]
                FinalDEGsL<-FinalDEGsL[order(FinalDEGsL[,8]),]
                cat(paste0("Low-regulated genes in the ",Second," in ",First," VS ",Second,"\n"))
        
                write.csv(FinalDEGsL, file = paste0("Low-regulated-",name,Second,"in",First,"VS",Second,".csv"))
                write.csv(FinalDEGsL, file = paste0("Up-regulated-",name,First,"in",First,"VS",Second,".csv"))
                DEGsS<-c(DEGsS,FinalDEGsL[,2])
                DEGsE<-c(DEGsE,as.character(FinalDEGsL[,3]))
            }
    
            if (length(FDRu)>0){
                mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
                genes <- siggenes.table$genes.up[,3]
                G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
                FinalDEGsU<-cbind(genes,siggenes.table$genes.up)
                FinalDEGsU<-merge(FinalDEGsU,G_list,by.x="genes",by.y="ensembl_gene_id")
                FinalDEGsU[,3]<-FinalDEGsU[,10]
                FinalDEGsU<-FinalDEGsU[,c(-1,-10)]
                FinalDEGsU<-FinalDEGsU[order(FinalDEGsU[,8]),]
                cat(paste0("Up-regulated genes in the ",Second," in ",First," VS ",Second,"\n"))
        
                write.csv(FinalDEGsU, file = paste0("Up-regulated-",name,Second,"in",First,"VS",Second,".csv"))
                write.csv(FinalDEGsU, file = paste0("Low-regulated-",name,First,"in",First,"VS",Second,".csv"))
                DEGsS<-c(DEGsS,FinalDEGsU[,2])
                DEGsE<-c(DEGsE,as.character(FinalDEGsU[,3]))
            }

        }else{
            DEGsTable[1,1]<-paste0(First," VS ",Second)
            DEGsTable[1,2]<-Second
            DEGsTable[1,3]<-NA
            DEGsTable[1,4]<-NA
            DEGsTable[1,5]<-NA
            DEGsTable[1,6]<-NA

            DEGsTable[2,1]<-paste0(First," VS ",Second)
            DEGsTable[2,2]<-First
            DEGsTable[2,3]<-NA
            DEGsTable[2,4]<-NA
            DEGsTable[2,5]<-NA
            DEGsTable[2,6]<-NA
        }
    
    cat("The results of DEGs are saved in your directory","\n")
    colnames(DEGsTable)<-c("Comparisons","Target cluster","Up-regulated genes","File name","Low-regulated genes","File name")
    write.csv(DEGsTable, file = "DEGsTable.csv")
    print(DEGsTable)
    sigDEG<-cbind(DEGsE,DEGsS)
    write.csv(sigDEG, file = "sigDEG.csv")
    return(list(sigDEG,DEGsTable))
}



setGeneric("FindOutliersKM", function(object,outminc=5,outlg=2,probthr=1e-3,thr=2**-(1:40),outdistquant=.75) standardGeneric("FindOutliersKM"))

setMethod("FindOutliersKM",
          signature = "PSCANseq ",
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
            axis(1,at=(y[1:(length(y)-1)] + y[-1])/2,lab=1:max(p))
            box()
            cat("The following cells are considered as outlier cells:",which(object@cpart>K),"\n")
            print(which(object@cpart>K))
            LL= which(object@cpart>K)  
            return(LL)
          }
        )








setGeneric("FindOutliersMB", function(object,outminc=5,outlg=2,probthr=1e-3,thr=2**-(1:40),outdistquant=.75) standardGeneric("FindOutliersMB"))

setMethod("FindOutliersMB",
          signature = "PSCANseq ",
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
            axis(1,at=(y[1:(length(y)-1)] + y[-1])/2,lab=1:max(p))
            box()
            cat("The following cells are considered as outlier cells:",which(object@cpart>K),"\n")
            print(which(object@cpart>K))
            LL= which(object@cpart>K)  
            return(LL)
          }
        )

binompval <- function(p,N,n){
  pval   <- pbinom(n,round(N,0),p,lower.tail=TRUE)
  pval[!is.na(pval) & pval > 0.5] <- 1-pval[!is.na(pval) & pval > 0.5]
  return(pval)
}



setGeneric("MBClustDiffGenes", function(object,fdr=.01) standardGeneric("MBClustDiffGenes"))
setMethod("MBClustDiffGenes",
          signature = "PSCANseq ",
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


setGeneric("KMClustDiffGenes", function(object,fdr=.01) standardGeneric("KMClustDiffGenes"))
setMethod("KMClustDiffGenes",
          signature = "PSCANseq ",
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









estimateSizeFactorsForMatrix = function (counts, locfunc = median){
	loggeomeans <- rowMeans(log(counts))
	apply(counts, 2, function(cnts) exp(locfunc((log(cnts) - loggeomeans)[is.finite(loggeomeans)])))
}

