#' @title Inference of outlier cells in K-means clustering
#' @description This functions performs the outlier identification
#' @param object \code{PSCANseq} class object.
#' @param outminc minimal transcript count of a gene in a clusters to be tested for being an outlier gene. Default is 5.
#' @param outlg Minimum number of outlier genes required for being an outlier cell. Default is 2.
#' @param probthr outlier probability threshold for a minimum of \code{outlg} genes to be an outlier cell. This probability is computed from a negative binomial 
#' background model of expression in a cluster. Default is 0.001.
#' @param thr probability values for which the number of outliers is computed in order to plot the dependence of the number of outliers on the probability threshold. Default is 2**-(1:40).
#' @param outdistquant Real number between zero and one. Outlier cells are merged to outlier clusters if their distance smaller than the outdistquant-quantile of
#' the distance distribution of  pairs of cells in the orginal clusters after outlier removal. Default is 0.75.
#' @param K Number of clusters to be used.
#' @param plot if `TRUE`, produces a plot of -log10prob per K
#' @param quiet if `TRUE`, intermediary output is suppressed
#' @importFrom stats coef pnbinom
#' @importFrom amap K
#' @rdname FindOutliersKM
#' @export
setMethod("FindOutliersKM",
          signature = "PSCANseq",
          definition = function(object, K, outminc,outlg,probthr,thr,outdistquant, plot = TRUE, quiet = FALSE) {
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
          if (plot) {
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
          }
          if (!quiet) {
            cat("The following cells are considered as outlier cells:",which(object@cpart>K),"\n")
            print(which(object@cpart>K))
          }
            LL= which(object@cpart>K)  
            return(LL)
          }
        )