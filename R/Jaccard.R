#' @title title
#' @description description 
#' @export
#' @param object object
#' @param Clustering Clustering
#' @param K K
#' @param plot if `TRUE`, plots the mean Jaccard similarities
#' @importFrom philentropy distance
#' @importFrom boot boot
#' @importFrom graphics barplot box
Jaccard <- function(object, Clustering, K, plot = TRUE) {
    JACCARD <-c()
    if ( ! ( Clustering %in% c( "K-means","MB") ) ) stop("Clustering has to be either K-means or MB")
    JS <- function(data, indices) {
      d <- data[indices, ] # allows boot to select sample 
      jac <- suppressMessages(distance(t(d), method = "jaccard"))
      jac1 <- 1 - jac
      JSmean <- mean(jac1)
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
    if (plot) {
        barplot(JACCARD,names.arg=1:length(JACCARD),ylab="Mean Jaccard's similarity values",xlab="Clusters",
        las=1,ylim=c(0,1),col=c("black","blue","green","red","yellow","gray"))
        box()
    }
    return(JACCARD)
}