#' @title Jaccard’s similarity
#' @description Robustness of the clusters can be assessed by Jaccard’s similarity, which reflects the reproducibility of individual clusters across bootstrapping runs.
#' Jaccard’s similarity is the intersect of two clusters divided by the union. 
#' @export
#' @param object \code{DISCBIO} class object.
#' @param Clustering Clustering has to be one of the following: ["K-means","MB"]. Default is "K-means"
#' @param K A numeric value of the number of clusters
#' @param plot if `TRUE`, plots the mean Jaccard similarities
#' @importFrom philentropy distance
#' @importFrom boot boot
#' @importFrom graphics barplot box
#' @return A plot of the mean Jaccard similarity coefficient per cluster.
#' @examples
#' sc <- DISCBIO(valuesG1msReduced)
#' sc <- Clustexp(sc, cln=3, quiet=TRUE) # K-means clustering
#' Jaccard(sc, Clustering="K-means", K=3)
Jaccard <- function(object, Clustering="K-means", K, plot = TRUE) {
    JACCARD <-c()
    
    # Validation
    if (length(object@kmeans$kpart) == 0) {
        stop("run Clustexp before Jaccard")
    }
    if (!(Clustering %in% c( "K-means","MB"))) {
        stop("Clustering has to be either K-means or MB")
    }

    JS <- function(data, indices) {
      d <- data[indices, ] # allows boot to select sample 
      jac <- suppressMessages(distance(t(d), method = "jaccard"))
      jac1 <- 1 - jac
      JSmean <- mean(jac1)
      return(JSmean)
    }
    for (i in 1:K) {
        # Optimize by avoiding if every loop. Only thing variable is data
        if (Clustering == "K-means") {
            results <- boot(
                data = object@fdata[, which(object@kmeans$kpart == i)],
                statistic = JS,
                R = 100,
                stype = "f"
            )
            # to get the mean of all bootstrappings (mean of mean Jaccard values)
            JACCARD[i] <- round(mean(results$t), digits = 3)
        }
        if (Clustering == "MB") {
            results <- boot(
                data = object@fdata[, which(object@MBclusters$clusterid == i)],
                statistic = JS,
                R = 100,
                stype = "f"
            )
            # to get the mean of all bootstrappings (mean of mean Jaccard values)
            JACCARD[i] <- round(mean(results$t), digits = 3)
        }
    }
    if (plot) {
        barplot(
            height = JACCARD,
            names.arg = 1:length(JACCARD),
            ylab = "Mean Jaccard's similarity values",
            xlab = "Clusters",
            las = 1,
            ylim = c(0,1),
            col = c("black", "blue", "green", "red", "yellow", "gray")
        )
        box()
    }
    return(JACCARD)
}