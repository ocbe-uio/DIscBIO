#' @title Jaccard’s similarity
#' @description Robustness of the clusters can be assessed by Jaccard’s
#'   similarity, which reflects the reproducibility of individual clusters
#'   across bootstrapping runs. Jaccard’s similarity is the intersect of two
#'   clusters divided by the union.
#' @export
#' @param object \code{DISCBIO} class object.
#' @param Clustering Clustering has to be one of the following:
#'   ["K-means","MB"]. Default is "K-means"
#' @param K A numeric value of the number of clusters
#' @param plot if `TRUE`, plots the mean Jaccard similarities
#' @param R number of bootstrap replicates
#' @importFrom philentropy distance
#' @importFrom graphics barplot box
#' @return A plot of the mean Jaccard similarity coefficient per cluster.
Jaccard <- function(object, Clustering = "K-means", K, plot = TRUE, R = 100) {
  JACCARD <- vector()

  # Validation
  if (!(Clustering %in% c("K-means", "MB"))) {
    stop("Clustering has to be either K-means or MB")
  }
  for (i in 1:K) {
    # Optimize by avoiding if every loop. Only thing variable is data
    if (Clustering == "K-means") {
      target_col <- object@kmeans$kpart
    } else if (Clustering == "MB") {
      target_col <- object@MBclusters$clusterid
    }
    results <- bootstrap(object@fdata[, which(target_col == i)], R)
    # to get the mean of all bootstrappings (mean of mean Jaccard values)
    JACCARD[i] <- round(mean(results), digits = 3)
  }
  if (plot) {
    barplot(
      height    = JACCARD,
      names.arg = seq_len(length(JACCARD)),
      ylab      = "Mean Jaccard's similarity values",
      xlab      = "Clusters",
      las       = 1,
      ylim      = c(0, 1),
      col       = c("black", "blue", "green", "red", "yellow", "gray")
    )
    box()
  }
  return(JACCARD)
}
