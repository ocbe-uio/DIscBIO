#' @title K-means runs
#' @param data data
#' @param krange krange
#' @param criterion criterion
#' @param iter.max iter.max
#' @param runs number of runs
#' @param scaledata scaledata
#' @param alpha alpha
#' @param critout critout
#' @param plot plot
#' @param method method
#' @param ... extra parameters to be passed to kmm
#' @importFrom amap Kmeans
#' @importFrom graphics pairs
#' @importFrom fpc cluster.stats calinhara dudahart2
Kmeansruns <- function (data, krange = 2:10, criterion = "ch", iter.max = 100, 
    runs = 100, scaledata = FALSE, alpha = 0.001, critout = FALSE, 
    plot = FALSE, method="euclidean", ...) 
{
    data <- as.matrix(data)
    if (criterion == "asw") 
        sdata <- dist(data)
    if (scaledata) 
        data <- scale(data)
    cluster1 <- 1 %in% krange
    crit <- numeric(max(krange))
    km <- list()
    for (k in krange) {
        if (k > 1) {
            minSS <- Inf
            kmopt <- NULL
            for (i in 1:runs) {
                options(show.error.messages = FALSE)
                repeat {
                  kmm <- try(Kmeans(data, k, iter.max = iter.max, method=method,
                    ...))
                  if (class(kmm) != "try-error") 
                    break
                }
                options(show.error.messages = TRUE)
                swss <- sum(kmm$withinss)
                if (swss < minSS) {
                  kmopt <- kmm
                  minSS <- swss
                }
                if (plot) {
                  par(ask = TRUE)
                  pairs(data, col = kmm$cluster, main = swss)
                }
            }
            km[[k]] <- kmopt
            crit[k] <- switch(criterion, asw = cluster.stats(sdata, 
                km[[k]]$cluster)$avg.silwidth, ch = calinhara(data, 
                km[[k]]$cluster))
            if (critout) 
                cat(k, " clusters ", crit[k], "\n")
        }
    }
    if (cluster1) 
        cluster1 <- dudahart2(data, km[[2]]$cluster, alpha = alpha)$cluster1
    k.best <- which.max(crit)
    if (cluster1) 
        k.best <- 1
    km[[k.best]]$crit <- crit
    km[[k.best]]$bestk <- k.best
    out <- km[[k.best]]
    out
}