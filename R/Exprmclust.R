#' @title title
#' @description description 
#' @export
#' @param data data
#' @param clusternum clusternum
#' @param modelNames modelNames
#' @param reduce reduce
#' @param cluster cluster
#' @importFrom mclust Mclust
#' @importFrom stats dist prcomp lm
#' @importFrom igraph graph.adjacency minimum.spanning.tree
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