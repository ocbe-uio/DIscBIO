#' @title Plotting clusters in a heatmap representation of the cell distances
#' @description  This functions plots a heatmap of the distance matrix grouped
#'   by clusters. Individual clusters are highlighted with rainbow colors along
#'   the x and y-axes.
#' @param object \code{DISCBIO} class object.
#' @param clustering_method either "k-means" or "model-based" ("k" and "mb" are also accepted)
#' @param hmethod  Agglomeration method used for determining the cluster order
#'   from hierarchical clustering of the cluster medoids. This should be one of
#'   "ward.D", "ward.D2", "single", "complete", "average". Default is "single".
#' @param quiet if `TRUE`, intermediary output is suppressed
#' @param rseed Random integer to fix random results.
#' @param plot if `TRUE`, plots the heatmap; otherwise, just prints cclmo
#' @return Unless otherwise specified, a heatmap and a vector of the underlying
#'   cluster order.
#' @importFrom stats hclust as.dist cor
setGeneric(
    "clustheatmap",
    function(
        object,
        clustering_method = "k-means",
        hmethod="single",
        rseed=NULL,
        quiet=FALSE,
        plot=TRUE)
	{
	    standardGeneric("clustheatmap")
	}
)

#' @export
#' @rdname clustheatmap
setMethod(
    "clustheatmap",
    signature = "DISCBIO",
    definition = function(
        object, clustering_method, hmethod, rseed, quiet, plot
    )
	{
        x <- object@fdata
		if (tolower(clustering_method) %in% c("k-means", "k")) {
			part <- object@kmeans$kpart
		} else if (tolower(clustering_method) %in% c("model-based", "mb")) {
			object@clusterpar$metric <- "pearson"
			y <- clustfun(
				object@fdata,
				clustnr = 20,
				bootnr = 50,
				metric = "pearson",
				do.gap = TRUE,
				SE.method = "Tibs2001SEmax",
				SE.factor = .25,
				B.gap = 50,
				cln = 0,
				rseed = rseed,
				quiet = quiet
			)
			object@distances <- as.matrix(y$di)
    		part <- object@MBclusters$clusterid
		}
        part <- object@kmeans$kpart
        na <- c()
        j <- 0
        for (i in 1:max(part)) {
            if (sum(part == i) == 0)
                next
            j <- j + 1
            na <- append(na, i)
            d <- x[, part == i]
            if (sum(part == i) == 1)
                cent <- d
            else
                cent <- apply(d, 1, mean)
            if (j == 1)
                tmp <- data.frame(cent)
            else
                tmp <- cbind(tmp, cent)
        }
        names(tmp) <- paste("cl", na, sep = ".")
        if (max(part) > 1)
            cclmo <-
            hclust(dist.gen(as.matrix(
                dist.gen(t(tmp), method = object@clusterpar$metric)
            )), method = hmethod)$order
        else
            cclmo <- 1
        q <- part
        for (i in 1:max(part)) {
            q[part == na[cclmo[i]]] <- i
        }
        part <- q
        di <- as.data.frame(as.matrix(dist.gen(t(object@distances))))
        pto <- part[order(part, decreasing = FALSE)]
        ptn <- c()
        for (i in 1:max(pto)) {
            pt <-
                names(pto)[pto == i]
            z <-
                if (length(pt) == 1)
                    pt
            else
                pt[hclust(as.dist(t(di[pt, pt])), method = hmethod)$order]
            ptn <- append(ptn, z)
        }
        col = c("black", "blue", "green", "red", "yellow", "gray")
        mi  <- min(di, na.rm = TRUE)
        ma  <- max(di, na.rm = TRUE)
        if (plot) {
            layout(
                matrix(
                    data = c(1, 3, 2, 4),
                    nrow = 2,
                    ncol = 2
                ),
                widths = c(5, 1, 5, 1),
                heights = c(5, 1, 1, 1)
            )
            ColorRamp <- colorRampPalette(
				brewer.pal(n = 7, name = "RdYlBu")
			)(100)
            ColorLevels <- seq(mi, ma, length = length(ColorRamp))
            if (mi == ma) {
                ColorLevels <- seq(
                    0.99 * mi, 1.01 * ma, length = length(ColorRamp)
                )
            }

            opar <- par(mar = c(3, 5, 2.5, 2))
            on.exit(par(opar))
            image(as.matrix(di[ptn, ptn]), col = ColorRamp, axes = FALSE)
            abline(0, 1)
            box()

            tmp <- c()
            for (u in 1:max(part)) {
                ol <- (0:(length(part) - 1) /
                    (length(part) - 1))[ptn %in% names(x)[part == u]]
                points(
                    rep(0, length(ol)),
                    ol,
                    col = col[cclmo[u]],
                    pch = 15,
                    cex = .75
                )
                points(
                    ol,
                    rep(0, length(ol)),
                    col = col[cclmo[u]],
                    pch = 15,
                    cex = .75
                )
                tmp <- append(tmp, mean(ol))
            }
            axis(1, at = tmp, labels = cclmo)
            axis(2, at = tmp, labels = cclmo)
            opar <- par(mar = c(3, 2.5, 2.5, 2))
            on.exit(par(opar))
            image(
                1,
                ColorLevels,
                matrix(
                    data = ColorLevels,
                    ncol = length(ColorLevels),
                    nrow = 1
                ),
                col = ColorRamp,
                xlab = "",
                ylab = "",
                las = 2,
                xaxt = "n"
            )
            layout(1)
        }
        return(cclmo)
    }
)