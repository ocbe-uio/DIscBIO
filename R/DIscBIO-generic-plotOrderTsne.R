#' @title Plotting the pseudo-time ordering in the t-SNE map
#' @description The tSNE representation can also be used to show the pseudo-time
#'   ordering.
#' @param object \code{DISCBIO} class object.
#' @return A plot of the pseudo-time ordering.
setGeneric("plotOrderTsne", function(object) {
  standardGeneric("plotOrderTsne")
})

#' @export
#' @rdname plotOrderTsne
setMethod(
  "plotOrderTsne",
  signature = "DISCBIO",
  definition = function(object) {
    ran_k <- length(object@tsne) > 0
    ran_m <- length(object@MBtsne) > 0
    if (ran_k) {
      total <- rbind(object@ndata, object@kordering)
      clustering_method <- "k-means"
      x <- object@tsne
    } else if (ran_m) {
      total <- rbind(object@ndata, object@MBordering)
      clustering_method <- "model-based"
      x <- object@MBtsne
    } else {
      stop("run comptsne before plotOrderTsne")
    }

    rownames(total)[nrow(total)] <- paste(
      "Pseudo-time ordering of", clustering_method, "clustering"
    )
    g <- rownames(total)[nrow(total)]
    n <- g[1]
    l <- apply(total[g, ] - .1, 2, sum) + .1

    mi <- min(l, na.rm = TRUE)
    ma <- max(l, na.rm = TRUE)
    ColorRamp <- colorRampPalette(
      rev(brewer.pal(n = 7, name = "RdYlBu"))
    )(100)
    ColorLevels <- seq(mi, ma, length = length(ColorRamp))
    v <- round((l - mi) / (ma - mi) * 99 + 1, 0)
    layout(
      matrix(
        data = c(1, 3, 2, 4),
        nrow = 2,
        ncol = 2
      ),
      widths = c(5, 1, 5, 1),
      heights = c(5, 1, 1, 1)
    )
    opar <- withr::local_par(mar = c(3, 5, 2.5, 2))
    on.exit(withr::local_par(opar))
    plot(
      x,
      xlab = "Dim 1",
      ylab = "Dim 2",
      main = n,
      pch = 20,
      cex = 0,
      col = "grey",
      las = 1
    )
    for (k in seq_len(length(v))) {
      points(
        x[k, 1],
        x[k, 2],
        col = ColorRamp[v[k]],
        pch = 20,
        cex = 1.5
      )
    }
    opar <- withr::local_par(mar = c(3, 2.5, 2.5, 2))
    on.exit(withr::local_par(opar))
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
      las = 1,
      xaxt = "n"
    )
    layout(1)
  }
)
