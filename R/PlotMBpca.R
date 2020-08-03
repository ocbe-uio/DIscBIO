#' @title Plotting pseudo-time ordering or gene expression in Model-based clustering in PCA
#' @description The PCA representation can either be used to show pseudo-time ordering or the gene expression of a particular gene.
#' @param object \code{DISCBIO} class object.
#' @param type either `order` to plot pseudo-time ordering or `exp` to plot gene expression
#' @param g  Individual gene name or vector with a group of gene names
#'   corresponding to a subset of valid row names of the \code{ndata} slot of
#'   the \code{DISCBIO} object. Ignored if `type="order"`.
#' @param n String of characters representing the title of the plot. Default is
#'   NULL and the first element of \code{g} is chosen. Ignored if
#'   `type="order"`.
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics layout par image
#' @return A plot of the PCA.
#' @export

PlotMBpca <- function(object, type="order", g=NULL, n=NULL) {
    # ==========================================================================
    # Validation
    # ==========================================================================
    data <- object@MBclusters
    if (type == "exp") {
        if (is.null(g)) {
            stop('g must be provided if type="exp"')
        }
        if (length(intersect(g, rownames(object@ndata))) < length(unique(g))) {
            stop(
                "second argument does not correspond to set of rownames slot",
                "ndata of SCseq object"
            )
        }
        if (is.null(n)) {
            n <- g[1]
        }
        l <- apply(object@ndata[g,] - .1, 2, sum) + .1
        x <- data$pcareduceres
    } else if (type == "order") {
        MBordertable <- cbind(data$pcareduceres, object@MBordering)
        l <- MBordertable[, 3]
        x <- MBordertable
    } else {
        stop("Invalid type. Valid alternatives as 'order' and 'exp'.")
    }
    # ==========================================================================
    # Plotting
    # ==========================================================================
    mi <- min(l, na.rm = TRUE)
    ma <- max(l, na.rm = TRUE)
    ColorRamp <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100)
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
    opar <- par(mar = c(5, 5, 2.5, 2))
    on.exit(par(opar))
    plot(
        x[, 1],
        x[, 2],
        xlab = "PC1",
        ylab = "PC2",
        pch = 20,
        cex = 0,
        col = "grey",
        las = 1,
        main = n
    )
    for (k in 1:length(v)) {
        points(
            x[k, 1],
            x[k, 2],
            col = ColorRamp[v[k]],
            pch = 20,
            cex = 2
        )
    }
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