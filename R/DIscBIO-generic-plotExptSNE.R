#' @title Highlighting gene expression in the t-SNE map
#' @description The t-SNE map representation can also be used to analyze
#'   expression of a gene or a group of genes, to investigate cluster specific
#'   gene expression patterns
#' @param object \code{DISCBIO} class object.
#' @param g  Individual gene name or vector with a group of gene names
#'   corresponding to a subset of valid row names of the \code{ndata} slot of
#'   the \code{DISCBIO} object.
#' @param n String of characters representing the title of the plot. Default is
#'   NULL and the first element of \code{g} is chosen.
#' @return t-SNE plot for one particular gene

setGeneric("plotExptSNE", function(object, g, n = NULL)
    standardGeneric("plotExptSNE")
)

#' @export
#' @rdname plotExptSNE
setMethod(
    "plotExptSNE",
    signature = "DISCBIO",
    definition = function(object, g, n = NULL) {
        # ======================================================================
        # Validation
        # ======================================================================
        ran_k <- length(object@tsne) > 0
        ran_m <- length(object@MBtsne) > 0
        if (ran_k) {
            x <- object@tsne
        } else if (ran_m) {
            x <- object@MBtsne
        } else {
            stop("run comptSNE before plotExptSNE")
        }
        if (length(intersect(g, rownames(object@ndata))) < length(unique(g)))
            stop(
                "second argument does not correspond to set of rownames slot",
                "ndata of SCseq object"
            )
        if (is.null(n)) n <- g[1]
        # ======================================================================
        # Plotting
        # ======================================================================
        logObj <- log(object@ndata)
        l <- apply(logObj[g, ] - .1, 2, sum) + .1
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
        opar <- par(mar = c(3, 5, 2.5, 2))
        on.exit(par(opar))
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
        for (k in 1:length(v)) {
            points(
                x[k, 1],
                x[k, 2],
                col = ColorRamp[v[k]],
                pch = 20,
                cex = 1.5
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
            las = 1,
            xaxt = "n"
        )
        layout(1)
    }
)
