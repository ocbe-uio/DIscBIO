#' @title Volcano Plot
#' @description Plotting differentially expressed genes (DEGs) in a particular
#'   cluster. Volcano plots are used to readily show the DEGs by plotting
#'   significance versus fold-change on the y and x axes, respectively.
#' @param object A data frame showing the differentially expressed genes (DEGs)
#'   in a particular cluster
#' @param value A numeric value of the false discovery rate. Default is 0.05..
#'   Default is 0.05
#' @param fc A numeric value of the fold change. Default is 0.5.
#' @param FS A numeric value of the font size. Default is 0.4.
#' @param name A string vector showing the name to be used on the plot title
#' @importFrom graphics title
#' @importFrom utils write.csv
#' @return A volcano plot
#' @export
VolcanoPlot <- function(object, value = 0.05, name = NULL, fc = 0.5, FS = .4) {
    if (length(object[1, ]) > 8) {
        object <- object[, -1]
    }
    NO0 <- object[, 8]
    NO0 <- NO0[which(NO0 != 0)]
    w <- which.min(NO0)
    adjV <- NO0[w] / 100
    object[, 8] <- ifelse(object[, 8] == 0 & length(adjV) > 0, adjV, object[, 8])
    if (all(object[, 8] == 0)) {
        message("All q-values are 0. Adjusting")
        object[, 8] <- object[, 8] + 1e-10
    }
    with(
        object,
        plot(
            abs(object[, 7]),
            -log10(object[, 8]),
            pch = 20,
            cex = 2,
            las = 1,
            xlab = "log2 Fold Change",
            ylab = "-log10 FDR",
            sub = paste("Volcano plot", name),
            font.sub = 4,
            col.sub = "black"
        )
    )
    FC <- subset(object, abs(object[, 7]) > fc) # Fold Change
    sigFC <- subset(
        object, object[, 8] < value & abs(object[, 7]) > fc
    ) # Significant genes
    with(FC, points(
        abs(FC[, 7]),
        -log10(FC[, 8]),
        pch = 20,
        cex = 2,
        col = "red"
    ))
    with(sigFC, points(
        abs(sigFC[, 7]),
        -log10(sigFC[, 8]),
        pch = 20,
        cex = 2,
        col = "blue"
    ))
    add_legend(
        "topleft",
        legend = c(
            paste0("DEGs (FC < ", fc, " - FDR> ", value, ")   "),
            paste0("DEGs (FC > ", fc, " - FDR> ", value, ")"),
            paste0("DEGs (FC > ", fc, " - FDR< ", value, ")   ")
        ),
        pch = 20,
        col = c("black", "red", "blue"),
        horiz = TRUE,
        bty = 'n',
        cex = FS
    )
}
