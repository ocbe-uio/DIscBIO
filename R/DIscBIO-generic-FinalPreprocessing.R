#' @title Final Preprocessing
#' @description This function generates the final filtered normalized dataset.
#' @param object \code{DISCBIO} class object.
#' @param GeneFlitering GeneFlitering has to be one of the followings:
#'   ["NoiseF","ExpF"]. Default is "NoiseF"
#' @param export A logical vector that allows writing the final gene list in
#'   excel file. Default is TRUE.
#' @param quiet if `TRUE`, intermediary output is suppressed
#' @param fileName File name for exporting (if `export = TRUE`)
#' @return The DISCBIO-class object input with the FinalGeneList slot filled.
#' @examples
#' sc <- DISCBIO(valuesG1msTest)
#' sc <- NoiseFiltering(sc, percentile=0.9, CV=0.2, export=FALSE)
#' sc <- FinalPreprocessing(sc, GeneFlitering="NoiseF", export=FALSE)
#'
setGeneric(
    "FinalPreprocessing",
    function(
        object,
        GeneFlitering = "NoiseF",
        export = FALSE,
        quiet = FALSE,
        fileName = "filteredDataset"
    )
    standardGeneric("FinalPreprocessing")
)

#' @export
#' @rdname FinalPreprocessing
setMethod(
    "FinalPreprocessing",
    signature = "DISCBIO",
    definition = function(object, GeneFlitering, export, quiet, fileName)
    {
        if (GeneFlitering == "NoiseF") {
            if (length(object@noiseF) < 1)
                stop("run NoiseFiltering before running FinalPreprocessing")
            if (nrow(object@fdata) < 1)
                stop("run Normalizedata before running FinalPreprocessing")
            gene_list <- object@noiseF
            gene_names <- rownames(object@fdata)
            idx_genes <- is.element(gene_names, gene_list)
            gene_names2 <- gene_names[idx_genes]
            filteredDataset <- object@fdata[gene_names2, ]
            object@fdata <- filteredDataset
            object@FinalGeneList <- rownames(filteredDataset)
        }
        if (GeneFlitering == "ExpF") {
            if (nrow(object@fdata) < 1)
                stop("run Normalizedata before running FinalPreprocessing")
            filteredDataset <- object@fdata
            object@FinalGeneList <- rownames(filteredDataset)
        }
        if (!quiet) {
            message(
                "The gene filtering method = Noise filtering\n\n",
                "The Filtered Normalized dataset contains:\n",
                "Genes: ", length(filteredDataset[, 1]), "\n",
                "cells: ", length(filteredDataset[1, ]), "\n\n"
            )
        }
        if (export) {
            fileNameExt <- paste0(fileName, ".Rdata")
            if (!quiet) {
                message(
                    "The Filtered Normalized dataset was saved as: ",
                    fileNameExt
                )
            }
            save(filteredDataset, file = fileNameExt)
        }
        return(object)
    }
)
