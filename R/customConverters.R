
#' Automatic Feature Id Conversion.
#'
#' Attempt to automatically convert non-ENSEMBL feature identifiers to
#' ENSEMBL identifiers. Features are included as rownames of the input
#' data.frame (or matrix). It is assumed that feature identifiers
#' (i.e., rownames of x) come from human or mouse genomes, and
#' are either OFFICIAL SYMBOLS or ENTREZIDS. If less than 20% of match
#' is identified, an error will be thrown.
#'
#' @param x data.frame or matrix including raw counts (typically, UMIs),
#'   where rows are features (genes) and rownames are feature identifiers
#'   (SYMBOLs or ENTREZIDs).
#' @param verbose logical, shall messages be printed to inform about
#'   conversion advances.
#'
#' @return a data.frame where rownames are ENSEMBL IDs. The new feature IDs are
#' automatically imputed based on existing feature IDs (SYMBOLs or ENTREZIDs).
#'
#' @export
customConvertFeats <- function(x, verbose = TRUE) {
    if (! (is.matrix(x) || is.data.frame(x)) ){
        stop("wrong input")
    }

    # uncomment
    myDict <- DIscBIO::HumanMouseGeneIds
    myDict <- myDict[!is.na(myDict[, "ENSEMBL"]), ]

    #
    xx <- rownames(x)
    keep.ei <- xx %in% myDict[, "ENTREZID"]
    keep.sy <- xx %in% myDict[, "SYMBOL"]
    keep.ez <- xx %in% myDict[, "ENSEMBL"]

    rat.ez <- sum(keep.ez) / length(xx)
    rat.ei <- sum(keep.ei) / length(xx)
    rat.sy <- sum(keep.sy) / length(xx)

    # Automatic selection and replacement
    if (rat.ez >= rat.sy && rat.ez >= rat.ei && rat.ez >= 0.2) {
        # nothing to do, return
        return(x)
    } else if (rat.ei >= rat.sy && rat.ei >= 0.2) {
        x <- x[!duplicated(xx), ]
        x <- x[rownames(x) %in% myDict[, "ENTREZID"], ]

        tmpDict <- myDict$ENSEMBL
        names(tmpDict) <- myDict$ENTREZID
        rownames(x) <- as.character(tmpDict[rownames(x)])
    } else if (rat.ei < rat.sy && rat.sy >= 0.2) {
        x <- x[!duplicated(xx), ]
        x <- x[rownames(x) %in% myDict[, "SYMBOL"], ]

        tmpDict <- myDict$ENSEMBL
        names(tmpDict) <- myDict$SYMBOL
        rownames(x) <- as.character(tmpDict[rownames(x)])
    } else {
        message("Feat ID Conversion could not be performed")
    }
    return(x)
}



#' Convert Single Cell Data Objects to DISCBIO.
#'
#' Initialize a DISCBIO-class object starting from a
#' SingleCellExperiment object or a Seurat object
#'
#' @param x an object of class Seurat or SingleCellExperiment
#' @param ... additional parameters to pass to the function
#'
#' @details Additional parameters to pass to `list` include, if x is
#'   a Seurat object, "assay", which is a string indicating the assay slot
#'   used to obtain data from (defaults to 'RNA')
#'
#' @return a DISCBIO-class object
#'
#' @export
#'
as.DISCBIO <- function(x, ...) {
    if ("Seurat" %in% class(x)) {
        # Get Arguments and parse out what we want
        all.args <- list(...)

        # Fetch arguments we care about
        if ("assay" %in% all.args) {
            assay <- all.args[["assay"]]
        } else {
            assay <- "RNA"
        }

        # Get feats and sample names
        all.feats <- base::as.character(
            rownames(x@assays[[assay]]@meta.features)
        )
        all.cells <- base::as.character(rownames(x@meta.data))

        # Get data
        all.counts <- base::data.frame(
            base::as.matrix(x@assays[[assay]]@counts)
        )

        # re-write row and colnames
        rownames(all.counts) <- all.feats
        colnames(all.counts) <- all.cells

        # prep output and return
        if (sum(grepl("^ENS", all.feats)) / length(all.feats) < 0.5) {
            all.counts <- customConvertFeats(all.counts)
        }
        y <- DISCBIO(all.counts)
    } else if ("SingleCellExperiment" %in% class(x)) {
        y <- DISCBIO(x)
    } else {
        stop("Conversion not supported")
    }
    return(y)
}


#' Convert a DISCBIO object to a SingleCellExperiment.
#'
#' Extract the SingleCellExperiment input data from the corresponding input slot
#' in a DISCBIO-class object
#'
#' @param x an object of class DISCBIO
#'
#' @return a SingleCellExperiment-class object
#'
#' @export
#' @examples
#' g1_disc <- DISCBIO(valuesG1msTest)
#' class(g1_disc)
#' g1_sce <- DISCBIO2SingleCellExperiment(g1_disc)
#' class(g1_sce)
#'
#'
DISCBIO2SingleCellExperiment <- function(x) {
    return(x@SingleCellExperiment)
}
