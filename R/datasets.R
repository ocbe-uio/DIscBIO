#' @title Single-cells data from a myxoid liposarcoma cell line
#'
#' @description A dataset of 30 single cells from a myxoid liposarcoma cell
#'   line. Columns refer to samples and rows refer to genes. The last 92 rows
#'   refer to external RNA controls consortium (ERCC) spike-ins. This dataset is
#'   part of a larger dataset containing 94 single cells. The complete dataset
#'   is fully compatible with this package and an rda file can be obtained at
#'   https://github.com/ocbe-uio/DIscBIO/blob/dev/data/valuesG1ms.rda
#' @name valuesG1msRed
#' @docType data
NULL

#' Human and Mouse Gene Identifiers.
#'
#' Data.frame including ENTREZID, SYMBOL, and ENSEMBL gene identifiers
#' of human and mouse genes.
#'
#' @source Data were imported, modified, and formatted from the
#' Mus.musculus (ver 1.3.1) and the Homo.sapiens (ver 1.3.1)
#' BioConductor libraries.
#'
#' @examples
#' data(HumanMouseGeneIds)
#' print(HumanMouseGeneIds[1:6,])
#' @name HumanMouseGeneIds
#' @docType data
NULL

#' @title Single-cells data from a myxoid liposarcoma cell line (K variant)
#'
#' @description Equivalent to valuesG1msReduced after treatment and K-means clustering.
#' @name valuesG1msReduced_treated_K
#' @docType data
NULL

#' @title Single-cells data from a myxoid liposarcoma cell line (MB variant)
#'
#' @description Equivalent to valuesG1msReduced after treatment and model-based clustering.
#' @name valuesG1msReduced_treated_MB
#' @docType data
NULL