#' @title Defining protein-protein interactions (PPI) over a list of genes,
#' @description This function uses STRING-api. The outcome of STRING analysis
#'   will be stored in tab separated values (TSV) files.
#' @export
#' @param data A gene list.
#' @param FileName A string vector showing the name to be used to save the
#'   resulted table. If null, no file will be exported
#' @param species The taxonomy name/id. Default is "9606" for Homo sapiens.
#' @importFrom httr content
#' @importFrom readr read_tsv
#' @return Either a TSV file stored in the user's file system and its
#' corresponding `data.frame` object in R or and R object containing that
#' information.
#' @examples
#' \donttest{
#' sc <- DISCBIO(valuesG1msReduced)
#' sc <- NoiseFiltering(sc, percentile=0.9, CV=0.2, export=FALSE)
#' sc <- Normalizedata(
#'     sc, mintotal=1000, minexpr=0, minnumber=0, maxexpr=Inf, downsample=FALSE,
#'     dsn=1, rseed=17000
#' )
#' sc <- FinalPreprocessing(sc, export=FALSE, quiet=TRUE)
#' sc <- Clustexp(sc, cln=3, quiet=TRUE) # K-means clustering
#' sc <- comptSNE(sc, max_iter=100, quiet=TRUE)
#' dff <- DEGanalysis2clust(sc, K=3, export=FALSE, quiet=TRUE, plot=FALSE)
#' PPI(dff$FinalDEGsU[, "Gene Name"])
#' }
PPI <- function(data, FileName = NULL, species = "9606") {
    # Save base enpoint as variable
    string_api_url <- "https://string-db.org/api/"
    output_format <- "tsv" #"json", "tsv-no-header", "tsv", "xml"
    method <- "network"
    your_identifiers <- ""
    optional_parameters <- ""
    # Construct API request
    genes <- data
    repos <- GET(
        url = paste0(
            string_api_url,
            output_format,
            '/',
            method,
            '?identifiers=',
            paste(as.character(data), collapse = "%0d"),
            "&",
            "species=",
            species
        )
    )
    cat(
        "Examine response components =",
        status_code(repos),
        "\t",
        "(200 means successful)",
        "\n"
    )
    # Process API request content
    repo_content <- content(repos)
    results  <- read_tsv(repo_content)
    if (!is.null(FileName)) {
        write.csv(results, file = paste0("PPI-", FileName, ".csv"))
    }
    return(results)
}