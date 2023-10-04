#' @title Defining protein-protein interactions (PPI) over a list of genes,
#' @description This function uses STRING-api. The outcome of STRING analysis
#'   will be stored in comma-separated values files.
#' @export
#' @param data A gene list.
#' @param FileName A string vector showing the name to be used to save the
#'   resulted table. If null, no file will be exported
#' @param species The taxonomy name/id. Default is "9606" for Homo sapiens.
#' @importFrom httr content
#' @importFrom utils read.table write.table
#' @return Either CSV files stored in the user's file system and its
#' corresponding `data.frame` object in R or and R object containing that
#' information.
PPI <- function(data, FileName = NULL, species = "9606") {
  repos <- retrieveURL(data, species, "tsv")
  # Process API request content
  repo_content <- content(repos)
  results <- data.frame(repo_content[, , 1])
  if (!is.null(FileName)) {
    write.csv(results, file = paste0("PPI-", FileName, ".csv"))
  }
  return(results)
}
