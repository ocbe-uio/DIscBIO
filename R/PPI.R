#' @title Defining protein-protein interactions (PPI) over a list of genes,
#' @description This function uses STRING-api. The outcome of STRING analysis will be stored in tab separated values (TSV) files.  
#' @export
#' @param data A gene list.
#' @param FileName A string vector showing the name to be used to save the resulted table.
#' @param species The taxonomy name/id. Default is "9606" for Homo sapiens.
#' @importFrom httr content
#' @importFrom readr read_tsv
PPI<-function(data,FileName){
	# Save base enpoint as variable
	string_api_url <- "https://string-db.org/api/"
	output_format <- "tsv" #"json", "tsv-no-header", "tsv", "xml"
	method <- "network"
	species <- "9606"
	your_identifiers <- ""
	optional_parameters <- ""
	# Construct API request
	genes <- data
	repos <- GET(url = paste0(string_api_url,output_format,'/',method,'?identifiers=',
                          paste(as.character(data), collapse = "%0d"),"&", "species=",species))
    cat("Examine response components =",status_code(repos),"\t","200 means successful","\n")
	# Process API request content 
	repo_content <- content(repos)
	results  <- read_tsv(repo_content)
	write.csv(results, file = paste0("PPI-",FileName,".csv"))
	return(results)
}