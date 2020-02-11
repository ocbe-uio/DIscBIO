#' @title Defining protein-protein interactions (PPI) over a list of genes,
#' @description This function uses STRING-api. The outcome of STRING analysis will be stored in tab separated values (TSV) files.  
#' @export
#' @param data A gene list.
#' @param FileName A string vector showing the name to be used to save the resulted table.
#' @param species The taxonomy name/id. Default is "9606" for Homo sapiens.
#' @importFrom httr content
#' @importFrom readr read_tsv
#' @return A TSV file stored in the user's file system and its corresponding `data.frame` object in R.
#' @examples
#' \dontrun{
#' sc <- DISCBIO(valuesG1msReduced)
#' sc <- NoiseFiltering(sc, percentile=0.9, CV=0.2)
#' sc <- Normalizedata(
#'     sc, mintotal=1000, minexpr=0, minnumber=0, maxexpr=Inf, downsample=FALSE,
#'     dsn=1, rseed=17000
#' )
#' sc <- FinalPreprocessing(sc, GeneFlitering="NoiseF")
#' sc <- Clustexp(sc, cln=3) # K-means clustering
#' sc <- comptSNE(sc, rseed=15555)
#' dff <- DEGanalysis2clust(sc, Clustering="K-means", K=3, fdr=0.1, name="Name")
#' DEGs <- dff[[2]][1, 6]
#' data <- read.csv(file=paste0(DEGs),head=TRUE,sep=",")
#' data <- data[,3]
#' FileName <- paste0(DEGs)
#' PPI(data, FileName)
#' }
PPI<-function(data,FileName,species="9606"){
    # Save base enpoint as variable
    string_api_url <- "https://string-db.org/api/"
    output_format <- "tsv" #"json", "tsv-no-header", "tsv", "xml"
    method <- "network"
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