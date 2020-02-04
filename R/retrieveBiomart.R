#' @title Retrieve data from BioMart
#' @description uses functions from the biomaRt package to retrieve dataframes from the BioMart Database
#' @details Since the BioMart database is not always accessible, this function envelops the requests to the database in a set of tryCatch functions to allow for multiple queries and easier feedback to the end user
#' @param gene_name gene signature
#' @param quiet if `TRUE`, suppresses messages
#' @param max_tries maximum number of times the function will try to reach the database
#' @return data.frame resulting from a successful call to getBM.
retrieveBiomart <- function(gene_name, quiet = FALSE, max_tries = 3) {
    # Generates a Mart object
    if (!quiet) {
        message("Retrieving mart object. Please wait.")
    }
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

    # Retrieve BioMart dataframe
    if (!quiet) {
            message(
            "Done.\n",
            "Accessing BioMart database. This operation may take minutes. ",
            "Please wait."
        )
    }
    G_list <- NULL
    tries <- 1
    while (is.null(G_list) & tries <= max_tries) {
        G_list <- tryCatch({
                if (quiet) {
                    suppressMessages(
                        getBM(
                            filters = "ensembl_gene_id",
                            attributes = c("ensembl_gene_id","hgnc_symbol"),
                            values = gene_name,
                            mart = mart,
                            useCache = FALSE
                        )
                    )
                } else {
                    getBM(
                        filters = "ensembl_gene_id",
                        attributes = c("ensembl_gene_id","hgnc_symbol"),
                        values = gene_name,
                        mart = mart,
                        useCache = FALSE
                    )
                }
            },
            error = function(err) {
                warning("The BioMart database could not be reached.")
                warning("Retrying (", tries, "/", max_tries, ")")
                warning("Here is the original error message:")
                warning(err)
                tries <- tries + 1
                return(NULL)
            }
        )
    }

    # Giving up
    if (is.null(G_list) & tries > 3) {
        warning(
            "The BioMart database could not be reached after",
            max_tries,
            "tries."
        )
        stop("Please try again later.")
    } else {
        if (!quiet) message("Done.")
    }
    return(G_list)
}