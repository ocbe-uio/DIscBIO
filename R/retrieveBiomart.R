#' @title Retrieve data from BioMart
#' @description uses functions from the biomaRt package to retrieve dataframes from the BioMart Database
#' @details Since the BioMart database is not always accessible, this function envelops the requests to the database in a set of tryCatch functions to allow for multiple queries and easier feedback to the end user
#' @param gene_name gene signature
#' @param quiet if `TRUE`, suppresses messages
#' @param max_tries maximum number of times the function will try to reach the database
#' @importFrom biomaRt useDataset useMart useEnsembl getBM
#' @return data.frame resulting from a successful call to getBM.
retrieveBiomart <- function(gene_name, quiet = FALSE, max_tries = 3) {
    # Generates a Mart object
    if (!quiet) message("Retrieving mart object. Please wait.")
    mart <- tryCatch({
            useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
        },
        error = function(err) {
            message(
                "The Mart object could not be retrieved.",
                "Will try again using a different function.",
                "Here is the original error message:\n",
                err
            )
            return(NULL)
        }
    )
    if (is.null(mart)) {
        mart <- tryCatch( {
                useDataset("hsapiens_gene_ensembl", useEnsembl("ensembl"))
            },
            error = function(err) {
                message(err)
                stop(
                    "The Mart object still could not be retrieved. ",
                    "The server may be down. Please try again later."
                )
                return(NULL)
            }
        )
    }

    # Retrieve BioMart dataframe
    if (!quiet) {
            message(
            "Done.\n",
            "Accessing the BioMart database.",
            "This operation may take a few minutes. Please wait."
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
                # TODO: check if this should really be wrapped around quiet. 
                # Warnings and errors should be visible no matter what!?
                if (!quiet) {
                    message(
                        "The BioMart database could not be reached. ",
                        "Retrying (", tries, "/", max_tries, "). ",
                        "Here is the original error message:\n",
                        err
                        )
                }
                tries <- tries + 1
                return(NULL)
            },
            warning = function(warn) {
                if (!quiet) {
                    message(
                        "The BioMart database could not be reached. ",
                        "Retrying (", tries, "/", max_tries, "). ",
                        "Here is the original warning:\n",
                        warn
                    )
                }
                tries <- tries + 1
                return(NULL)
            }
        )
    }

    # Giving up or returning output
    if (is.null(G_list) & tries > 3) {
        stop(
            "The BioMart database could not be reached after",
            max_tries,
            "tries. ",
            "Please try again later."
        )
    } else {
        if (!quiet) message("Done.")
    }

    return(G_list)
}