#' @title Prepare Example Dataset
#' @description Internal function that prepares a pre-treated dataset for use in
#' several examples
#' @param clustering determines the clustering method (either "K" or "MB").
#' @note This is an internal function, reachable through the terminal by using `DIscBIO:::prepExampleDataset`
#' @return rda file
prepExampleDataset <- function(clustering) {
	# ==========================================================================
	# Initial data treatment
	# ==========================================================================
	message("Treating dataset")
	sc <- DISCBIO(valuesG1msReduced)
	sc <- NoiseFiltering(
		sc, percentile=0.9, CV=0.2, export=FALSE, plot=FALSE, quiet=TRUE
	)
	sc <- Normalizedata(sc)
	sc <- FinalPreprocessing(sc, export=FALSE, quiet=TRUE)
	# ==========================================================================
	# Clustering
	# ==========================================================================
	message("Clustering")
	if (clustering == "K") {
		sc <- Clustexp(sc, cln=3, quiet=TRUE)
		sc <- comptSNE(sc, quiet=TRUE)
		output_name <- "valuesG1msReduced_treated_K"
	} else if (clustering == "MB") {
		sc <- Exprmclust(sc, quiet=TRUE)
		sc <- comptsneMB(sc, rseed=15555, quiet=TRUE)
		output_name <- "valuesG1msReduced_treated_MB"
	}
	# ==========================================================================
	# Output
	# ==========================================================================
	message("Saving")
	save(sc, file = paste0("data/", output_name, ".rda"))
}