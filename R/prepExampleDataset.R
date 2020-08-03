#' @title Prepare Example Dataset
#' @description Internal function that prepares a pre-treated dataset for use in
#' several examples
#' @param dataset Dataset used for transformation
#' @param save save results?
#' @details This function serves the purpose of treating datasets such as
#' valuesG1msReduced to reduce examples of other functions by bypassing some
#' analysis steps covered in the vignettes.
#' @return Two rda files, ones for K-means clustering and another for
#' Model-based clustering.
#' @author Waldir Leoncio
prepExampleDataset <- function(dataset, save=TRUE) {
	# ==========================================================================
	# Initial data treatment
	# ==========================================================================
	message("Treating dataset")
	sc <- DISCBIO(dataset)
	sc <- NoiseFiltering(
		sc, percentile=0.9, CV=0.2, export=FALSE, plot=FALSE, quiet=TRUE
	)
	sc <- Normalizedata(sc)
	sc <- FinalPreprocessing(sc, export=FALSE, quiet=TRUE)
	# ==========================================================================
	# Clustering
	# ==========================================================================
	message("K-means clustering")
	sc_k <- Clustexp(sc, cln=3, quiet=TRUE)
	sc_k <- comptSNE(sc_k, quiet=TRUE)
	valuesG1msReduced_treated_K <- sc_k
	message("Model-based clustering")
	sc_mb <- Exprmclust(sc, quiet=TRUE)
	sc_mb <- comptSNE(sc_mb, rseed=15555, quiet=TRUE)
	valuesG1msReduced_treated_MB <- sc_mb
	# ==========================================================================
	# Output
	# ==========================================================================
	message("Saving datasets")
	if (save) {
		save(
			valuesG1msReduced_treated_K,
			file = "data/valuesG1msReduced_treated_K.rda"
		)
		save(
			valuesG1msReduced_treated_MB,
			file = "data/valuesG1msReduced_treated_MB.rda"
		)
	} else {
		message("Not saving dataset because (save == FALSE)")
	}
}