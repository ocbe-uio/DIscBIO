context("Stress-testing DEG analyses")
# ==============================================================================
# Determining contants
# ==============================================================================
n_genes <- 1e3
n_clust <- 3
manual_test <- FALSE
# ==============================================================================
# Defining datasets to be used
# ==============================================================================
if (manual_test) {
    data(pan_indrop_matrix_8000cells_18556genes)
    sc <- DISCBIO(pan_indrop_matrix_8000cells_18556genes[, seq_len(n_genes)])
} else {
    sc <- DISCBIO(valuesG1msRed)
    sc <- NoiseFiltering(sc, plot=FALSE, export=FALSE, quiet=TRUE)
}

# ==============================================================================
# Performing operations
# ==============================================================================
sc <- Normalizedata(sc)
sc <- FinalPreprocessing(sc, GeneFlitering="ExpF", export=FALSE, quiet=TRUE)

# ------------------------------------------------------------------------------
# DEG Analysis
# ------------------------------------------------------------------------------
sc <- Clustexp(
    sc, clustnr=n_clust, quiet=TRUE, bootnr = 2, B.gap = 2, rseed=17000
)
cdiff <- DEGanalysis(
    sc, Clustering="K-means", K=n_clust, fdr=0.10, name="all_clusters",
    export=FALSE, quiet=TRUE, plot=FALSE, nresamp=2, nperms=2
)
test_that("DEG analysis doesn't freeze", {
    expect_output(str(cdiff), "List of 2")
})