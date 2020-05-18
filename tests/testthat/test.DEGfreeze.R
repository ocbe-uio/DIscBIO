# library(devtools)#TEMP
# library(testthat)#TEMP
# source("Aux/samr.R")#TEMP
# source("Aux/resample.R")#TEMP
# source("Aux/rankcol.R")#TEMP
# source("Aux/samr/R/samr.morefuns.R")#TEMP
# source("Aux/samr/R/samr.R")#TEMP
# library(DIscBIO)#TEMP
context("Stress-testing DEG analyses")
# ==============================================================================
# Loading dataset
# ==============================================================================
data(pan_indrop_matrix_8000cells_18556genes)

# ==============================================================================
# Determining contants
# ==============================================================================
n_genes <- 1e3
n_clust <- 3

# ==============================================================================
# Defining datasets to be used
# ==============================================================================
message(Sys.time(), " - Subsetting")
sc <- DISCBIO(pan_indrop_matrix_8000cells_18556genes[, seq_len(n_genes)])

# ==============================================================================
# Performing operations
# ==============================================================================
# MIinExp <- mean(rowMeans(dataset, na.rm=TRUE))
# MinNumber <- round(length(dataset[1, ]) / 10) # 2b xpressed in > 10% cells.
# sc <- Normalizedata(
#     sc, mintotal=1000, minexpr=MIinExp, minnumber=MinNumber, maxexpr=Inf,
#     downsample=FALSE, dsn=1, rseed=17000
# )
message("Filtering")
sc <- Normalizedata(sc, mintotal=1000, maxexpr=Inf, downsample=FALSE, dsn=1)
sc <- FinalPreprocessing(sc, GeneFlitering="ExpF", export=FALSE, quiet=TRUE)
sc <- Clustexp(sc, clustnr=n_clust, quiet=TRUE, bootnr = 2, B.gap = 2)

# ------------------------------------------------------------------------------
# DEG Analysis
# ------------------------------------------------------------------------------
message("Analysing DEG")
cdiff <- DEGanalysis(
    sc, Clustering="K-means", K=n_clust, fdr=0.10, name="all_clusters",
    export=FALSE, quiet=TRUE, plot=FALSE, nresamp=2, nperms=2
)
# sc <- Exprmclust(sc, K=K, reduce=TRUE, quiet=TRUE)
# sc <- MB_Order(sc, quiet=TRUE, export=FALSE)
# MBcdiffBinomial <- MBClustDiffGenes(
# 	sc, K=K, pValue=0.05, fdr=0.15, export=FALSE, quiet=FALSE
# )
# sc <- Clustexp(sc, cln=K, quiet=TRUE)
# MBcdiff <- DEGanalysis(
# 	sc, Clustering="MB", K=K, fdr=0.05, name="all_clusters", export = FALSE,
# 	quiet=TRUE
# )
test_that("DEG analysis doesn't freeze", {
    expect_output(str(cdiff), "List of 2")
})