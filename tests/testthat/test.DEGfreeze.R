library(devtools)#TEMP
context("Stress-testing DEG analyses")
# ==============================================================================
# Loading dataset
# ==============================================================================
install()#TEMP
library(DIscBIO)#TEMP
data(pan_indrop_matrix_8000cells_18556genes)

# ==============================================================================
# Determining contants
# ==============================================================================
n_genes <- 500
n_clust <- 2

# ==============================================================================
# Defining datasets to be used
# ==============================================================================
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
sc <- FinalPreprocessing(sc, GeneFlitering="ExpF", export=FALSE, quiet=TRUE)
sc <- Clustexp(sc, cln=n_clust, quiet=TRUE)

# ------------------------------------------------------------------------------
# DEG Analysis
# ------------------------------------------------------------------------------
cdiff <- DEGanalysis(
    sc, Clustering="K-means", K=n_clust, fdr=0.10, name="all_clusters",
    export=FALSE, quiet=FALSE, plot=FALSE, nresamp=5, nperms=10
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