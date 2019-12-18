# This test makes sure the package works with respect to the interactive notebook

context("Reproducing Jupyter notebook: loading and pre-processing")

# Loading datasets =============================================================
test_that("Loading datasets generate the expected output", {
    expect_equal(dim(valuesG1ms), c(59838, 94))
    expect_equal(dim(MLSrawWithoutERCC), c(59746, 94))
})

# 1. Data Pre-processing =======================================================

percentile <- 0.8
CV <- 0.3 
Object <- valuesG1ms
gene_names <- GeneNames

# Remove .x from gene ID's, since they can't be handled later on
gene_names <- as.list(sub("*\\..*", "", unlist(gene_names)))

gene_names2 <- head(gene_names, -92)

geneCol  <- "yellow"
FgeneCol <- "black"
erccCol  <- "blue"
noiseF   <- NoiseFiltering(
    Object, percentile, CV, gene_names2, geneCol, FgeneCol, erccCol, Val=T,
    plot = FALSE, export = FALSE, quiet = TRUE
)       # Val=F  will plot all the ERCC spike-ins

sc <- PSCANseq(MLSrawWithoutERCC) # TODO: understand how PSCANseq can be called as a function (this could cause problems in a clean installation)
sc <- Normalizedata(sc, mintotal = 1000)
gene_list <- noiseF
gene_names <- rownames(sc@ndata)
idx_genes <- is.element(gene_names, gene_list)
gene_names2 <- gene_names[idx_genes]
filteredDataset <- sc@ndata[gene_names2, ]

# ASK: should these cats be part of some function output?
# cat("The gene filtering method= Noise filtering","\n","\n") 
# cat("The Filtered Normalized dataset is called: filteredDataset","\n","\n") 
# cat("The filtered Normalized dataset contains:","\n","Genes:",length(filteredDataset[,1]),"\n","cells:",length(filteredDataset[1,]),"\n","\n")
# save(filteredDataset,file="filteredDataset.Rdata")
sc@fdata <- filteredDataset

# ############# Generating a filtered dataset with raw data for DEG analysis
gene_list <- noiseF
gene_names <- rownames(MLSrawWithoutERCC)
idx_genes <- is.element(gene_names, gene_list)
gene_names2 <- gene_names[idx_genes]
LipoNoisFilteredRawDataset <- MLSrawWithoutERCC[gene_names2, ]

# ASK: should these cats be part of some function output?  
# cat("The Noise Filtered Raw dataset contains:","\n","Genes:",length(LipoNoisFilteredRawDataset[,1]),"\n","cells:",length(LipoNoisFilteredRawDataset[1,]))
# # save(LipoNoisFilteredRawDataset,file="LipoNoisFilteredRawDataset.Rdata")

test_that("Data pre-processing results are reproduced", {
    expect_equal(class(noiseF), "character")
    expect_equal(length(noiseF), 5684)
    expect_equal(noiseF[1234], "ENSG00000105852")
    expect_equal(dim(filteredDataset), c(5684, 94))
    expect_equal(dim(LipoNoisFilteredRawDataset), c(5684, 94))
})

# 2. Cellular Clustering and Pseudo Time ordering ==============================

context("Reproducing Jupyter notebook: cell clusters and pseudo-time ordering")

K <- 3  # Number of Clusters
sc_temp <- Clustexp(
    object = sc, clustnr = 20, bootnr = 50, metric = "pearson", do.gap = T,
    SE.method = "Tibs2001SEmax", SE.factor = .25, B.gap = 50, cln = K,
    rseed = 17000, quiet = TRUE
)
sc <- comptSNE(sc_temp, rseed = 15555, quiet = TRUE)
Clusters <- sc@kmeans$kpart
KmeansClusters <- Clusters # To be used for defining DEGs
AllClusters <- sc@cpart
sc@cpart <- Clusters
Clustering <- "K-means" # Jaccard of k-means clusters

# The cell will be considered as an outlier if it has a minimum of 0.5% of the number of filtered genes as outlier genes.
outlg <- round(length(sc@fdata[, 1]) / 200)
Outliers <- FindOutliersKM(
    object = sc, K = K, outminc = 5, outlg = outlg, probthr = .5 * 1e-3,
    thr = 2**-(1:40), outdistquant = 1, plot = FALSE, quiet = TRUE
)

# RemovingOutliers=FALSE     
# RemovingOutliers=TRUE                    # Removing the defined outlier cells based on K-means Clustering

# ASK: should this be part of some function?
# if(RemovingOutliers==TRUE){
#     names(Outliers)=NULL
#     Outliers
#     valuesG1ms=valuesG1ms[-Outliers]
#     MLSrawWithoutERCC=MLSrawWithoutERCC[-Outliers]

#     dim(valuesG1ms)
#     dim(MLSrawWithoutERCC)

#     colnames(valuesG1ms)
#     colnames(MLSrawWithoutERCC)
#     cat("Outlier cells were removed, now you need to start from the beginning")
# }

sampleNames <- colnames(sc@fdata)
Order <- KmeanOrder(
    sc@fdata, Clusters, sampleNames, quiet = TRUE, export = FALSE
)
sc@kordering <- Order
sc@ndata <- rbind(sc@ndata, sc@kordering)
rownames(sc@ndata)[nrow(sc@ndata)] <- "Pseudo-time ordering of Kmeans clustering"  # in case the user wanted to repeat this step the name: "Pseudo-time ordering of Kmeans clustering" should be change each time.
g <- rownames(sc@ndata)[nrow(sc@ndata)]
# plotExptSNE(sc, g)

test_that("Reproduced cellular clustering and pseudo time ordering", {
    expect_output(str(sc_temp@expdata), "59746 obs. of  94 variables")
    expect_equal(as.numeric(table(Clusters)), c(29, 50, 15))
    # ASK: are these differences justifiable?
    expect_equivalent(
        object = Jaccard(sc, Clustering, K, plot = FALSE),
        expected = c(.671, .768, .681),
        tolerance = 0.01)
    expect_equal(as.numeric(Outliers), c(1, 13, 22))
    expect_equal(KMclustheatmap(sc, plot = FALSE), c(3, 1, 2))
})