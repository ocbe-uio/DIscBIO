# ---------------------------------------------------------------------------- #
#                              Data pre-processing                             #
# ---------------------------------------------------------------------------- #

context("Data loading and pre-processing")

sc <- DISCBIO(valuesG1msReduced)  # Reduced dataset used for testing

test_that("Loading datasets generate the expected output", {
    expect_equal(dim(valuesG1msReduced), c(1092, 30))
})

test_that("Data signature changes", {
    expect_equal(class(sc)[1], "DISCBIO")
    expect_equal(attr(class(sc), "package"), "DIscBIO")
})

# This function will be used only if the dataset has ERCC
sc <- NoiseFiltering(sc, plot=FALSE, export=FALSE, quiet=TRUE)

test_that("Noise filtering is added", {
    expect_equal(length(sc@noiseF), 341)
})

# In this case this function is used to normalize the reads
sc <- Normalizedata(
    sc, mintotal=1000, minexpr=0, minnumber=0, maxexpr=Inf, downsample=FALSE,
    dsn=1, rseed=17000
)

test_that("Data is normalized", {
    expect_equal(class(sc@fdata), "data.frame")
    expect_output(str(sc@fdata), "1000 obs. of  30 variables")
})

# This function can be used for: 1- filtering and normalizing the dataset that has no ERCC. 2- to normalize and filter genes and cells after the noise filtering.
sc <- FinalPreprocessing(sc, GeneFlitering="NoiseF", export=FALSE, quiet=TRUE)

test_that("Data is normalized", {
    expect_equal(dim(sc@fdata), c(341, 30))
})

# ---------------------------------------------------------------------------- #
#                              K-means clustering                              #
# ---------------------------------------------------------------------------- #

context("K-means clustering")

sc <- Clustexp(sc, cln=2, quiet=TRUE, rseed=17000) # K-means clustering
sc <- comptSNE(sc, rseed=15555, quiet=TRUE, max_iter=10)

test_that("tSNE is computed", {
    expect_equal(class(sc@tsne), "data.frame")
    expect_output(str(sc@tsne), "30 obs. of  2 variables")
})

test_that("Cluster plots output is as expexted", {
    expect_equivalent(
        object = Jaccard(sc, Clustering="K-means", K=2, plot=FALSE, R=5),
        expected = c(.790, .653),
        tolerance = .01
    )
    expect_equal(
        object = KMclustheatmap(sc, hmethod = "single", plot = FALSE),
        expected = c(1, 2)
    )
})

# --------------------------------- Outliers --------------------------------- #

context("Outliers")

Outliers <- FindOutliersKM(
    sc, K=2, outminc=5, outlg=2, probthr=.5*1e-3, thr=2 ^ (-1:-40),
    outdistquant=.75, plot = FALSE, quiet = TRUE
)

# Adjusting outliers
outlg <- round(length(sc@fdata[, 1]) / 200) # The cell will be considered as an outlier if it has a minimum of 0.5% of the number of filtered genes as outlier genes.
Outliers2 <- FindOutliersKM(
    sc, K=2, outminc=5, outlg=outlg, probthr=.5*1e-3, thr=2 ^ (-1:-40),
    outdistquant=.75, plot = FALSE, quiet = TRUE
)
Order <- KmeanOrder(sc, quiet = TRUE, export = FALSE)

test_that("Outliers are the expected", {
    expect_equivalent(Outliers, c(3, 7, 19))
    expect_equivalent(Outliers2, c(3, 7, 19))
    expect_equivalent(
        object = Order@kordering,
        expected = c(
            23, 20, 6, 21, 27, 26, 24, 28, 10, 19, 15, 25, 16, 8, 14, 13, 22, 4,
            17, 2, 3, 18, 11, 29, 9, 5, 12, 1, 30, 7
        )
    )
})

# --------------------- Differential Expression Analysis --------------------- #

context("Differential Expression Analysis")

# Binomial differential expression analysis
cdiff1 <- KMClustDiffGenes(sc, K=1, fdr=.2, export=FALSE, quiet=TRUE)

# differential expression analysis between all clusters
cdiff2 <- DEGanalysis(
    sc, Clustering="K-means", K=2, fdr=.2, name="Name", export=FALSE,
    quiet=TRUE, plot=FALSE, nperms=5, nresamp=2
)

# differential expression analysis between two particular clusters.
cdiff3 <- DEGanalysis2clust(
    sc, Clustering="K-means", K=2, fdr=.15, name="Name", First="CL1",
    Second="CL2", export=FALSE, quiet=TRUE, plot=FALSE, nperms=5, nresamp=2
)

test_that("DEGs are calculated", {
    expect_identical(
        object = sapply(cdiff1, function(x) class(x)[1]),
        expected = c("matrix", "data.frame")
    )
    expect_identical(
        object = sapply(cdiff2, function(x) class(x)[1]),
        expected = c("matrix", "data.frame")
    )
    expect_equivalent(
        object = sapply(cdiff3, function(x) class(x)[1]),
        expected = c("matrix", "data.frame", "data.frame", "data.frame")
    )
})

# Decision tree
sigDEG <- cdiff3[[1]]

DATAforDT <- ClassVectoringDT(
    sc, Clustering="K-means", K=2, First="CL1", Second="CL2", sigDEG,
    quiet = TRUE
)

j48dt <- J48DT(DATAforDT, quiet = TRUE, plot = FALSE)
j48dt_eval <- J48DTeval(
    DATAforDT, num.folds=10, First="CL1", Second="CL2", quiet=TRUE
)
rpartDT <- RpartDT(DATAforDT, quiet = TRUE, plot = FALSE)
rpartEVAL <- RpartEVAL(
    DATAforDT, num.folds=10, First="CL1", Second="CL2", quiet = TRUE
)

test_that("Decision tree elements are defined", {
    expect_output(str(DATAforDT), "3 obs. of  30 variables")
    expect_s3_class(j48dt, "J48")
    expect_identical(
        attributes(j48dt)$class, c("J48", "Weka_tree", "Weka_classifier")
    )
    expect_identical(j48dt_eval, c(TP = 14, FN = 4, FP = 3, TN = 9))
    expect_s3_class(rpartDT, "rpart")
    expect_identical(rpartEVAL, c(TP = 16, FN = 2, FP = 4, TN = 8))
})

# ---------------------------------------------------------------------------- #
#                            Model-based clustering                            #
# ---------------------------------------------------------------------------- #

context("Model-based clustering")

# Technically, this should be done before Clustexp, but it's ok in practice to
# apply it after K-means because it uses different slots.
sc <- Exprmclust(sc, K=2, quiet=TRUE)

test_that("Model-based clustering elements are OK", {
    expect_identical(
        object = names(sc@MBclusters),
        expected = c("pcareduceres", "MSTtree", "clusterid", "clucenter")
    )
})

sc <- comptsneMB(sc, rseed=15555, quiet=TRUE, max_iter=100)

test_that("tSNE clustering works fine", {
    expect_equal(dim(sc@MBtsne), c(30, 2))
})

# --------------------------------- Outliers --------------------------------- #

context("MB outliers")

sc <- Clustexp(sc, cln=2, quiet=TRUE, rseed=17000)

Outliers <- FindOutliersMB(
    sc, K=2, outminc=5, outlg=2, probthr=.5*1e-3, thr=2 ^ (-1:-40),
    outdistquant=.75, plot = FALSE, quiet = TRUE
)
outlg <- round(length(sc@fdata[, 1]) / 200) # The cell will be considered as an outlier if it has a minimum of 0.5% of the number of filtered genes as outlier genes.
Outliers2 <- FindOutliersMB(
    sc, K=2, outminc=5, outlg=outlg, probthr=.5*1e-3, thr=2 ^ (-1:-40),
    outdistquant=.75, plot = FALSE, quiet = TRUE
)

test_that("MB clustering and outliers work as expected", {
    expect_equivalent(
        object = Jaccard(sc, Clustering="MB", K=2, plot = FALSE, R=5),
        expected = c(.819, .499),
        tolerance = 0.01
    )
    expect_equivalent(Outliers, c(3, 4, 7))
    expect_equal(outlg, 2)
    expect_equal(Outliers2, c("G1_12" = 3, "G1_18" = 4, "G1_21" = 7))
})

sc <- MB_Order(sc, quiet = TRUE, export = FALSE)
mb_heat <- MBclustheatmap(sc, hmethod="single", plot = FALSE, quiet = TRUE)

test_that("More MB things are OK", {
    expect_equal(
        object = sc@MBordering,
        expected = c(
            8, 28, 27, 18, 21, 10, 29, 26, 17, 25, 5, 13, 12, 19, 16, 15, 23,
            20, 30, 14, 7, 9, 24, 22, 3, 2, 4, 6, 1, 11
        )
    )
    expect_equal(mb_heat, c(1, 2))
})

# ----------------------------------- DEGs ----------------------------------- #

context("MB DEGs")

# Binomial DE analysis
cdiff1 <- MBClustDiffGenes(sc, K=2, fdr=.2, export=FALSE, quiet=TRUE)

# DE analysis between all clusters
cdiff2 <- DEGanalysis(
    sc, Clustering="MB", K=2, fdr=.2, name="Name", export=FALSE, quiet=TRUE,
    plot = FALSE, nperms=5, nresamp=2
)

# differential expression analysis between particular clusters.
cdiff3 <- DEGanalysis2clust(
    sc, Clustering="MB", K=2, fdr=.15, name="Name", First="CL1", Second="CL2",
    export = FALSE, plot = FALSE, quiet = TRUE, nperms=5, nresamp=2
)

test_that("DEGs are calculated", {
    expect_identical(
        object = sapply(cdiff1, function(x) class(x)[1]),
        expected = c("matrix", "data.frame")
    )
    expect_identical(
        object = sapply(cdiff2, function(x) class(x)[1]),
        expected = c("matrix", "data.frame")
    )
    expect_equivalent(
        object = sapply(cdiff3, function(x) class(x)[1]),
        expected = c("matrix", "data.frame", "data.frame", "data.frame")
    )
})

# Decision tree
sigDEG <- cdiff3[[1]]
DATAforDT <- ClassVectoringDT(
    sc, Clustering="MB", K=2, First="CL1", Second="CL2", sigDEG, quiet = TRUE
)
j48dt <- J48DT(DATAforDT, quiet = TRUE, plot = FALSE)
j48dt_eval <- J48DTeval(
    DATAforDT, num.folds=10, First="CL1", Second="CL2", quiet = TRUE
)
rpartDT <- RpartDT(DATAforDT, quiet = TRUE, plot = FALSE)
rpartEVAL <- RpartEVAL(
    DATAforDT, num.folds=10, First="CL1", Second="CL2", quiet = TRUE
)

test_that("Decision tree elements are defined", {
    expect_output(str(DATAforDT), "29 obs. of  30 variables") # used to be 31
    expect_s3_class(j48dt, "J48")
    expect_identical(
        attributes(j48dt)$class, c("J48", "Weka_tree", "Weka_classifier")
    )
    expect_identical(j48dt_eval, c(TP = 21, FN = 2, FP = 4, TN = 3))
    expect_s3_class(rpartDT, "rpart")
    expect_identical(rpartEVAL, c(TP = 20, FN = 3, FP = 4, TN = 3))
})
