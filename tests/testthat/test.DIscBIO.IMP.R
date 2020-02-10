# ---------------------------------------------------------------------------- #
#                              Data pre-processing                             #
# ---------------------------------------------------------------------------- #

context("Data loading and pre-processing")

# The "valuesG1ms" is the only needed dataset
sc <- DISCBIO(valuesG1ms)

test_that("Loading datasets generate the expected output", {
    expect_equal(dim(valuesG1ms), c(59838, 94))
})

test_that("Data signature changes", {
    expect_equal(class(sc)[1], "DISCBIO")
    expect_equal(attr(class(sc), "package"), "DIscBIO")
})

# This function will be used only if the dataset has ERCC
sc <- NoiseFiltering(sc, plot=FALSE, export=FALSE, quiet=TRUE)
test_that("Noise filtering is added", {
    expect_equal(length(sc@noiseF), 5684)
})

# In this case this function is used to normalize the reads
sc <- Normalizedata(
    sc, mintotal=1000, minexpr=0, minnumber=0, maxexpr=Inf, downsample=FALSE, 
    dsn=1, rseed=17000
)

test_that("Data is normalized", {
    expect_equal(class(sc@fdata), "data.frame")
    expect_output(str(sc@fdata), "59746 obs. of  94 variables")
})

# This function can be used for: 1- filtering and normalizing the dataset that has no ERCC. 2- to normalize and filter genes and cells after the noise filtering.
sc <- FinalPreprocessing(sc, GeneFlitering="NoiseF", export=FALSE, quiet=TRUE)

test_that("Data is normalized", {
    expect_equal(dim(sc@fdata), c(5684, 94))
})

# ---------------------------------------------------------------------------- #
#                              K-means clustering                              #
# ---------------------------------------------------------------------------- #

context("K-means clustering")

sc <- Clustexp(sc, cln=3, quiet=TRUE) # K-means clustering
sc <- comptSNE(sc, rseed=15555, quiet=TRUE)

test_that("tSNE is computed", {
    expect_equal(class(sc@tsne), "data.frame")
    expect_output(str(sc@tsne), "94 obs. of  2 variables")
})

test_that("Cluster plots output is as expexted", {
    expect_equal(
        object = Jaccard(sc, Clustering="K-means", K=1, plot = FALSE),
        expected = .68
    )
    expect_equal(
        object = KMclustheatmap(sc, hmethod = "single", plot = FALSE),
        expected = c(3, 1, 2)
    )
})

# --------------------------------- Outliers --------------------------------- #

context("Outliers")

Outliers <- FindOutliersKM(
    sc, K=3, outminc=5, outlg=2, probthr=.5*1e-3, thr=2**-(1:40),
    outdistquant=.75, plot = FALSE, quiet = TRUE
)

# Adjusting outliers
outlg <- round(length(sc@fdata[, 1]) / 200) # The cell will be considered as an outlier if it has a minimum of 0.5% of the number of filtered genes as outlier genes.
Outliers2 <- FindOutliersKM(
    sc, K=3, outminc=5, outlg=outlg, probthr=.5*1e-3, thr=2**-(1:40),
    outdistquant=.75, plot = FALSE, quiet = TRUE
)
Order <- KmeanOrder(sc, quiet = TRUE, export = FALSE)

test_that("Outliers are the expected", {
    expect_equivalent(
        object = Outliers,
        expected = c(
            1, 3:9, 11:15, 17:25, 27, 29, 35, 37, 40, 43, 45, 52, 58, 61:63, 66,
            71:75, 78:79, 82, 85, 92
        )
    )
    expect_equivalent(Outliers2, c(1, 13, 22))
    expect_equivalent(
        object = Order@kordering, 
        expected = c(
            23, 72, 29, 11, 19, 26, 48, 77, 80, 83, 82, 16, 88, 5, 13, 24, 86, 
            42, 34, 45, 87, 15, 65, 27, 10, 58, 89, 35, 71, 62, 32, 44, 38, 56, 
            55, 50, 53, 2, 49, 59, 54, 61, 60, 1, 43, 9, 31, 6, 18, 64, 46, 41, 
            22, 20, 33, 37, 63, 90, 52, 70, 66, 92, 39, 51, 30, 93, 12, 14, 75, 
            78, 47, 28, 91, 68, 21, 73, 94, 25, 17, 8, 74, 57, 69, 4, 84, 7, 40,
            81, 85, 79, 67, 3, 76, 36
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
    sc, Clustering="K-means", K=3, fdr=.15, name="Name", First="CL1",
    Second="CL2", export=FALSE, quiet=TRUE, plot=FALSE
)

test_that("DEGs are calculated", {
    expect_identical(sapply(cdiff1, class), c("matrix", "data.frame"))
    expect_identical(sapply(cdiff2, class), c("matrix", "data.frame"))
    expect_identical(sapply(cdiff3, class), c("matrix", "data.frame"))
})

# Decision tree
sigDEG <- cdiff3[[1]]

DATAforDT <- ClassVectoringDT(
    sc, Clustering="K-means", K=3, First="CL1", Second="CL2", sigDEG,
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
    expect_output(str(DATAforDT), "5 obs. of  79 variables")
    expect_s3_class(j48dt, "J48")
    expect_s3_class(summary(j48dt), "Weka_classifier_evaluation")
    expect_identical(j48dt_eval, c(TP = 18, FN = 11, FP = 6, TN = 44))
    expect_s3_class(rpartDT, "rpart")
    expect_identical(rpartEVAL, c(TP = 19, FN = 10, FP = 7, TN = 43))
})

# ---------------------------------------------------------------------------- #
#                            Model-based clustering                            #
# ---------------------------------------------------------------------------- #

context("Model-based clustering")

# Technically, this should be done before Clustexp, but it's ok in practice to 
# apply it after K-means because it uses different slots.
sc <- Exprmclust(sc, quiet = TRUE)

test_that("Model-based clustering elements are OK", {
    expect_identical(
        object = names(sc@MBclusters),
        expected = c("pcareduceres", "MSTtree", "clusterid", "clucenter")
    )
})

sc <- comptsneMB(sc, rseed=15555, quiet = TRUE)

test_that("tSNE clustering works fine", {
    expect_equal(dim(sc@MBtsne), c(94, 2))
})

# --------------------------------- Outliers --------------------------------- #

context("MB outliers")

sc <- Clustexp(sc, cln=3, quiet=TRUE)

Outliers <- FindOutliersMB(
    sc, K=3, outminc=5, outlg=2, probthr=.5*1e-3, thr=2**-(1:40),
    outdistquant=.75, plot = FALSE, quiet = TRUE
)
outlg <- round(length(sc@fdata[, 1]) / 200) # The cell will be considered as an outlier if it has a minimum of 0.5% of the number of filtered genes as outlier genes. 
Outliers2 <- FindOutliersMB(
    sc, K=3, outminc=5, outlg=outlg, probthr=.5*1e-3, thr=2**-(1:40),
    outdistquant=.75, plot = FALSE, quiet = TRUE
)

test_that("MB clustering and outliers work as expected", {
    expect_equal(
        object = Jaccard(sc, Clustering="MB", K=1, plot = FALSE),
        expected = .66
    )
    expect_equivalent(
        object = Outliers,
        expected = c(
            1, 3, 4, 5, 7, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 22, 23, 
            24, 25, 27, 29, 35, 40, 43, 49, 54, 58, 61, 63, 68, 72, 73, 75, 78, 
            79, 82, 85, 94
        )
    )
    expect_equal(outlg, 28)
    expect_equal(Outliers2, c("G1" = 1, "G1.12" = 13, "G1.21" = 22))
})

sc <- MB_Order(sc, quiet = TRUE, export = FALSE)
mb_heat <- MBclustheatmap(sc, hmethod="single", plot = FALSE, quiet = TRUE)

test_that("More MB things are OK", {
    expect_equal(
        object = sc@MBordering,
        expected = c(
            39, 29, 56, 35, 25, 52, 54, 23, 37, 26, 48, 46, 40, 38, 49, 33, 18, 
            44, 47, 70, 32, 36, 57, 42, 53, 6, 27, 71, 41, 34, 88, 64, 94, 86, 
            76, 72, 51, 59, 68, 60, 84, 78, 89, 79, 61, 63, 81, 58, 92, 80, 93, 
            85, 83, 62, 69, 77, 67, 75, 90, 82, 65, 73, 19, 17, 50, 91, 7, 11, 
            15, 12, 28, 66, 20, 87, 21, 4, 3, 45, 43, 2, 8, 5, 1, 55, 9, 10, 24,
            14, 74, 16, 13, 22, 30, 31
        )
    )
    expect_equal(mb_heat, c(2, 1, 3))
})

# ----------------------------------- DEGs ----------------------------------- #

context("MB DEGs")

# Binomial DE analysis
cdiff1 <- MBClustDiffGenes(sc, K=1, fdr=.2, export=FALSE, quiet=TRUE)

# DE analysis between all clusters
cdiff2 <- DEGanalysis(
    sc, Clustering="MB", K=2, fdr=.2, name="Name", export=FALSE, quiet=TRUE,
    plot = FALSE
)

# differential expression analysis between particular clusters.
cdiff3 <- DEGanalysis2clust(
    sc, Clustering="MB", K=2, fdr=.15, name="Name", First="CL1", Second="CL2",
    export = FALSE, plot = FALSE, quiet = TRUE
)

test_that("DEGs are calculated", {
    expect_identical(sapply(cdiff1, class), c("matrix", "data.frame"))
    expect_identical(sapply(cdiff2, class), c("matrix", "data.frame"))
    expect_identical(sapply(cdiff3, class), c("matrix", "data.frame"))
})

# Decision tree   
sigDEG <- cdiff3[[1]]
DATAforDT <- ClassVectoringDT(
    sc, Clustering="MB", K=3, First="CL1", Second="CL2", sigDEG, quiet = TRUE
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
    expect_output(str(DATAforDT), "340 obs. of  65 variables") # used to be 31
    expect_s3_class(j48dt, "J48")
    expect_s3_class(summary(j48dt), "Weka_classifier_evaluation")
    expect_identical(j48dt_eval, c(TP = 40, FN = 8, FP = 10, TN = 7))
    expect_s3_class(rpartDT, "rpart")
    expect_identical(rpartEVAL, c(TP = 38, FN = 10, FP = 14, TN = 3))
})