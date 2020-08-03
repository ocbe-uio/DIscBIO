# ---------------------------------------------------------------------------- #
#                              Data pre-processing                             #
# ---------------------------------------------------------------------------- #

context("Data loading and pre-processing")

sc <- DISCBIO(valuesG1msRed)  # Reduced dataset used for testing

test_that("Loading datasets generate the expected output", {
    expect_equal(dim(valuesG1msRed), c(1000, 21))
})

test_that("Data signature changes", {
    expect_equal(class(sc)[1], "DISCBIO")
    expect_equal(attr(class(sc), "package"), "DIscBIO")
})

# This function will be used only if the dataset has ERCC
sc <- NoiseFiltering(sc, plot=FALSE, export=FALSE, quiet=TRUE)

test_that("Noise filtering is added", {
    expect_equal(length(sc@noiseF), 272)
})

# In this case this function is used to normalize the reads
sc <- Normalizedata(
    sc, mintotal=1000, minexpr=0, minnumber=0, maxexpr=Inf, downsample=FALSE,
    dsn=1, rseed=17000
)

test_that("Data is normalized", {
    expect_equal(class(sc@fdata), "data.frame")
    expect_output(str(sc@fdata), "908 obs. of  21 variables")
})

# This function can be used for: 1- filtering and normalizing the dataset that has no ERCC. 2- to normalize and filter genes and cells after the noise filtering.
sc <- FinalPreprocessing(sc, GeneFlitering="NoiseF", export=FALSE, quiet=TRUE)

test_that("Data is normalized", {
    expect_equal(dim(sc@fdata), c(272, 21))
})

# ---------------------------------------------------------------------------- #
#                              K-means clustering                              #
# ---------------------------------------------------------------------------- #

context("K-means clustering")

sc <- Clustexp(sc, clustnr=2, cln=2, bootnr=10, quiet=TRUE, rseed=17000)
sc <- comptSNE(sc, rseed=15555, quiet=TRUE, max_iter=1, epoch=10)

test_that("tSNE is computed", {
    expect_equal(class(sc@tsne), "data.frame")
    expect_output(str(sc@tsne), "21 obs. of  2 variables")
})

test_that("Cluster plots output is as expexted", {
    expect_equivalent(
        object = Jaccard(sc, Clustering="K-means", K=2, plot=FALSE, R=5),
        expected = c(.327, .425),
        tolerance = .01
    )
    expect_equal(
        object = clustheatmap(sc, hmethod = "single", plot = FALSE),
        expected = c(1, 2)
    )
})

# --------------------------------- Outliers --------------------------------- #

context("Outliers")

Outliers <- FindOutliers(
    sc, K=2, outminc=5, outlg=2, probthr=.5*1e-3, thr=2 ^ (-1:-40),
    outdistquant=.75, plot = FALSE, quiet = TRUE
)
Order <- KmeanOrder(sc, quiet = TRUE, export = FALSE)

test_that("Outliers are the expected", {
    expect_equivalent(Outliers, c(3, 7, 14))
    expect_equivalent(
        object = Order@kordering,
        expected = c(
            15, 13, 5, 14, 20, 18, 16, 21, 17, 9, 8, 12, 3, 10, 2, 11, 7, 19,
            6, 4, 1
        )
    )
})

# --------------------- Differential Expression Analysis --------------------- #

context("Differential Expression Analysis")

# differential expression analysis between two particular clusters.
cdiff3 <- suppressMessages(
    DEGanalysis2clust(
        sc, Clustering="K-means", K=2, fdr=.2, name="Name", First="CL1",
        Second="CL2", export=FALSE, quiet=TRUE, plot=FALSE, nperms=2, nresamp=2
    )
)

test_that("DEGs are calculated", {
    expect_equivalent(
        object = sapply(cdiff3, function(x) class(x)[1]),
        expected = c("matrix", "data.frame", "data.frame", "data.frame")
    )
})
