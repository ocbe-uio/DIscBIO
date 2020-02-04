rm(list = ls())#TEMP

# ---------------------------------------------------------------------------- #
#                              Data pre-processing                             #
# ---------------------------------------------------------------------------- #

context("Data pre-processing")

# The "valuesG1ms" is the only needed dataset
sc <- DISCBIO(valuesG1ms)
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
        object = Jaccard(sc, Clustering="K-means", K=3, plot = FALSE),
        expected = c(.680, .766, .681)
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
cdiff1 <- KMClustDiffGenes(sc, K=3, fdr=.2, export=FALSE, quiet=TRUE)

# differential expression analysis between all clusters
cdiff2 <- DEGanalysis(
    sc, Clustering="K-means", K=3, fdr=.2, name="Name", export=FALSE,
    quiet=TRUE, plot=FALSE
)

# differential expression analysis between two particular clusters.
cdiff3 <- DEGanalysis2clust(
    sc, Clustering="K-means", K=3, fdr=.2, name="Name", First="CL1",
    Second="CL2", export=FALSE, quiet=TRUE, plot=FALSE
)

test_that("DEGs are calculated", {
    expect_identical(sapply(cdiff1, class), c("matrix", "data.frame"))
    expect_identical(sapply(cdiff2, class), c("matrix", "data.frame"))
    expect_identical(sapply(cdiff3, class), c("matrix", "data.frame"))
})

# Plotting the DEGs
cdiff <- cdiff3
name <- cdiff[[2]][1, 6] # From the DE analysis table between all cluster pairs
# U <- read.csv(file = paste0(name), head=TRUE, sep=",")
# Vplot <- VolcanoPlot(U, value=0.05, name=name, adj=FALSE, FS=.4)
# TODO: get an example for VolcanoPlot which doesn't involve importing files

# Decision tree
sigDEG <- cdiff[[1]]

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
    expect_output(str(DATAforDT), "31 obs. of  79 variables")
    expect_s3_class(j48dt, "J48")
    expect_s3_class(summary(j48dt), "Weka_classifier_evaluation")
    expect_identical(j48dt_eval, c(TP = 20, FN = 9, FP = 7, TN = 43))
    expect_s3_class(rpartDT, "rpart")
    expect_identical(rpartEVAL, c(TP = 15, FN = 14, FP = 6, TN = 44))
})

# ----------------------------- Network analysis ----------------------------- #

context("Network analysis")
DEGs <- cdiff[[2]][1, 6] # From the DE analysis table between all cluster pairs
# data<-read.csv(file=paste0(DEGs),head=TRUE,sep=",")
# data<-data[,3]
# FileName <- paste0(DEGs)
# ppi <- PPI(data, FileName)
# ppi
# networking<-NetAnalysis(ppi)
# FileName="Up.DownDEG" 
# network<-Networking(data,FileName)
# TODO: find reproducible example that doesn't depend on importing data

# ---------------------------------------------------------------------------- #
#                            Model-based clustering                            #
# ---------------------------------------------------------------------------- #

context("Model-based clustering")
sc <- Exprmclust(sc, K = 3,reduce = T, quiet = TRUE) # ASK: shouldn't this be done before Clustexp?

test_that("Model-based clustering elements are OK", {
    # TODO: add test for Exprmclust
})
###############################################################################
############################# ORIGINAL CODE BELOW #############################
###############################################################################


# ########## plotting the clusters in PCA
# PlotmclustMB(sc)  # TODO: check if this works after scope change
# PCAplotSymbols(sc)

# # Plotting the model-based clusters in tSNE maps
# sc<- comptsneMB(sc,rseed=15555)
# plottsneMB(sc)
# plotMBLabelstSNE(sc)

# # Silhouette of Model-based clusters
# par(mar=c(6,2,4,2))
# plotsilhouetteMB(sc,K=3)


# # Jaccard Similarity
# Jaccard(sc,Clustering="MB", K=3, plot = TRUE) 


# #################Defining outliers in Model-based clustering

# Outliers<- FindOutliersMB(sc, K=3, outminc=5,outlg=2,probthr=.5*1e-3,thr=2**-(1:40),outdistquant=.75, plot = TRUE, quiet = FALSE)

# ########## In case the user would like to adjust the outlg
# outlg<-round(length(sc@fdata[,1])/200)     # The cell will be considered as an outlier if it has a minimum of 0.5% of the number of filtered genes as outlier genes. 
# Outliers<- FindOutliersMB(sc, K=3, outminc=5,outlg=outlg,probthr=.5*1e-3,thr=2**-(1:40),outdistquant=.75, plot = TRUE, quiet = FALSE)


# ###############Cellular pseudo-time ordering based on Model-based clusters
# sc<-MB_Order(sc,quiet = FALSE, export = TRUE)
# PlotMBorderPCA(sc)

# ##### Plotting the expression of a particular gene in tSNE map
# sc<-comptsneMB(sc)
# g='ENSG00000001460'
# plotexptsneMB(sc,g)

# ######## Plotting the Model-based clusters in heatmap
# MBclustheatmap(sc,hmethod="single", plot = TRUE) # TODO: check if scope changes didn't break this 

# ###################
# cdiff<-MBClustDiffGenes(sc,K=3,fdr=.01)    ##########3 Binomial differential expression analysis

# cdiff<-DEGanalysis(sc,Clustering="MB",K=3,fdr=0.1,name="Name",export = TRUE)   ####### differential expression analysis between all clusters
# cdiff<-DEGanalysis2clust(sc,Clustering="MB",K=3,fdr=0.1,name="Name",First="CL1",Second="CL2",export = TRUE)     ####### differential expression analysis between particular clusters.

# ############ Plotting the DEGs
# name<-cdiff[[2]][1,6]     # From the table of the differential expression analysis between all pairs of clusters
# U<-read.csv(file=paste0(name),head=TRUE,sep=",")
# Vplot<-VolcanoPlot(U,value=0.05,name=name,adj=FALSE,FS=.4)


# ### Decision tree
# sigDEG<-cdiff[[1]]     
# DATAforDT<-ClassVectoringDT(sc,Clustering="MB",K=3,First="CL1",Second="CL2",sigDEG)
# j48dt<-J48DT(DATAforDT)
# summary(j48dt) 
# j48dt<-J48DTeval(DATAforDT,num.folds=10,First="CL1",Second="CL2")
# rpartDT<-RpartDT(DATAforDT)
# rpartEVAL<-RpartEVAL(DATAforDT,num.folds=10,First="CL1",Second="CL2")

# ########## Networking analysis
# DEGs<-cdiff[[2]][1,6]     # From the table of the differential expression analysis between all pairs of clusters
# data<-read.csv(file=paste0(DEGs),head=TRUE,sep=",")
# data<-data[,3]

# FileName=paste0(DEGs)

# ppi<-PPI(data,FileName)
# ppi

# networking<-NetAnalysis(ppi)
# networking                

# FileName="Up.DownDEG" 
# network<-Networking(data,FileName)
