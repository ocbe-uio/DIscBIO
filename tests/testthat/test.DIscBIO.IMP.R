rm(list = ls())#TEMP
context("Data pre-processing")

# The "valuesG1ms" is the only needed dataset
sc <- DISCBIO(valuesG1ms)
test_that("Data signature changes", {
    expect_equal(class(sc)[1], "DISCBIO") # TODO: change to DIscBIO
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

context("K-means clustering")

sc <- Clustexp(sc, cln=3, quiet=TRUE) # K-means clustering
sc <- comptSNE(sc, rseed=15555, quiet=TRUE)

test_that("tSNE is computed", {
    expect_equal(class(sc@tsne), "data.frame")
    expect_output(str(sc@tsne), "94 obs. of  2 variables")
})

# Jaccard Similarity
test_that("Similarities are as expexted", {
    expect_equal(
        object = Jaccard(sc, Clustering="K-means", K=3, plot = FALSE),
        expected = c(.680, .766, .681)
    )
})

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

test_that("Outliers are the expected", {
    expect_equivalent(Outliers2, c(1, 13, 22))
})

Order <- KmeanOrder(sc, quiet = TRUE, export = FALSE)
# plotOrderKMtsne(sc)

###############################################################################
############################# ORIGINAL CODE BELOW #############################
###############################################################################

# ####Cellular pseudo-time ordering based on k-means clusters 

# ##### Plotting the expression of a particular gene in tSNE map
# g='ENSG00000001460'
# plotExptSNE(sc,g)

# ######## Plotting the K-means clusters in heatmap
# KMclustheatmap(sc,hmethod="single", plot = TRUE) 

# cdiff<-KMClustDiffGenes(sc,K=3,fdr=.01)    ##########3 Binomial differential expression analysis

# cdiff<-DEGanalysis(sc,Clustering="K-means",K=3,fdr=0.1,name="Name",export = TRUE)   ####### differential expression analysis between all clusters
# cdiff<-DEGanalysisM(sc,Clustering="K-means",K=3,fdr=0.1,name="Name",First="CL1",Second="CL2",export = TRUE)     ####### differential expression analysis between particular clusters.

# ############ Plotting the DEGs
# name<-cdiff[[2]][1,6]     # From the table of the differential expression analysis between all pairs of clusters
# U<-read.csv(file=paste0(name),head=TRUE,sep=",")
# Vplot<-VolcanoPlot(U,value=0.05,name=name,adj=FALSE,FS=.4)


# ### Decision tree
# sigDEG<-cdiff[[1]]     
# DATAforDT<-ClassVectoringDT(sc,Clustering="K-means",K=3,First="CL1",Second="CL2",sigDEG)
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


# ##########################################################################################################################
# ##########################################################################################################################
# #################################################        Model-Based Clustering     ##########################################
# ##########################################################################################################################
# ##########################################################################################################################

# sc <- ExprmclustMB(sc,clusternum =3,reduce = T,quiet = FALSE)    ########## #TODO: Maybe here it is better to change the "clusternum" into "K"

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
# cdiff<-DEGanalysisM(sc,Clustering="MB",K=3,fdr=0.1,name="Name",First="CL1",Second="CL2",export = TRUE)     ####### differential expression analysis between particular clusters.

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