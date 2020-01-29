

sc<- PSCANseq(valuesG1ms)    ############### The "valuesG1ms" is the only needed dataset
sc<-NoiseFiltering(sc)       ############### This function will be used only if the dataset has ERCC
sc<-Normalizedata(sc, mintotal=1000, minexpr=0, minnumber=0, maxexpr=Inf, downsample=FALSE, dsn=1, rseed=17000) #### In this case this function is used to normalize the reads.
sc<-Normalizedata(sc, mintotal=?, minexpr=?, minnumber=?, maxexpr=Inf, downsample=FALSE, dsn=1, rseed=17000) #### In this case this function is used to filter out genes and cells, numbers should be used for the ?.
                                                                                                             #### This function can be used for: 1- filtering and normalizing the dataset that has no ERCC. 2- to normalize and filter genes and cells after the noise filtering.

sc<-FinalPreprocessing(sc,GeneFlitering="NoiseF",export = TRUE) ### The GeneFiltering can be either "NoiseF" or"ExpF"




##########################################################################################################################
##########################################################################################################################
#################################################        K-means Clustering     ##########################################
##########################################################################################################################
##########################################################################################################################


sc<- Clustexp(sc,cln=3)    #### K-means clustering
plotGap(sc)        ### Plotting gap statistics

############ Plotting K-means clusters
sc<- comptSNE(sc,rseed=15555)
plottSNE(sc)
plotKmeansLabelstSNE(sc) # To plot the the ID of the cells in eacj cluster
plotSymbolstSNE(sc,types=sub("(\\_\\d+)$","", names(sc@ndata)))    


# Silhouette of k-means clusters
par(mar=c(6,2,4,2))
plotSilhouette(sc,K=3)  

# Jaccard Similarity
Jaccard(sc,Clustering="K-means", K=3, plot = TRUE) 


#################Defining outliers in K-means clustering

Outliers<- FindOutliersKM(sc, K=3, outminc=5,outlg=2,probthr=.5*1e-3,thr=2**-(1:40),outdistquant=.75, plot = TRUE, quiet = FALSE)

########## In case the user would like to adjust the outlg
outlg<-round(length(sc@fdata[,1])/200)     # The cell will be considered as an outlier if it has a minimum of 0.5% of the number of filtered genes as outlier genes. 
Outliers<- FindOutliersKM(sc, K=3, outminc=5,outlg=outlg,probthr=.5*1e-3,thr=2**-(1:40),outdistquant=.75, plot = TRUE, quiet = FALSE)


####Cellular pseudo-time ordering based on k-means clusters 
Order<-KmeanOrder(sc,quiet = FALSE, export = TRUE)
plotOrderKMtsne(sc)

##### Plotting the expression of a particular gene in tSNE map
g='ENSG00000001460'
plotExptSNE(sc,g)



######## Plotting the K-means clusters in heatmap
KMclustheatmap(sc,hmethod="single", plot = TRUE) 

cdiff<-KMClustDiffGenes(sc,K=3,fdr=.01)    ##########3 Binomial differential expression analysis

cdiff<-DEGanalysis(sc,Clustering="K-means",K=3,fdr=0.1,name="Name",export = TRUE)   ####### differential expression analysis between all clusters
cdiff<-DEGanalysisM(sc,Clustering="K-means",K=3,fdr=0.1,name="Name",First="CL1",Second="CL2",export = TRUE)     ####### differential expression analysis between particular clusters.

############ Plotting the DEGs
name<-cdiff[[2]][1,6]     # From the table of the differential expression analysis between all pairs of clusters
U<-read.csv(file=paste0(name),head=TRUE,sep=",")
Vplot<-VolcanoPlot(U,value=0.05,name=name,adj=FALSE,FS=.4)


### Decision tree
sigDEG<-cdiff[[1]]     
DATAforDT<-ClassVectoringDT(sc,Clustering="K-means",K=3,First="CL1",Second="CL2",sigDEG)
j48dt<-J48DT(DATAforDT)
summary(j48dt) 
j48dt<-J48DTeval(DATAforDT,num.folds=10,First="CL1",Second="CL2")
rpartDT<-RpartDT(DATAforDT)
rpartEVAL<-RpartEVAL(DATAforDT,num.folds=10,First="CL1",Second="CL2")

########## Networking analysis
DEGs<-cdiff[[2]][1,6]     # From the table of the differential expression analysis between all pairs of clusters
data<-read.csv(file=paste0(DEGs),head=TRUE,sep=",")
data<-data[,3]

FileName=paste0(DEGs)

ppi<-PPI(data,FileName)
ppi

networking<-NetAnalysis(ppi)
networking                

FileName="Up.DownDEG" 
network<-Networking(data,FileName)







##########################################################################################################################
##########################################################################################################################
#################################################        Model-Based Clustering     ##########################################
##########################################################################################################################
##########################################################################################################################

sc <- ExprmclustMB(sc,clusternum =3,reduce = T,quiet = FALSE)    ########## Maybe here it is better to change the "clusternum" into "K"

########## plotting the clusters in PCA
PlotmclustMB(sc)
PCAplotSymbols(sc)

# Plotting the model-based clusters in tSNE maps
sc<- comptsneMB(sc,rseed=15555)
plottsneMB(sc)
plotMBLabelstSNE(sc)

# Silhouette of Model-based clusters
par(mar=c(6,2,4,2))
plotsilhouetteMB(sc,K=3)


# Jaccard Similarity
Jaccard(sc,Clustering="MB", K=3, plot = TRUE) 


#################Defining outliers in Model-based clustering

Outliers<- FindOutliersMB(sc, K=3, outminc=5,outlg=2,probthr=.5*1e-3,thr=2**-(1:40),outdistquant=.75, plot = TRUE, quiet = FALSE)

########## In case the user would like to adjust the outlg
outlg<-round(length(sc@fdata[,1])/200)     # The cell will be considered as an outlier if it has a minimum of 0.5% of the number of filtered genes as outlier genes. 
Outliers<- FindOutliersMB(sc, K=3, outminc=5,outlg=outlg,probthr=.5*1e-3,thr=2**-(1:40),outdistquant=.75, plot = TRUE, quiet = FALSE)


###############Cellular pseudo-time ordering based on Model-based clusters
sc<-MB_Order(sc,quiet = FALSE, export = TRUE)
PlotMBorderPCA(sc)




##### Plotting the expression of a particular gene in tSNE map
sc<-comptsneMB(sc)
g='ENSG00000001460'
plotexptsneMB(sc,g)



######## Plotting the Model-based clusters in heatmap
MBclustheatmap(sc,hmethod="single", plot = TRUE) 

###################
cdiff<-MBClustDiffGenes(sc,K=3,fdr=.01)    ##########3 Binomial differential expression analysis

cdiff<-DEGanalysis(sc,Clustering="MB",K=3,fdr=0.1,name="Name",export = TRUE)   ####### differential expression analysis between all clusters
cdiff<-DEGanalysisM(sc,Clustering="MB",K=3,fdr=0.1,name="Name",First="CL1",Second="CL2",export = TRUE)     ####### differential expression analysis between particular clusters.

############ Plotting the DEGs
name<-cdiff[[2]][1,6]     # From the table of the differential expression analysis between all pairs of clusters
U<-read.csv(file=paste0(name),head=TRUE,sep=",")
Vplot<-VolcanoPlot(U,value=0.05,name=name,adj=FALSE,FS=.4)


### Decision tree
sigDEG<-cdiff[[1]]     
DATAforDT<-ClassVectoringDT(sc,Clustering="MB",K=3,First="CL1",Second="CL2",sigDEG)
j48dt<-J48DT(DATAforDT)
summary(j48dt) 
j48dt<-J48DTeval(DATAforDT,num.folds=10,First="CL1",Second="CL2")
rpartDT<-RpartDT(DATAforDT)
rpartEVAL<-RpartEVAL(DATAforDT,num.folds=10,First="CL1",Second="CL2")

########## Networking analysis
DEGs<-cdiff[[2]][1,6]     # From the table of the differential expression analysis between all pairs of clusters
data<-read.csv(file=paste0(DEGs),head=TRUE,sep=",")
data<-data[,3]

FileName=paste0(DEGs)

ppi<-PPI(data,FileName)
ppi

networking<-NetAnalysis(ppi)
networking                

FileName="Up.DownDEG" 
network<-Networking(data,FileName)






