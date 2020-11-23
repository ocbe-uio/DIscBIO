library(DIscBIO)

FileName<-"CTCdataset"        # Name of the dataset
#CSV=TRUE                     # If the dataset has ".csv", the user shoud set CSV to TRUE
CSV=FALSE                     # If the dataset has ".rda", the user shoud set CSV to FALSE

if (CSV==TRUE){
    DataSet <- read.csv(file = paste0(FileName,".csv"), sep = ",",header=T)
    rownames(DataSet)<-DataSet[,1]
    DataSet<-DataSet[,-1]
} else{
    load(paste0(FileName,".rda"))
    DataSet<-get(FileName)
}
cat(paste0("The ", FileName," contains:","\n","Genes: ",length(DataSet[,1]),"\n","cells: ",length(DataSet[1,]),"\n"))

sc<- DISCBIO(DataSet)       # The DISCBIO class is the central object storing all information generated throughout the pipeline 

# Estimating a value for the "mintotal" parameters
# As a common practice, mintotal is set to 1000
S1<-summary(colSums(DataSet,na.rm=TRUE))            # It gives an idea about the number of reads across cells
print(S1) 

# Estimating a value for the "minexpr" parameter
S2<-summary(rowMeans(DataSet,na.rm=TRUE))            # It gives an idea about the overall expression of the genes
print(S2)                                                 
minexpr= S2[3]                                       # S2[3] is referring to the median whereas S2[4] is referring to the mean

# Estimating a value for the "minnumber" parameters
minnumber= round(length(DataSet[1,])/10)                             # To be expressed in at 10% of the cells.
print(minnumber)

sc<-Normalizedata(sc, mintotal=1000, minexpr=minexpr, minnumber=minnumber, maxexpr=Inf, downsample=FALSE, dsn=1, rseed=17000) 
sc<-FinalPreprocessing(sc,GeneFlitering="ExpF",export = TRUE)        # The GeneFiltering should be set to "ExpF"

load("SC.RData")           # Loading the "SC" object that has include the data of the k-means clustering 
load("Ndata.RData")        # Loading the "Ndata" object and stored in the @ndata will be used to plot the expression of genes 
load("expdata.RData")      # Loading the "expdata" object and stored in the @expdata will be used to plot the expression of genes 
sc<-SC                     # Storing the data of SC in the sc 
sc@ndata<-Ndata
sc@expdata<-expdata

########## Removing the unneeded objects
rm(Ndata)
rm(expdata)
rm(DataSet)
rm(SC)
                                                    
#sc<- Clustexp(sc,cln=4,quiet=T,clustnr=6,rseed=17000)    
plotGap(sc)                                               ### Plotting gap statistics

outlg<-round(length(sc@fdata[,1]) * 0.05)     # The cell will be considered as an outlier if it has a minimum of 5% of the number of filtered genes as outlier genes.
Outliers<- FindOutliers(sc, K=4, outminc=5,outlg=outlg,plot = TRUE, quiet = FALSE)

# Silhouette plot
options(repr.plot.width=12, repr.plot.height=25)
plotSilhouette(sc,K=4)       # K is the number of clusters

# Jaccard Plot
options(repr.plot.width=10, repr.plot.height=12)
Jaccard(sc,Clustering="K-means", K=4, plot = TRUE)     # Jaccard 

############ Plotting the clusters
plottSNE(sc)

sc<-pseudoTimeOrdering(sc,quiet = TRUE, export = FALSE)
plotOrderTsne(sc)

g='ENSG00000104413'                   #### Plotting the log expression of  ESRP1
plotExptSNE(sc,g)

g='ENSG00000251562'                   #### Plotting the log expression of  MALAT1
plotExptSNE(sc,g)
