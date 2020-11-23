library(DIscBIO)
library(partykit)
library(enrichR)

load("SC.RData")           # Loading the "SC" object that has include the data of the k-means clustering 
load("Ndata.RData")        # Loading the "Ndata" object and stored in the @ndata will be used to plot the expression of genes 
load("expdata.RData")      # Loading the "expdata" object and stored in the @expdata will be used to plot the expression of genes 
sc<-SC                     # Storing the data of SC in the sc 
sc@ndata<-Ndata
sc@expdata<-expdata

########## Removing the unneeded objects
rm(Ndata)
rm(expdata)
rm(SC)
####### differential expression analysis between cluster 1 and cluster 4 using FDR of 0.05
cdiff<-DEGanalysis2clust(sc,Clustering="K-means",K=4,fdr=0.05,name="Name",First="CL1",Second="CL4",export = TRUE,quiet=T)     

#### To show the result table
head(cdiff[[1]])                  # The first component 
head(cdiff[[2]])                  # The second component 

cdiffBinomial<-ClustDiffGenes(sc,K=4,export = T,fdr=.01,quiet=T)    ########## Binomial differential expression analysis
#### To show the result table
head(cdiffBinomial[[1]])                  # The first component 
head(cdiffBinomial[[2]])                  # The second component

name<-cdiffBinomial[[2]][1,6]    ############ Selecting the DEGs' ############## Down-DEG-cluster1.csv
U<-read.csv(file=paste0(name),head=TRUE,sep=",")
Vplot<-VolcanoPlot(U,value=0.0001,name=name,FS=0.7,fc=0.75)

###################### Finding biomarker genes between cluster 1 and cluster 4
First="CL1"
Second="CL4"
load("DATAforDT.RData")

j48dt<-J48DT(DATAforDT)           #J48 Decision Tree
summary(j48dt) 
rm(j48dt)

rpartDT<-RpartDT(DATAforDT)
rm(rpartDT)

DEGs="All_DEGs"
FileName=paste0(DEGs)

data<-cdiffBinomial[[1]] [1:200,2]       # DEGs gene list from Binomial analysis (taking only the firat 200 genes)

ppi<-PPI(data,FileName)

networking<-NetAnalysis(ppi)
networking                            ##### In case the Examine response components = 200 and an error "linkmat[i, ]" appeared, that means there are no PPI.

data=networking[1:25,1]              # plotting the network of the top 25 highly connected genes 
network<-Networking(data,FileName,plot_width = 25, plot_height = 10)


############ Selecting the DEGs' table  ##############
DEGs=cdiff[[2]][1,4]             # Up-regulated genes in cluster 4 (from SAMseq)
FileName=paste0(DEGs)

data<-read.csv(file=paste0(DEGs),head=TRUE,sep=",")
data<-data[,3]
ppi<-PPI(data,FileName)
networking<-NetAnalysis(ppi)
networking                            ##### In case the Examine response components = 200 and an error "linkmat[i, ]" appeared, that means there are no PPI.

data=networking[1:25,1]              # plotting the network of the top 25 highly connected genes 
network<-Networking(data,FileName,plot_width = 25, plot_height = 10)

dbs <- listEnrichrDbs()
head(dbs)
#print(dbs)

############ Selecting the DEGs' table  ##############
DEGs=cdiff[[2]][1,4]             # Up-regulated genes in cluster 4 (from SAMseq)
FileName=paste0(DEGs)

data<-read.csv(file=paste0(DEGs),head=TRUE,sep=",")
data<-as.character(data[,3])

dbs <- c("KEGG_2019_Human","GO_Biological_Process_2015")
enriched <- enrichr(data, dbs)
KEGG_2019_Human<-enriched[[1]][,c(1,2,3,9)]
GO_Biological_Process_2015<-enriched[[2]][,c(1,2,3,9)]

GEA<-rbind(KEGG_2019_Human,GO_Biological_Process_2015)
GEA


