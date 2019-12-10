#' @title title
#' @description description 
#' @export
#' @param object object
#' @param Cluster_ID Cluster_ID
#' @param K K
#' @param First First
#' @param Second Second
#' @param sigDEG sigDEG
#' @importFrom biomaRt useDataset useMart getBM
ClassVectoringDT<- function(object,Cluster_ID,K,First,Second,sigDEG){
    SC<-PSCANseq(object)
    SC<- Normalizedata(SC)
    DatasetForDT<- SC@fdata
    Nam <-colnames(DatasetForDT)
    num<-c(1:K)
    num1<- paste("CL", num, sep="")
    for (n in num){
        Nam<-ifelse((Cluster_ID==n),num1[n],Nam)
    }
    colnames(DatasetForDT)<-Nam

    #cat("The dataset is ready for decision tree analysis","\n")
    #dim(DatasetForDT)
    sg1 <- DatasetForDT[,which(colnames(DatasetForDT)==First | colnames(DatasetForDT)==Second)]
    dim(sg1)
    # Creating a dataset that includes only the DEGs
    gene_list<-sigDEG[,1]
    gene_names<-rownames(DatasetForDT)
    idx_genes <- is.element(gene_names, gene_list)
    gene_names2 <- gene_names[idx_genes]
    DEGsfilteredDataset<- sg1[gene_names2,]
    cat("The DEGs filtered normalized dataset contains:","\n","Genes:",length(DEGsfilteredDataset[,1]),"\n","cells:",length(DEGsfilteredDataset[1,]))
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    
    genes <- rownames(DEGsfilteredDataset)
    G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
    
    DATAforDT<-cbind(genes,DEGsfilteredDataset)
    DATAforDT<-merge(DATAforDT,G_list,by.x="genes",by.y="ensembl_gene_id")
    DATAforDT
    DATAforDT[,1]<-DATAforDT[,length(DATAforDT[1,])]
    DATAforDT<-DATAforDT[!duplicated(DATAforDT[,1]), ]

    rownames(DATAforDT)<-DATAforDT[,1]
    DATAforDT<-DATAforDT[,c(-1,-length(DATAforDT[1,]))]
    sg <- factor(gsub(paste0("(",First,"|",Second,").*"), "\\1", colnames(DATAforDT)), levels = c(paste0(First), paste0(Second)))
    sg<-sg[!is.na(sg)]
    colnames(DATAforDT)<-sg
    return(DATAforDT)
}