#' @export
#' @title title
#' @rdname MBClustDiffGenes
setGeneric("MBClustDiffGenes", function(object,fdr=.01) standardGeneric("MBClustDiffGenes"))
#' @export
#' @rdname MBClustDiffGenes
#' @param object object
#' @param fdr fdr
setMethod("MBClustDiffGenes",
          signature = "PSCANseq",
          definition = function(object,fdr){
            if ( ! is.numeric(fdr) ) stop("pvalue has to be a number between 0 and 1") else if (  fdr < 0 | fdr > 1 ) stop("fdr has to be a number between 0 and 1")
            cdiff <- list()
            x     <- object@ndata
            y     <- object@expdata[,names(object@ndata)]
            part  <- object@MBclusters$clusterid
            for ( i in 1:max(part) ){
              if ( sum(part == i) == 0 ) next
              m <- apply(x,1,mean)
              n <- if ( sum(part == i) > 1 ) apply(x[,part == i],1,mean) else x[,part == i]
              no <- if ( sum(part == i) > 1 ) median(apply(y[,part == i],2,sum))/median(apply(x[,part == i],2,sum)) else sum(y[,part == i])/sum(x[,part == i])
              m <- m*no
              n <- n*no
              pv <- binompval(m/sum(m),sum(n),n)
              d <- data.frame(mean.all=m,mean.cl=n,fc=n/m,pv=pv)[order(pv,decreasing=FALSE),]
              #cdiff[[paste("cl",i,sep=".")]] <- d[d$pv < pvalue,]
              cdiff[[i]] <- d[d$pv < fdr,]
            }
            DEGsE<-c()
            DEGsS<-c()
            DEGsTable<-data.frame()

            for (n in 1:K){
                p.adj<- p.adjust(cdiff[[n]][,4],method="bonferroni")
                out<-cbind(cdiff[[n]], p.adj)
                out<-subset(out,out[,5]<fdr)
                Regulation<-c()
                for (i in 1:length(out[,1])){
                    if(out[i,1]>out[i,2]){
                        Regulation[i]="Down"
                    }else{
                        Regulation[i]="Up"
                    }
                }
                out<-cbind(out,Regulation)
                mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
                genes <- rownames(out)
                G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
                Final<-cbind(genes,out)
                Final<-merge(Final,G_list,by.x="genes",by.y="ensembl_gene_id")
                Final<-Final[!duplicated(Final[,8]), ]
                rownames(Final)<-Final[,1]
                Final[,1]<-Final[,8]
                Final<-Final[,-8]
    
                DEGsS<-c(DEGsS,Final[,1])
                DEGsE<-c(DEGsE,as.character(rownames(Final)))
                
                Up<-subset(Final,Final[,7]=="Up")
                Up<-select(Up, "Regulation","genes","pv","mean.all", "mean.cl","fc","p.adj")
                Up[,3]<-rownames(Up)
                Up[,6]<-log2(Up[,6])
                Up[,1]<-Up[,2]
                colnames(Up)<-c("Genes","genes","E.genes","mean.all", "mean.cl","log2.fc","p.adj")
                write.csv(Up, file = paste0("Up-DEG-cluster",n,".csv"))

                Down<-subset(Final,Final[,7]=="Down")
                Down<-select(Down, "Regulation","genes","pv","mean.all", "mean.cl","fc","p.adj")
                Down[,3]<-rownames(Down)
                Down[,6]<-log2(Down[,6])
                Down[,1]<- Down[,2]
                colnames(Down)<-c("Genes","genes","E.genes","mean.all", "mean.cl","log2.fc","p.adj")
                write.csv(Down, file = paste0("Down-DEG-cluster",n,".csv"))
    
    
                sigDEG<-cbind(DEGsE,DEGsS)
                write.csv(sigDEG, file = "binomial-sigDEG.csv")
                
                DEGsTable[n,1]<-paste0("Cluster ",n)
                DEGsTable[n,2]<-"Remaining Clusters"
                DEGsTable[n,3]<-length(Up[,1])
                DEGsTable[n,4]<-paste0("Up-DEG-cluster",n,".csv")
                DEGsTable[n,5]<-length(Down[,1])
                DEGsTable[n,6]<- paste0("Down-DEG-cluster",n,".csv")
            }
            colnames(DEGsTable)<-c("Target Cluster","VS","Up-regulated genes","File name","Low-regulated genes","File name")
            write.csv(DEGsTable, file = "binomial-DEGsTable.csv")

            return(list(sigDEG,DEGsTable))
                          
          }
          )