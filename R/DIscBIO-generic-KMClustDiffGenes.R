#' @title ClustDiffGenes
#' @description description
#' @param object \code{DISCBIO} class object.
#' @param K A numeric value of the number of clusters.
#' @param fdr A numeric value of the false discovery rate. Default is 0.01.
#' @param export A logical vector that allows writing the final gene list in excel file. Default is TRUE.
#' @param quiet if `TRUE`, suppresses intermediate text output
#' @importFrom dplyr select
#' @importFrom stats pbinom median
#' @rdname KMClustDiffGenes
#' @export
#' @examples 
#' sc <- DISCBIO(valuesG1ms)
#' sc <- Clustexp(sc, cln=3, quiet=TRUE) # K-means clustering
#' KMClustDiffGenes(sc, K=3, fdr=.3)
setGeneric("KMClustDiffGenes", function(object,K,fdr=.01, export=TRUE, quiet=FALSE) standardGeneric("KMClustDiffGenes"))
#' @export
#' @rdname KMClustDiffGenes
setMethod("KMClustDiffGenes",
          signature = "DISCBIO",
          definition = function(object,K,fdr,export=TRUE, quiet=FALSE){
            # Validation
            if (!is.numeric(fdr)) {
                stop("pvalue has to be a number between 0 and 1")
            } else if (fdr < 0 | fdr > 1) {
                stop("fdr has to be a number between 0 and 1")
            }
            if (length(object@kmeans$kpart) == 0) {
				stop("run Clustexp before KMClustDiffGenes")
			}

            cdiff <- list()
            x     <- object@ndata
            y     <- object@expdata[,names(object@ndata)]
            part  <- object@kmeans$kpart
			binompval <- function(p,N,n){
				pval   <- pbinom(n,round(N,0),p,lower.tail=TRUE)
				pval[!is.na(pval) & pval > 0.5] <- 1-pval[!is.na(pval) & pval > 0.5]
				return(pval)
			}
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
                if (quiet) {
                    suppressMessages(
                        G_list <- getBM(
                            filters="ensembl_gene_id",
                            attributes=c("ensembl_gene_id","hgnc_symbol"),
                            values=genes,
                            mart= mart,
                            verbose=!quiet
                        )
                    )
                } else {
                    G_list <- getBM(
                        filters="ensembl_gene_id",
                        attributes=c("ensembl_gene_id","hgnc_symbol"),
                        values=genes,
                        mart= mart,
                        verbose=!quiet
                    )
                }
                Final<-cbind(genes,out)
                Final<-merge(Final,G_list,by.x="genes",by.y="ensembl_gene_id")
                Final<-Final[!duplicated(Final[,8]), ]

                rownames(Final) <- Final[, 1] # FIXME: contains duplicate values
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
                if (export) {
                    write.csv(Up, file = paste0("Up-DEG-cluster",n,".csv"))
                }

                Down<-subset(Final,Final[,7]=="Down")
                Down<-select(Down, "Regulation","genes","pv","mean.all", "mean.cl","fc","p.adj")
                Down[,3]<-rownames(Down)
                Down[,6]<-log2(Down[,6])
                Down[,1]<- Down[,2]
                colnames(Down)<-c("Genes","genes","E.genes","mean.all", "mean.cl","log2.fc","p.adj")
                if (export) {
                    write.csv(Down, file = paste0("Down-DEG-cluster",n,".csv"))
                }
    
                sigDEG<-cbind(DEGsE,DEGsS)
                if (export) {
                    write.csv(sigDEG, file = "binomial-sigDEG.csv")
                }
                
                DEGsTable[n,1]<-paste0("Cluster ",n)
                DEGsTable[n,2]<-"Remaining Clusters"
                DEGsTable[n,3]<-length(Up[,1])
                DEGsTable[n,4]<-paste0("Up-DEG-cluster",n,".csv")
                DEGsTable[n,5]<-length(Down[,1])
                DEGsTable[n,6]<- paste0("Down-DEG-cluster",n,".csv")
            }
            colnames(DEGsTable)<-c("Target Cluster","VS","Up-regulated genes","File name","Low-regulated genes","File name")
            if (export) {
                write.csv(DEGsTable, file = "binomial-DEGsTable.csv")
            }

            return(list(sigDEG,DEGsTable))
                          
          }
          )