#' @title Networking analysis.
#' @description This function checks the connectivity degree and the betweenness centrality, which reflect the communication flow in the defined PPI networks
#' @export
#' @param data Protein-protein interaction data frame resulted from running the PPI function. 
#' @param export if `TRUE`, exports the analysis table as a csv file
#' @param FileName suffix for the file name (if export = TRUE)
#' @importFrom igraph graph.data.frame as_adjacency_matrix distance_table average.path.length get.adjacency V E mean_distance betweenness
#' @importFrom NetIndices GenInd
#' @examples 
#' \dontrun{
#' sc <- DISCBIO(valuesG1msReduced)
#' sc <- NoiseFiltering(sc)
#' sc <- Normalizedata(
#'     sc, mintotal=1000, minexpr=0, minnumber=0, maxexpr=Inf, downsample=FALSE,
#'     dsn=1, rseed=17000
#' )
#' sc <- FinalPreprocessing(sc, GeneFlitering="NoiseF")
#' sc <- Clustexp(sc, cln=3) # K-means clustering
#' sc <- comptSNE(sc, rseed=15555)
#' dff <- DEGanalysis2clust(sc, Clustering="K-means", K=3, fdr=0.1, name="Name")
#' DEGs <- dff[[2]][1, 6]
#' data <- read.csv(file=paste0(DEGs),head=TRUE,sep=",")
#' data <- data[,3]
#' FileName <- paste0(DEGs)
#' ppi <- PPI(data, FileName)
#' NetAnalysis(ppi)
#' }
NetAnalysis<-function(data, export=TRUE, FileName="1"){
    if ( length(data[,1])<1 ) stop( "No Protein-Protein Interactions" )
    df<-data[,-c(1,2)]
    gg <- graph.data.frame(df)
    betweenness<-betweenness(gg)
    betweenness.table<- data.frame(betweenness)
    names <- rownames(betweenness.table)
    rownames(betweenness.table) <- NULL
    degree<-degree(gg)
    degree.table<- data.frame(degree)
    names <- rownames(degree.table)
    rownames(degree.table) <- NULL
    AnalysisTable <- cbind(names,degree.table,betweenness.table)

    if (export) {
        write.csv(
            AnalysisTable,
            file = paste0("NetworkAnalysisTable-",FileName,".csv")
        )
    }

    test.graph.adj<-get.adjacency(gg,sparse=F)
    test.graph.properties<-GenInd(test.graph.adj)
    cat("Number of nodes: ",test.graph.properties$N,"\n")
    V(gg)
    cat("Number of links: ",test.graph.properties$Ltot,"\n")
    E(gg)
    cat("Link Density: ",test.graph.properties$LD,"\n")
    cat("The connectance of the graph: ",test.graph.properties$C,"\n")
    cat("Mean Distences",mean_distance(gg),"\n")
    cat("Average Path Length",average.path.length(gg),"\n","\n")
	AnalysisTable<-AnalysisTable[order(AnalysisTable[,2],decreasing=T),]
    return(AnalysisTable)    
}


