#' @title title
#' @description description 
#' @export
#' @param object object
#' @importFrom igraph graph.data.frame as_adjacency_matrix distance_table average.path.length get.adjacency V E mean_distance
#' @importFrom NetIndices GenInd
NetAnalysis<-function(object){
    if ( length(object[,1])<1 ) stop( "No Protein-Protein Interactions" )
    df<-object[,-c(1,2)]
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
    write.csv(AnalysisTable, file = paste0("NetworkAnalysisTable-",FileName,".csv"))

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
    return(AnalysisTable)    
}