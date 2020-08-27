#' @title Networking analysis.
#' @description This function checks the connectivity degree and the betweenness
#'   centrality, which reflect the communication flow in the defined PPI
#'   networks
#' @export
#' @param data Protein-protein interaction data frame resulted from running the
#'   PPI function.
#' @param export if `TRUE`, exports the analysis table as a csv file
#' @param FileName suffix for the file name (if export = TRUE)
#' @importFrom igraph graph.data.frame as_adjacency_matrix distance_table
#'   average.path.length get.adjacency V E mean_distance betweenness
#' @importFrom NetIndices GenInd
#' @return A network analysis table
NetAnalysis <- function(
    data, export = FALSE, FileName = "NetworkAnalysisTable-1"
) {
    if (length(data[, 1]) < 1)
        stop("No Protein-Protein Interactions")
    df <- data[, -c(1, 2)]
    gg <- graph.data.frame(df)
    betweenness <- betweenness(gg)
    betweenness.table <- data.frame(betweenness)
    names <- rownames(betweenness.table)
    rownames(betweenness.table) <- NULL
    degree <- degree(gg)
    degree.table <- data.frame(degree)
    names <- rownames(degree.table)
    rownames(degree.table) <- NULL
    AnalysisTable <- cbind(names, degree.table, betweenness.table)

    if (export) {
        write.csv(AnalysisTable, file = paste0(FileName, ".csv")
        )
    }

    test.graph.adj <- get.adjacency(gg, sparse = FALSE)
    test.graph.properties <- GenInd(test.graph.adj)
    message("Number of nodes: ", test.graph.properties$N)
    message("Number of links: ", test.graph.properties$Ltot)
    message("Link Density: ", test.graph.properties$LD)
    message("The connectance of the graph: ", test.graph.properties$C)
    message("Mean Distences", mean_distance(gg))
    message("Average Path Length", average.path.length(gg), "\n")
    AnalysisTable <-
        AnalysisTable[order(AnalysisTable[, 2], decreasing = TRUE), ]
    return(AnalysisTable)
}
