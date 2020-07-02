#' @title RPART Decision Tree
#' @description The decision tree analysis is implemented over a training
#'   dataset, which consisted of the DEGs obtained by either SAMseq or the
#'   binomial differential expression.
#' @param data The exact output of the exprmclust function.
#' @param quiet If `TRUE`, suppresses intermediary output
#' @param plot If `FALSE`, suppresses plot output
#' @export
#' @importFrom rpart rpart
#' @importFrom rpart.plot rpart.plot
#' @return Information about the model and, by default, a plot of the decision
#'   tree.
#' @examples
#' data(valuesG1msReduced_treated_K)  # details: DIscBIO:::prepExampleDataset
#' sc <- valuesG1msReduced_treated_K
#' cdiff <- DEGanalysis2clust(
#'     sc, Clustering="K-means", K=3, fdr=.2, name="Name", First="CL1",
#'     Second="CL2", export=FALSE
#' )
#' sigDEG <- cdiff[[1]]
#' DATAforDT <- ClassVectoringDT(
#'     sc, Clustering="K-means", K=3, First="CL1", Second="CL2", sigDEG,
#' )
#' RpartDT(DATAforDT)
RpartDT <- function(data, quiet = FALSE, plot = TRUE) {
    exp.df <- as.data.frame(t(data))
    classVector <- factor(colnames(data))
    model <- rpart(
        classVector ~ .,
        exp.df,
        method = "class",
        minsplit = 1,
        cp = -1
    )
    if (!quiet) print(model)
    if (plot) rpart.plot(model, type = 4, extra = 101)
    return(model)
}
