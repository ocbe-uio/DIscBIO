#' @title J48 Decision Tree
#' @description The decision tree analysis is implemented over a training dataset, which consisted of the DEGs obtained by either SAMseq or the binomial differential expression.  
#' @export
#' @param data A data frame resulted from running the function ClassVectoringDT.
#' @param quiet If `TRUE`, suppresses intermediary output
#' @param plot If `FALSE`, suppresses plot output
#' @importFrom RWeka J48
#' @importFrom graphics plot
#' @importFrom partykit as.party
#' @importFrom grid gpar
#' @return Information about the J48 model and, by default, a plot of the decision tree.
#' @examples
#' sc <- DISCBIO(valuesG1msReduced)
#' sc <- NoiseFiltering(sc, percentile=0.9, CV=0.2, export=FALSE)
#' sc <- Normalizedata(
#'     sc, mintotal=1000, minexpr=0, minnumber=0, maxexpr=Inf, downsample=FALSE,
#'     dsn=1, rseed=17000
#' )
#' sc <- FinalPreprocessing(sc, GeneFlitering="NoiseF", export=FALSE)
#' sc <- Clustexp(sc, cln=3) # K-means clustering
#' sc <- comptSNE(sc, rseed=15555)
#' cdiff <- DEGanalysis2clust(
#'     sc, Clustering="K-means", K=3, fdr=.2, name="Name", First="CL1",
#'     Second="CL2", export=FALSE
#' )
#' sigDEG <- cdiff[[1]]
#' DATAforDT <- ClassVectoringDT(
#'     sc, Clustering="K-means", K=3, First="CL1", Second="CL2", sigDEG,
#' )
#' J48DT(DATAforDT)
J48DT<-function(data, quiet = FALSE, plot = TRUE){
		msg <- NULL
		if ( ! is.data.frame(data) ){
            msg <- c(msg, "input data must be data.frame")
        }else if ( nrow(data) < 2 ){
            msg <- c(msg, "input data must have more than one row")
        }else if ( ncol(data) < 2 ){
            msg <- c(msg, "input data must have more than one column")
        }else if (sum( apply( is.na(data),1,sum ) ) > 0 ){
            msg <- c(msg, "NAs are not allowed in input data")
        }else if (sum( apply( data,1,min ) ) < 0 ){
            msg <- c(msg, "negative values are not allowed in input data")
        }
        if (is.null(msg)) TRUE
        else msg
		
		exp.df<-as.data.frame(t(data))
		classVector<- factor(colnames(data))
		j48.model<-J48(classVector~.,exp.df)
		if (!quiet) print(j48.model)
		if (plot) {
            plot(
                as.party(j48.model),
                gp = gpar(
                    cex = 0.65, col = "black", lty = "solid", lwd = 1.5, 
                    fontsize = 12
                )
            )
        }
		return(j48.model)
}