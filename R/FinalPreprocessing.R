#' @title Final Preprocessing
#' @description This function generates the final filtered normalized dataset.
#' @param object \code{PSCANseq} class object.
#' @param GeneFlitering GeneFlitering has to be one of the followings: ["NoiseF","ExpF"]. Default is "NoiseF"
#' @param export A logical vector that allows writing the final gene list in excel file. Default is TRUE. 
#' @export
#' @rdname FinalPreprocessing

setGeneric("FinalPreprocessing", function(object,GeneFlitering="NoiseF",export = TRUE) standardGeneric("FinalPreprocessing"))

setMethod("FinalPreprocessing",
          signature = "PSCANseq",
          definition = function(object,GeneFlitering,export = TRUE){
			if (GeneFlitering == "NoiseF") {
				if ( length(object@noiseF) < 1 ) stop("run NoiseFiltering before running FinalPreprocessing")
				if ( nrow(object@fdata) < 1 ) stop("run Normalizedata before running FinalPreprocessing")
				gene_list<-object@noiseF
				gene_names<-rownames(object@fdata)
				idx_genes <- is.element(gene_names, gene_list)
				gene_names2 <- gene_names[idx_genes]
				filteredDataset<- object@fdata[gene_names2,]
				object@fdata<-filteredDataset
				object@FinalGeneList<-rownames(filteredDataset)
				cat("The gene filtering method= Noise filtering","\n","\n") 
				cat("The Filtered Normalized dataset contains:","\n","Genes:",length(filteredDataset[,1]),"\n","cells:",length(filteredDataset[1,]),"\n","\n")
				if (export) {
					cat("The Filtered Normalized dataset was saved as: filteredDataset.Rdata\n")
					save(filteredDataset,file="filteredDataset.Rdata")
				}
			}	

			if (GeneFlitering == "ExpF") {
				if ( nrow(object@fdata) < 1 ) stop("run Normalizedata before running FinalPreprocessing")
				filteredDataset<-object@fdata
				object@FinalGeneList<-rownames(filteredDataset)
				cat("The gene filtering method= Expression filtering","\n","\n") 
				cat("The Filtered Normalized dataset contains:","\n","Genes:",length(filteredDataset[,1]),"\n","cells:",length(filteredDataset[1,]),"\n","\n")
				if (export) {
					cat("The Filtered Normalized dataset was saved as: filteredDataset.Rdata\n")
					save(filteredDataset,file="filteredDataset.Rdata")
				}
			}
            return(object)
})
