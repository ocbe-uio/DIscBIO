#' DISCBIO
#' @title The DISCBIO Class
#' @description The DISCBIO class is the central object storing all information generated throughout the pipeline. 
#' @slot expdata The raw expression data matrix with cells as columns and genes as rows in sparse matrix format. It does not contain ERCC spike-ins.
#' @slot expdataAll The raw expression data matrix with cells as columns and genes as rows in sparse matrix format. It can contain ERCC spike-ins.
#' @slot ndata Data with expression normalized to one for each cell.     
#' @slot fdata Filtered data with expression normalized to one for each cell.   
#' @slot distances A distance matrix.
#' @slot tsne A data.frame with coordinates of two-dimensional tsne layout for the K-means clustering.    
#' @slot background A list storing the polynomial fit for the background model of gene expression variability. It is used for outlier identification. 
#' @slot out A list storing information on outlier cells used for the prediction of rare cell types.  
#' @slot cpart A vector containing the final clustering partition computed by K-means.   
#' @slot fcol A vector contaning the colour scheme for the clusters.    
#' @slot filterpar A list containing the parameters used for cell and gene filterung based on expression.
#' @slot clusterpar A list containing the parameters used for the K-means clustering.
#' @slot outlierpar A list containing the parameters used for outlier identification.
#' @slot kmeans A list containing the results of running the Clustexp() function.    
#' @slot MBclusters A vector containing the final clustering partition computed by Model-based clustering.
#' @slot kordering  A vector containing the Pseudo-time ordering based on k-means clusters.
#' @slot MBordering A vector containing the Pseudo-time ordering based on Model-based clusters.
#' @slot MBtsne A data.frame with coordinates of two-dimensional tsne layout for the Model-based clustering. 
#' @slot noiseF A vector containing the gene list resulted from running the noise filtering.
#' @slot FinalGeneList  A vector containing the final gene list resulted from running the noise filtering or/and the expression filtering.
#' @importFrom methods new validObject
#' @name DISCBIO
#' @rdname DISCBIO
#' @aliases DISCBIO-class, DISCBIO-class
#' @exportClass DISCBIO
#' @examples 
#' class(valuesG1msReduced)
#' G1_reclassified <- DISCBIO(valuesG1msReduced)
#' class(G1_reclassified)
#' str(G1_reclassified, max.level=2) 
#' identical(G1_reclassified@expdataAll, valuesG1msReduced)
#' @export
DISCBIO <- setClass(
    Class = "DISCBIO",
    slots = c(
        expdata         = "data.frame",
        expdataAll      = "data.frame",
        ndata           = "data.frame",
        fdata           = "data.frame", 
        distances       = "matrix",
        tsne            = "data.frame",
        background      = "list",
        out             = "list", 
        cpart           = "vector",
        fcol            = "vector",
        filterpar       = "list",
        clusterpar      = "list", 
        outlierpar      = "list",
        kmeans          = "list",
        MBclusters      = "vector",
        kordering       = "vector",
        MBordering      = "vector",
        MBtsne          = "data.frame",
        noiseF          = "vector",
        FinalGeneList   = "vector"
    )
)

#' validity function for DISCBIO
#'
#' @param object An DISCBIO object.
#' @name DISCBIO
#' @export
setValidity("DISCBIO",
    function(object) {
        msg <- NULL
        if ( ! is.data.frame(object@expdata) ){
            msg <- c(msg, "input data must be data.frame")
        } else if ( nrow(object@expdata) < 2 ){
            msg <- c(msg, "input data must have more than one row")
        } else if ( ncol(object@expdata) < 2 ){
            msg <- c(msg, "input data must have more than one column")
        } else if (sum( apply( is.na(object@expdata),1,sum ) ) > 0 ){
            msg <- c(msg, "NAs are not allowed in input data")
        } else if (sum( apply( object@expdata,1,min ) ) < 0 ){
            msg <- c(msg, "negative values are not allowed in input data")
        }
        if (is.null(msg)) TRUE else msg
    }
)

setMethod("initialize",
    signature = "DISCBIO",
    definition = function(.Object, expdataAll ) {
        .Object@expdataAll <- expdataAll
        shortNames <- substr(rownames(expdataAll), 1, 4)
        geneTypes <- factor(c(ENSG = "ENSG", ERCC = "ERCC" )[shortNames])
        expdata <- expdataAll[which(geneTypes == "ENSG"), ]
        .Object@expdata <- expdata
        .Object@ndata <- expdata
        .Object@fdata <- expdata
        validObject(.Object)
        return(.Object)
    }
)