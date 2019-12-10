#' @title title
#' @export
#' @rdname Clustexp
#' @docType methods
#' @param object object
#' @param clustnr clustnr
#' @param bootnr bootnr
#' @param metric metric
#' @param do.gap do.gap
#' @param SE.method SE.method
#' @param SE.factor SE.factor
#' @param B.gap B.gap
#' @param cln cln
#' @param rseed rseed
setGeneric("Clustexp", function(object, clustnr = 20, bootnr = 50,
                                metric = "pearson", do.gap = TRUE,
                                SE.method = "Tibs2001SEmax", SE.factor = .25,
                                B.gap = 50, cln = 0, rseed = 17000) {
        standardGeneric("Clustexp")
})

#' @rdname Clustexp
setMethod(
    f = "Clustexp",
    signature = "PSCANseq",
    definition = function(object, clustnr, bootnr, metric, do.gap, SE.method,
                          SE.factor, B.gap, cln, rseed) {
        # Validation
        if ( ! is.numeric(clustnr) ) stop("clustnr has to be a positive integer") else if ( round(clustnr) != clustnr | clustnr <= 0 ) stop("clustnr has to be a positive integer")
        if ( ! is.numeric(bootnr) ) stop("bootnr has to be a positive integer") else if ( round(bootnr) != bootnr | bootnr <= 0 ) stop("bootnr has to be a positive integer")
        if ( ! ( metric %in% c( "spearman","pearson","kendall","euclidean","maximum","manhattan","canberra","binary","minkowski") ) ) stop("metric has to be one of the following: spearman, pearson, kendall, euclidean, maximum, manhattan, canberra, binary, minkowski")
        if ( ! ( SE.method %in% c( "firstSEmax","Tibs2001SEmax","globalSEmax","firstmax","globalmax") ) ) stop("SE.method has to be one of the following: firstSEmax, Tibs2001SEmax, globalSEmax, firstmax, globalmax")
        if ( ! is.numeric(SE.factor) ) stop("SE.factor has to be a non-negative integer") else if  ( SE.factor < 0 )  stop("SE.factor has to be a non-negative integer")
        if ( ! ( is.numeric(do.gap) | is.logical(do.gap) ) ) stop( "do.gap has to be logical (TRUE/FALSE)" )
        if ( ! is.numeric(B.gap) ) stop("B.gap has to be a positive integer") else if ( round(B.gap) != B.gap | B.gap <= 0 ) stop("B.gap has to be a positive integer")
        if ( ! is.numeric(cln) ) stop("cln has to be a non-negative integer") else if ( round(cln) != cln | cln < 0 ) stop("cln has to be a non-negative integer")          
        if ( ! is.numeric(rseed) ) stop("rseed has to be numeric")
        if ( !do.gap & cln == 0 ) stop("cln has to be a positive integer or do.gap has to be TRUE")
        
        # Operations
        object@clusterpar <- list(clustnr=clustnr,bootnr=bootnr,metric=metric,do.gap=do.gap,SE.method=SE.method,SE.factor=SE.factor,B.gap=B.gap,cln=cln,rseed=rseed)
        y <- clustfun(object@fdata,clustnr,bootnr,metric,do.gap,SE.method,SE.factor,B.gap,cln,rseed)
        object@kmeans    <- list(kpart=y$clb$result$partition, jaccard=y$clb$bootmean, gap=y$gpr)
        object@distances <- as.matrix( y$di )
        set.seed(111111)
        object@fcol <- sample(rainbow(max(y$clb$result$partition)))
        return(object)
})