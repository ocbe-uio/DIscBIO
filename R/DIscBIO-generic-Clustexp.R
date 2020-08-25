#' @title Clustering of single-cell transcriptome data
#' @description This functions performs the initial clustering of the RaceID
#'   algorithm.
#' @docType methods
#' @param object \code{DISCBIO} class object.
#' @param clustnr Maximum number of clusters for the derivation of the cluster
#'   number by the saturation of mean within-cluster-dispersion. Default is 20.
#' @param bootnr A numeric value of booststrapping runs for \code{clusterboot}.
#'   Default is 50.
#' @param metric Is the method to transform the input data to a distance object.
#'   Metric has to be one of the following: ["spearman", "pearson", "kendall",
#'   "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"].
#' @param do.gap A logical vector that allows generating the number of clusters
#'   based on the gap statistics. Default is TRUE.
#' @param SE.method The SE.method determines the first local maximum of the gap
#'   statistics. The SE.method has to be one of the following:["firstSEmax",
#'   "Tibs2001SEmax", "globalSEmax", "firstmax", "globalmax"]. Default is
#'   "Tibs2001SEmax"
#' @param SE.factor A numeric value of the fraction of the standard deviation by
#'   which the local maximum is required to differ from the neighboring points
#'   it is compared to. Default is 0.25.
#' @param B.gap Number of bootstrap runs for the calculation of the gap
#'   statistics. Default is 50
#' @param cln Number of clusters to be used. Default is \code{NULL} and the
#'   cluster number is inferred by the saturation criterion.
#' @param rseed Random integer to enforce reproducible clustering results.
#' @param quiet if `TRUE`, intermediate output is suppressed
#' @importFrom stats as.dist cor kmeans
#' @importFrom cluster clusGap maxSE
#' @importFrom graphics pairs
#' @importFrom methods is
#' @return The DISCBIO-class object input with the cpart slot filled.
#' @examples
#' sc <- DISCBIO(valuesG1msTest) # changes signature of data
#' sc <- Clustexp(sc, cln=2)
setGeneric("Clustexp", function(object, clustnr = 3, bootnr = 50,
    metric = "pearson", do.gap = TRUE, SE.method = "Tibs2001SEmax",
    SE.factor = .25, B.gap = 50, cln = 0, rseed = NULL, quiet = FALSE)
    {
        standardGeneric("Clustexp")
    }
)

#' @export
#' @rdname Clustexp
setMethod(
    f = "Clustexp",
    signature = "DISCBIO",
    definition = function(object, clustnr, bootnr, metric, do.gap, SE.method,
        SE.factor, B.gap, cln, rseed, quiet)
    {
        if (!is.numeric(clustnr))
            stop("clustnr has to be a positive integer")
        else if (round(clustnr) != clustnr | clustnr <= 0)
            stop("clustnr has to be a positive integer")
        if (!is.numeric(bootnr))
            stop("bootnr has to be a positive integer")
        else if (round(bootnr) != bootnr | bootnr <= 0)
            stop("bootnr has to be a positive integer")
        if (!(
            metric %in% c(
                "spearman", "pearson", "kendall", "euclidean", "maximum",
                "manhattan", "canberra", "binary", "minkowski"
            )
        ))
            stop(
                "metric has to be one of the following: spearman, ",
                "pearson, kendall, euclidean, maximum, manhattan, ",
                "canberra, binary, minkowski"
            )
        if (!(
            SE.method %in% c(
                "firstSEmax", "Tibs2001SEmax", "globalSEmax", "firstmax",
                "globalmax"
            )
        ))
            stop(
                "SE.method has to be one of the following: ",
                "firstSEmax, Tibs2001SEmax, globalSEmax, ",
                "firstmax, globalmax"
            )
        if (!is.numeric(SE.factor))
            stop("SE.factor has to be a non-negative integer")
        else if (SE.factor < 0)
            stop("SE.factor has to be a non-negative integer")
        if (!(is.numeric(do.gap) | is.logical(do.gap)))
            stop("do.gap has to be logical (TRUE/FALSE)")
        if (!is.numeric(B.gap))
            stop("B.gap has to be a positive integer")
        else if (round(B.gap) != B.gap | B.gap <= 0)
            stop("B.gap has to be a positive integer")
        if (!is.numeric(cln))
            stop("cln has to be a non-negative integer")
        else if (round(cln) != cln | cln < 0)
            stop("cln has to be a non-negative integer")
        if (!is.null(rseed) & !is.numeric(rseed))
            stop("rseed has to be numeric or NULL")
        if (!do.gap & cln == 0)
            stop("cln has to be a positive integer or do.gap has to be TRUE")

        # Operations
        object@clusterpar <-
            list(
                clustnr = clustnr,
                bootnr = bootnr,
                metric = metric,
                do.gap = do.gap,
                SE.method = SE.method,
                SE.factor = SE.factor,
                B.gap = B.gap,
                cln = cln,
                rseed = rseed
            )
        y <- clustfun(
            object@fdata,
            clustnr,
            bootnr,
            metric,
            do.gap,
            SE.method,
            SE.factor,
            B.gap,
            cln,
            rseed = rseed,
            quiet = quiet
        )
        object@kmeans <- list(
            kpart = y$clb$result$partition,
            jaccard = y$clb$bootmean,
            gap = y$gpr
        )
        object@distances <- as.matrix(y$di)
        object@fcol <- rainbow(max(y$clb$result$partition))
        object@cpart <- object@kmeans$kpart
        return(object)
    }
)