#' @title Normalizing and filtering
#' @description This function allows filtering of genes and cells to be used in
#'   the downstream analysis.
#' @param object \code{DISCBIO} class object.
#' @param mintotal minimum total transcript number required. Cells with less
#'   than \code{mintotal} transcripts are filtered out. Default is 1000.
#' @param minexpr minimum required transcript count of a gene in at least
#'   \code{minnumber} cells. All other genes are filtered out. Default is 0.
#' @param minnumber minimum number of cells that are expressing each gene at
#'   minexpr transcripts. Default is 0.
#' @param maxexpr maximum allowed transcript count of a gene in at least a
#'   single cell after normalization or downsampling. All other genes are
#'   filtered out. Default is Inf.
#' @param downsample A logical vector. Default is FALSE. If downsample is set to
#'   TRUE, then transcript counts are downsampled to mintotal transcripts per
#'   cell, instead of the normalization. Downsampled versions of the transcript
#'   count data are averaged across dsn samples
#' @param dsn A numeric value of the number of samples to be used to average the
#'   downsampled versions of the transcript count data. Default is 1 which means
#'   that sampling noise should be comparable across cells. For high numbers of
#'   dsn the data will become similar to the median normalization.
#' @param rseed Random integer to enforce reproducible clustering.
#'   results
#' @include DIscBIO-classes.R
#' @return The DISCBIO-class object input with the ndata and fdata slots filled.
#' @examples
#' sc <- DISCBIO(valuesG1msTest) # changes signature of data
#'
#' # In this case this function is used to normalize the reads
#' sc_normal <- Normalizedata(
#'     sc, mintotal=1000, minexpr=0, minnumber=0, maxexpr=Inf, downsample=FALSE,
#'     dsn=1, rseed=17000
#' )
#' summary(sc_normal@fdata)
#'
setGeneric(
    "Normalizedata",
    function(
        object, mintotal = 1000, minexpr = 0, minnumber = 0, maxexpr = Inf,
        downsample = FALSE, dsn = 1, rseed = NULL
    )
    standardGeneric("Normalizedata")
)

#' @export
#' @rdname Normalizedata
setMethod(
    "Normalizedata",
    signature = "DISCBIO",
    definition = function(
        object, mintotal, minexpr, minnumber, maxexpr, downsample, dsn, rseed) {
        # Validation
        if (!is.numeric(mintotal))
            stop("mintotal has to be a positive number")
        else if (mintotal <= 0)
            stop("mintotal has to be a positive number")
        if (!is.numeric(minexpr))
            stop("minexpr has to be a non-negative number")
        else if (minexpr < 0)
            stop("minexpr has to be a non-negative number")
        if (!is.numeric(minnumber))
            stop("minnumber has to be a non-negative integer number")
        else if (round(minnumber) != minnumber | minnumber < 0)
            stop("minnumber has to be a non-negative integer number")
        if (!(is.numeric(downsample) | is.logical(downsample)))
            stop("downsample has to be logical (TRUE/FALSE)")
        if (!is.numeric(dsn))
            stop("dsn has to be a positive integer number")
        else if (round(dsn) != dsn | dsn <= 0)
            stop("dsn has to be a positive integer number")
        object@filterpar <- list(
            mintotal = mintotal,
            minexpr = minexpr,
            minnumber = minnumber,
            maxexpr = maxexpr,
            downsample = downsample,
            dsn = dsn
        )
        cols <- apply(object@expdata, 2, sum, na.rm = TRUE) >= mintotal
        object@ndata <- object@expdata[, cols]
        if (downsample) {
            set.seed(rseed)
            object@ndata <- downsample(object@expdata, n = mintotal, dsn = dsn)
        } else{
            x <- object@ndata
            object@ndata <- as.data.frame(t(t(x) / apply(x, 2, sum)) *
                median(apply(x, 2, sum, na.rm = TRUE)) + .1)
        }
        x <- object@ndata
        object@fdata <-
            x[apply(x >= minexpr, 1, sum, na.rm = TRUE) >= minnumber, ]
        x <- object@fdata
        object@fdata <- x[apply(x, 1, max, na.rm = TRUE) < maxexpr, ]
        return(object)
    }
)