#' @title Normalizing and filtering
#' @description This function allows filtering of genes and cells to be used in the downstream analysis.
#' @param object \code{PSCANseq} class object.
#' @param mintotal minimum total transcript number required. Cells with less than \code{mintotal} transcripts are filtered out. Default is 1000.
#' @param minexpr minimum required transcript count of a gene in at least \code{minnumber} cells. All other genes are filtered out. Default is 0.
#' @param minnumber minimum number of cells that are expressing each gene at minexpr transcripts. Default is 0.
#' @param maxexpr maximum allowed transcript count of a gene in at least a single cell after normalization or downsampling. All other genes are filtered out. Default is Inf.
#' @param downsample A logical vector. Default is FALSE. If downsample is set to TRUE, then transcript counts are downsampled to mintotal transcripts per cell, instead of 
#' the normalization. Downsampled versions of the transcript count data are averaged across dsn samples
#' @param dsn A numeric value of the number of samples to be used to average the downsampled versions of the transcript count data. Default is 1 which means that sampling
#' noise should be comparable across cells. For high numbers of dsn the data will become similar to the median normalization.
#' @param rseed Integer number. Random seed to enforce reproducible clustering results. Default is 17000.
#' @include DIscBIO-classes.R
setGeneric("Normalizedata", function(object, mintotal=1000, minexpr=0, minnumber=0, maxexpr=Inf, downsample=FALSE, dsn=1, rseed=17000) standardGeneric("Normalizedata"))

#' @export
#' @rdname Normalizedata
setMethod("Normalizedata",
          signature = "PSCANseq",
          definition = function(object,mintotal,minexpr,minnumber,maxexpr,downsample,dsn,rseed) {
            if ( ! is.numeric(mintotal) ) stop( "mintotal has to be a positive number" ) else if ( mintotal <= 0 ) stop( "mintotal has to be a positive number" )
            if ( ! is.numeric(minexpr) ) stop( "minexpr has to be a non-negative number" ) else if ( minexpr < 0 ) stop( "minexpr has to be a non-negative number" )
            if ( ! is.numeric(minnumber) ) stop( "minnumber has to be a non-negative integer number" ) else if ( round(minnumber) != minnumber | minnumber < 0 ) stop( "minnumber has to be a non-negative integer number" )
            if ( ! ( is.numeric(downsample) | is.logical(downsample) ) ) stop( "downsample has to be logical (TRUE/FALSE)" )
            if ( ! is.numeric(dsn) ) stop( "dsn has to be a positive integer number" ) else if ( round(dsn) != dsn | dsn <= 0 ) stop( "dsn has to be a positive integer number" )
            object@filterpar <- list(mintotal=mintotal, minexpr=minexpr, minnumber=minnumber, maxexpr=maxexpr, downsample=downsample, dsn=dsn)
            object@ndata <- object@expdata[,apply(object@expdata,2,sum,na.rm=TRUE) >= mintotal]
                
            if ( downsample ){
                downsample <- function(x,n,dsn){
				x <- round( x[,apply(x,2,sum,na.rm=TRUE) >= n], 0)
				nn <- min( apply(x,2,sum) )
				for ( j in 1:dsn ){
					z  <- data.frame(GENEID=rownames(x))
					rownames(z) <- rownames(x)
					initv <- rep(0,nrow(z))
					for ( i in 1:dim(x)[2] ){
						y <- aggregate(rep(1,nn),list(sample(rep(rownames(x),x[,i]),nn)),sum)
						na <- names(x)[i]
						names(y) <- c("GENEID",na)
						rownames(y) <- y$GENEID
						z[,na] <- initv
						k <- intersect(rownames(z),y$GENEID)
						z[k,na] <- y[k,na]
						z[is.na(z[,na]),na] <- 0
					}
					rownames(z) <- as.vector(z$GENEID)
					ds <- if ( j == 1 ) z[,-1] else ds + z[,-1]
				}
				ds <- ds/dsn + .1
				return(ds)
                }
                set.seed(rseed)
                object@ndata <- downsample(object@expdata,n=mintotal,dsn=dsn)
            }else{
              x <- object@ndata
              object@ndata <- as.data.frame( t(t(x)/apply(x,2,sum))*median(apply(x,2,sum,na.rm=TRUE)) + .1 )
            }
            x <- object@ndata
            object@fdata <- x[apply(x>=minexpr,1,sum,na.rm=TRUE) >= minnumber,]
            x <- object@fdata
            object@fdata <- x[apply(x,1,max,na.rm=TRUE) < maxexpr,]
            return(object)
    }
    )