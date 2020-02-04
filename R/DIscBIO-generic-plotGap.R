#' @title Plotting Gap Statistics
#' @export
#' @rdname plotGap
#' @param object \code{DISCBIO} class object.
#' @examples
#' sc <- DISCBIO(valuesG1ms) # changes signature of data
#' sc <- Clustexp(sc, cln=3) # data must be clustered before plotting
#' plotGap(sc)
setGeneric("plotGap", function(object) standardGeneric("plotGap"))

#' @rdname plotGap
setMethod("plotGap",
          signature = "DISCBIO",
          definition = function(object){
            if ( length(object@kmeans$kpart) == 0 ) stop("run clustexp before plotgap")
            plot(object@kmeans$gap,ylim=c(0.1,0.5),las=1,main="Gap Statistics")
          }
          )		 