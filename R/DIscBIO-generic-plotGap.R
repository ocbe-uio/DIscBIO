#' @title title
#' @export
#' @rdname plotGap
#' @param object object
setGeneric("plotGap", function(object) standardGeneric("plotGap"))

#' @rdname plotGap
setMethod("plotGap",
          signature = "PSCANseq",
          definition = function(object){
            if ( length(object@kmeans$kpart) == 0 ) stop("run clustexp before plotgap")
            plot(object@kmeans$gap,ylim=c(0.1,0.5),las=1,main="Gap Statistics")
          }
          )