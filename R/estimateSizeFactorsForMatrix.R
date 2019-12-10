#' @title title
#' @description description 
#' @export
#' @param counts counts
#' @param locfunc locfunc
#' @importFrom stats median
estimateSizeFactorsForMatrix = function (counts, locfunc = median){
	loggeomeans <- rowMeans(log(counts))
	apply(counts, 2, function(cnts) exp(locfunc((log(cnts) - loggeomeans)[is.finite(loggeomeans)])))
}