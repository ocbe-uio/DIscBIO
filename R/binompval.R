#' @title title
#' @description description 
#' @export
#' @param p p
#' @param N N
#' @param n n
#' @importFrom stats pbinom
binompval <- function(p,N,n){
  pval   <- pbinom(n,round(N,0),p,lower.tail=TRUE)
  pval[!is.na(pval) & pval > 0.5] <- 1-pval[!is.na(pval) & pval > 0.5]
  return(pval)
}