#' @title J48 Decision Tree
#' @description The decision tree analysis is implemented over a training
#'   dataset, which consisted of the DEGs obtained by either SAMseq or the
#'   binomial differential expression.
#' @export
#' @param data A data frame resulted from running the function ClassVectoringDT.
#' @param quiet If `TRUE`, suppresses intermediary output
#' @param plot If `FALSE`, suppresses plot output
#' @importFrom RWeka J48
#' @importFrom graphics plot
#' @return Information about the J48 model and, by default, a plot of the
#'   decision tree.
J48DT <- function(data, quiet = FALSE, plot = TRUE) {
  msg <- NULL
  if (!is.data.frame(data)) {
    msg <- c(msg, "input data must be data.frame")
  } else if (nrow(data) < 2) {
    msg <- c(msg, "input data must have more than one row")
  } else if (ncol(data) < 2) {
    msg <- c(msg, "input data must have more than one column")
  } else if (sum(apply(is.na(data), 1, sum)) > 0) {
    msg <- c(msg, "NAs are not allowed in input data")
  } else if (sum(apply(data, 1, min)) < 0) {
    msg <- c(msg, "negative values are not allowed in input data")
  }
  if (is.null(msg)) TRUE else msg

  exp.df <- as.data.frame(t(data))
  exp.df$classVector <- factor(colnames(data))
  j48.model <- J48(classVector ~ ., exp.df)
  if (!quiet) print(j48.model)
  if (plot) plot(j48.model)
  return(j48.model)
}
