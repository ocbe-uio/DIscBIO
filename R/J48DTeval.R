#' @title Evaluating the performance of the J48 decision tree.
#' @description This function evaluates the performance of the generated trees
#'   for error estimation by ten-fold cross validation assessment.
#' @export
#' @param data The resulted data from running the function J48DT.
#' @param num.folds A numeric value of the number of folds for the cross
#'   validation assessment. Default is 10.
#' @param First A string vector showing the first target cluster.  Default is
#'   "CL1"
#' @param Second A string vector showing the second target cluster.  Default is
#'   "CL2"
#' @param quiet If `TRUE`, suppresses intermediary output
#' @importFrom stats predict
#' @return Statistics about the J48 model
J48DTeval <- function(
  data, num.folds = 10, First = "CL1", Second = "CL2", quiet = FALSE
) {
  exp.imput.df <- as.data.frame(t(data))
  num.instances <- nrow(exp.imput.df)
  indices <- 1:num.instances
  classVector <- factor(colnames(data))
  cv.segments <- split(
    sample(indices), rep(1:num.folds, length = num.instances)
  )
  j48.performance <- cross.val(
    exp.imput.df, classVector, cv.segments, j48.performance, "J48", quiet
  )
  if (!quiet) print(j48.performance)

  j48.confusion.matrix <- matrix(j48.performance, nrow = 2)
  rownames(j48.confusion.matrix) <- c(
    paste0("Predicted", First), paste0("Predicted", Second)
  )
  colnames(j48.confusion.matrix) <- c(First, Second)
  if (!quiet) print(j48.confusion.matrix)
  j48.sn <- round(SN(j48.confusion.matrix), digits = 2)
  j48.sp <- round(SP(j48.confusion.matrix), digits = 2)
  j48.acc <- round(ACC(j48.confusion.matrix), digits = 2)
  j48.mcc <- round(MCC(j48.confusion.matrix), digits = 2)

  if (!quiet) {
    message(
      "J48 SN: ", j48.sn, "\n",
      "J48 SP: ", j48.sp, "\n",
      "J48 ACC: ", j48.acc, "\n",
      "J48 MCC: ", j48.mcc, "\n"
    )
  }
  return(j48.performance)
}
