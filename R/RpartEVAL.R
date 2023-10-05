#' @title Evaluating the performance of the RPART Decision Tree.
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
#' @return Performance statistics of the model
RpartEVAL <- function(data, num.folds = 10, First = "CL1", Second = "CL2",
                      quiet = FALSE) {
  exp.imput.df <- as.data.frame(t(data))
  num.instances <- nrow(exp.imput.df)
  indices <- 1:num.instances
  classVector <- factor(colnames(data))
  cv.segments <- split(
    sample(indices), rep(1:num.folds, length = num.instances)
  )
  Rpart.performance <- cross.val(
    exp.imput.df, classVector, cv.segments, Rpart.performance, "rpart", quiet
  )
  if (!quiet) print(Rpart.performance)
  Rpart.confusion.matrix <- matrix(Rpart.performance, nrow = 2)
  rownames(Rpart.confusion.matrix) <- c(
    paste0("Predicted", First), paste0("Predicted", Second)
  )
  colnames(Rpart.confusion.matrix) <- c(First, Second)
  if (!quiet) print(Rpart.confusion.matrix)
  Rpart.sn <- round(SN(Rpart.confusion.matrix), digits = 2)
  Rpart.sp <- round(SP(Rpart.confusion.matrix), digits = 2)
  Rpart.acc <- round(ACC(Rpart.confusion.matrix), digits = 2)
  Rpart.mcc <- round(MCC(Rpart.confusion.matrix), digits = 2)

  if (!quiet) {
    message(
      "Rpart SN: ", Rpart.sn, "\n",
      "Rpart SP: ", Rpart.sp, "\n",
      "Rpart ACC: ", Rpart.acc, "\n",
      "Rpart MCC: ", Rpart.mcc, "\n"
    )
  }
  return(Rpart.performance)
}
