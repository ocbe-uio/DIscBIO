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
        cross.val <- function(exp.df, class.vec, segments, performance,
                                class.algo) {
            #Start cross validation loop
            class1 <- levels(class.vec)[1]
            for (fold in 1:length(segments)) {
                if (!quiet) message("Fold ", fold, " of ", length(segments))
                #Define training and test set
                test.ind <- segments[[fold]]
                training.set <- exp.df[-test.ind, ]
                training.class <- class.vec[-test.ind]
                test.set <- exp.df[test.ind, , drop = FALSE]
                test.class <- class.vec[test.ind]
                #Train J48 on training set
                if (class.algo == "J48") {
                    cv.model <- J48(training.class ~ ., training.set)
                    pred.class <- predict(cv.model, test.set)
                } else if (class.algo == "rpart") {
                    cv.model <- rpart(
                        training.class ~ ., training.set, method = "class"
                    )
                    pred.class <- predict(
                        cv.model, test.set, type = "class"
                    )
                } else{
                    stop("Unknown classification algorithm")
                }
                #Evaluate model on test set
                performance <- eval.pred(
                    pred.class, test.class, class1, performance
                )
            }
            return(performance)
        }

        cv.segments <- split(
            sample(indices), rep(1:num.folds, length = num.instances)
        )
        Rpart.performance <- c(
            "TP" = 0,
            "FN" = 0,
            "FP" = 0,
            "TN" = 0
        )
        Rpart.performance <- cross.val(exp.imput.df,
            classVector,
            cv.segments,
            Rpart.performance,
            "rpart"
        )
        if (!quiet) print(Rpart.performance)
        Rpart.confusion.matrix <- matrix(Rpart.performance, nrow = 2)
        rownames(Rpart.confusion.matrix) <- c(
            paste0("Predicted", First), paste0("Predicted", Second)
        )
        colnames(Rpart.confusion.matrix) <- c(First, Second)
        if (!quiet) print(Rpart.confusion.matrix)
        Rpart.sn <- SN(Rpart.confusion.matrix)
        Rpart.sp <- SP(Rpart.confusion.matrix)
        Rpart.acc <- ACC(Rpart.confusion.matrix)
        Rpart.mcc <- MCC(Rpart.confusion.matrix)

        if (!quiet) {
            message(
                "Rpart SN: ", Rpart.sn, "\n",
                "Rpart SP: ", Rpart.sp, "\n",
                "Rpart ACC: ", Rpart.acc, "\n",
                "Rpart MCC: ", Rpart.mcc, "\n",
            )
        }
        return(Rpart.performance)
    }