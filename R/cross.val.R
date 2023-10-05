cross.val <- function(
  exp.df, class.vec, segments, performance, class.algo, quiet = TRUE
) {
  # Validation
  if (!(class.algo %in% c("J48", "rpart"))) {
    stop("Unknown classification algorithm")
  }
  # Start cross validation loop
  class1 <- levels(class.vec)[1]
  for (fold in seq_len(length(segments))) {
    if (!quiet) message("Fold ", fold, " of ", length(segments))
    # Define training and test set
    test.ind <- segments[[fold]]
    training.set <- exp.df[-test.ind, ]
    test.set <- exp.df[test.ind, , drop = FALSE]
    test.set$training.class <- class.vec[-test.ind]
    test.class <- class.vec[test.ind]
    # Train J48 on training set
    if (class.algo == "J48") {
      cv.model <- J48(training.class ~ ., training.set)
      pred.class <- predict(cv.model, test.set)
    } else {
      cv.model <- rpart(training.class ~ ., training.set, method = "class")
      pred.class <- predict(cv.model, test.set, type = "class")
    }
    # Evaluate model on test set
    performance <- eval.pred(
      pred.class, test.class, class1, performance
    )
  }
  return(performance)
}
