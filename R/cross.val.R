#' @title title
#' @description description 
#' @export
#' @param exp.df exp.df
#' @param class.vec class.vec
#' @param segments segments
#' @param performance performance
#' @param class.algo class.algo
#' @importFrom stats predict
#' @importFrom rpart rpart
cross.val <- function(exp.df, class.vec, segments, performance, class.algo){
	
	#Start cross validation loop
	class1 <- levels(class.vec)[1]
	for(fold in 1:length(segments)){
		cat("Fold", fold, "of", length(segments), "\n")

		#Define training and test set
		test.ind <- segments[[fold]]
		training.set <- exp.df[-test.ind,]
		training.class <- class.vec[-test.ind]
		test.set <- exp.df[test.ind,, drop=FALSE]
		test.class <- class.vec[test.ind]
		
		#Train J48 on training set
		if(class.algo == "J48"){
			cv.model <- J48(training.class ~ ., training.set)
            pred.class <- predict(cv.model, test.set)
		} else if(class.algo == "rpart"){
			cv.model <- rpart(training.class ~ ., training.set,method="class")
		pred.class <- predict(cv.model, test.set,type="class")
		} else{
			stop("Unknown classification algorithm")
		}

		#Evaluate model on test set
#		pred.class <- predict(cv.model, test.set)
#		pred.class <- predict(cv.model, test.set,type="class")
		performance <- eval.pred(pred.class, test.class, class1, performance)
	}
	return(performance)

}