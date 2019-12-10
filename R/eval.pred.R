#Function for counting TPs, FNs, FPs and TNs
eval.pred <- function(pred.class, true.class, class1, performance){
	for(index in 1:length(pred.class)){
		pred <- pred.class[index]
		true <- true.class[index]
		if(pred == true && true == class1){
			performance["TP"] <- performance["TP"] + 1
		} else if(pred != true && true == class1){
			performance["FN"] <- performance["FN"] + 1
		} else if(pred != true && true != class1){
			performance["FP"] <- performance["FP"] + 1
		} else if(pred == true && true != class1){
			performance["TN"] <- performance["TN"] + 1
		}
	}
	return(performance)
}