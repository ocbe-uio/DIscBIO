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
#' @examples
#' sc <- DISCBIO(valuesG1msReduced)
#' sc <- NoiseFiltering(sc, percentile=0.9, CV=0.2, export=FALSE)
#' sc <- Normalizedata(
#'     sc, mintotal=1000, minexpr=0, minnumber=0, maxexpr=Inf, downsample=FALSE,
#'     dsn=1, rseed=17000
#' )
#' sc <- FinalPreprocessing(sc, GeneFlitering="NoiseF", export=FALSE)
#' sc <- Clustexp(sc, cln=3) # K-means clustering
#' sc <- comptSNE(sc, rseed=15555)
#' cdiff <- DEGanalysis2clust(
#'     sc, Clustering="K-means", K=3, fdr=.2, name="Name", First="CL1",
#'     Second="CL2", export=FALSE
#' )
#' sigDEG <- cdiff[[1]]
#' DATAforDT <- ClassVectoringDT(
#'     sc, Clustering="K-means", K=3, First="CL1", Second="CL2", sigDEG,
#' )
#' RpartEVAL(DATAforDT,num.folds=10,First="CL1",Second="CL2")

RpartEVAL <-
    function(data,
             num.folds = 10,
             First = "CL1",
             Second = "CL2",
             quiet = FALSE) {
        exp.imput.df <- as.data.frame(t(data))
        num.instances <- nrow(exp.imput.df)
        indices <- 1:num.instances
        classVector <- factor(colnames(data))
        
        cross.val <-
            function(exp.df,
                     class.vec,
                     segments,
                     performance,
                     class.algo) {
                #Start cross validation loop
                class1 <- levels(class.vec)[1]
                for (fold in 1:length(segments)) {
                    if (!quiet)
                        cat("Fold", fold, "of", length(segments), "\n")
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
                        cv.model <- rpart(training.class ~ ., training.set, method = "class")
                        pred.class <-
                            predict(cv.model, test.set, type = "class")
                    } else{
                        stop("Unknown classification algorithm")
                    }
                    #Evaluate model on test set
                    
                    eval.pred <-
                        function(pred.class,
                                 true.class,
                                 class1,
                                 performance) {
                            for (index in 1:length(pred.class)) {
                                pred <- pred.class[index]
                                true <- true.class[index]
                                if (pred == true && true == class1) {
                                    performance["TP"] <- performance["TP"] + 1
                                } else if (pred != true && true == class1) {
                                    performance["FN"] <- performance["FN"] + 1
                                } else if (pred != true && true != class1) {
                                    performance["FP"] <- performance["FP"] + 1
                                } else if (pred == true && true != class1) {
                                    performance["TN"] <- performance["TN"] + 1
                                }
                            }
                            return(performance)
                        }
                    performance <-
                        eval.pred(pred.class, test.class, class1, performance)
                }
                return(performance)
            }
        
        cv.segments <-
            split(sample(indices), rep(1:num.folds, length = num.instances))
        Rpart.performance <- c(
            "TP" = 0,
            "FN" = 0,
            "FP" = 0,
            "TN" = 0
        )
        Rpart.performance <-
            cross.val(exp.imput.df,
                      classVector,
                      cv.segments,
                      Rpart.performance,
                      "rpart")
        if (!quiet)
            print(Rpart.performance)
        Rpart.confusion.matrix <- matrix(Rpart.performance, nrow = 2)
        rownames(Rpart.confusion.matrix) <-
            c(paste0("Predicted", First), paste0("Predicted", Second))
        colnames(Rpart.confusion.matrix) <- c(First, Second)
        if (!quiet)
            print(Rpart.confusion.matrix)
        
        SN <- function(con.mat) {
            TP <- con.mat[1, 1]
            FN <- con.mat[2, 1]
            return(TP / (TP + FN))
        }
        SP <- function(con.mat) {
            TN <- con.mat[2, 2]
            FP <- con.mat[1, 2]
            return(TN / (TN + FP))
        }
        ACC <- function(con.mat) {
            TP <- con.mat[1, 1]
            FN <- con.mat[2, 1]
            TN <- con.mat[2, 2]
            FP <- con.mat[1, 2]
            return((TP + TN) / (TP + FN + TN + FP))
        }
        MCC <- function(con.mat) {
            TP <- con.mat[1, 1]
            FN <- con.mat[2, 1]
            TN <- con.mat[2, 2]
            FP <- con.mat[1, 2]
            denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
            denom <- ifelse(denom == 0, NA, denom)
            return((TP * TN - FP * FN) / denom)
        }
        
        Rpart.sn <- SN(Rpart.confusion.matrix)
        Rpart.sp <- SP(Rpart.confusion.matrix)
        Rpart.acc <- ACC(Rpart.confusion.matrix)
        Rpart.mcc <- MCC(Rpart.confusion.matrix)
        
        if (!quiet) {
            cat(
                "Rpart SN: ",
                Rpart.sn,
                "\n",
                "Rpart SP: ",
                Rpart.sp,
                "\n",
                "Rpart ACC: ",
                Rpart.acc,
                "\n",
                "Rpart MCC: ",
                Rpart.mcc,
                "\n",
                sep = ""
            )
        }
        
        return(Rpart.performance)
    }