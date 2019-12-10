SN <- function(con.mat){
    TP <- con.mat[1,1]
    FN <- con.mat[2,1]
    return(TP/(TP+FN))
}
SP <- function(con.mat){
    TN <- con.mat[2,2]
    FP <- con.mat[1,2]
    return(TN/(TN+FP))
}
ACC <- function(con.mat){
    TP <- con.mat[1,1]
    FN <- con.mat[2,1]
    TN <- con.mat[2,2]
    FP <- con.mat[1,2]
    return((TP+TN)/(TP+FN+TN+FP))
}
MCC <- function(con.mat){
    TP <- con.mat[1,1]
    FN <- con.mat[2,1]
    TN <- con.mat[2,2]
    FP <- con.mat[1,2]
    denom <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    denom <- ifelse(denom==0, NA, denom)
    return((TP*TN-FP*FN)/denom)
}