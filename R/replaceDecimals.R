#' @title Replace Decimals
#' @description Replaces decimals separators between comma and periods on a
#' character vector
#' @note This function was especially designed to be used with retormatSiggenes
#' @param x vector of characters
#' @param from decimal separator on input file
#' @param to decimal separator for output file
#' @seealso reformatSiggenes
replaceDecimals <- function(x, from=",", to=".") {
    x <- gsub(",", ".", x)
    return(x)
}