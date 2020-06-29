#' @title Replace Decimals
#' @description Replaces decimals separators between comma and periods
#' @note This function was especially designed to be used with retormatSiggenes
#' @seealso reformatSiggenes
replaceDecimals <- function(x, from=",", to=".") {
    x <- gsub(",", ".", x)
    return(x)
}