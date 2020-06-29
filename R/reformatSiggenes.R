#' @title Reformat Siggenes Table
#' @description Reformats the Siggenes table output from the SAMR package
#' @param table output from `samr::samr.compute.siggenes.table`
#' @seealso replaceDecimals
reformatSiggenes <- function(table) {
    if (is.null(table)) return(table)
    table <- as.data.frame(table)
    # ==========================================================================
    # Replacing decimal separators
    # ==========================================================================
    table[, "Score(d)"] <- replaceDecimals(table[, "Score(d)"])
    table[, "Numerator(r)"] <- replaceDecimals(table[, "Numerator(r)"])
    table[, "Denominator(s+s0)"] <- replaceDecimals(table[, "Denominator(s+s0)"])
    table[, "Fold Change"] <- replaceDecimals(table[, "Fold Change"])
    table[, "q-value(%)"] <- replaceDecimals(table[, "q-value(%)"])
    # ==========================================================================
    # Changing vector classes
    # ==========================================================================
    table[, "Row"] <- as.numeric(table[, "Row"])
    table[, "Score(d)"] <- as.numeric(table[, "Score(d)"])
    table[, "Numerator(r)"] <- as.numeric(table[, "Numerator(r)"])
    table[, "Denominator(s+s0)"] <- as.numeric(table[, "Denominator(s+s0)"])
    table[, "Fold Change"] <- as.numeric(table[, "Fold Change"])
    table[, "q-value(%)"] <- as.numeric(table[, "q-value(%)"])
    # ==========================================================================
    # Returning output
    # ==========================================================================
    return(table)
}