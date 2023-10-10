#' @title tSNE map for K-means clustering with symbols
#' @description Visualizing the K-means clusters using tSNE maps
#' @param object \code{DISCBIO} class object.
#' @param types If types=NULL then the names of the cells will be grouped
#'   automatically. Default is NULL
#' @param legloc A keyword from the list "bottomright", "bottom", "bottomleft",
#'   "left", "topleft", "top", "topright", "right" and "center". Default is
#'   "bottomright"
setGeneric(
  "plotSymbolstSNE",
  function(object, types = NULL, legloc = "bottomright") {
    standardGeneric("plotSymbolstSNE")
  }
)

#' @export
#' @return Plot of tsne objet slot, grouped by gene.
#' @rdname plotSymbolstSNE
setMethod(
  "plotSymbolstSNE",
  signature = "DISCBIO",
  definition = function(object, types, legloc = "bottomright") {
    if (is.null(types)) {
      types <- names(object@fdata)
    }
    if (length(object@tsne) == 0) {
      stop("run comptsne before plotSymbolstSNE")
    }
    if (length(types) != ncol(object@fdata)) {
      stop(
        "types argument has wrong length.",
        "Length has to equal to the column number of object@ndata"
      )
    }

    coloc <- rainbow(length(unique(types)))
    syms <- vector()
    plot(
      object@tsne,
      xlab = "Dim 1",
      ylab = "Dim 2",
      pch = 20,
      col = "grey"
    )
    for (i in seq_len(length(unique(types)))) {
      f <- types == sort(unique(types))[i]
      syms <- append(syms, ((i - 1) %% 25) + 1)
      points(
        object@tsne[f, 1],
        object@tsne[f, 2],
        col = coloc[i],
        pch = ((i - 1) %% 25) + 1,
        cex = 1
      )
    }
    legend(legloc, legend = sort(unique(types)), col = coloc, pch = syms)
  }
)
