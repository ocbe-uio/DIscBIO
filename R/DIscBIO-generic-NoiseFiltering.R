#' @title Noise Filtering
#' @description Given a matrix or data frame of count data, this function
#'   estimates the size factors as follows: Each column is divided by the
#'   geometric means of the rows. The median (or, if requested, another location
#'   estimator) of these ratios (skipping the genes with a # geometric mean of
#'   zero) is used as the size factor for this column. Source: DESeq package.
#' @param object \code{DISCBIO} class object.
#' @param percentile A numeric value of the percentile. It is used to validate
#'   the ERCC spik-ins. Default is 0.8.
#' @param CV A numeric value of the coefficient of variation. It is used to
#'   validate the ERCC spik-ins. Default is 0.5.
#' @param geneCol Color of the genes that did not pass the filtration.
#' @param FgeneCol Color of the genes that passt the filtration.
#' @param erccCol Color of the ERCC spik-ins.
#' @param Val A logical vector that allows plotting only the validated ERCC
#'   spike-ins. Default is TRUE. If Val=FALSE will plot all the ERCC spike-ins.
#' @param plot A logical vector that allows plotting the technical noise.
#'   Default is TRUE.
#' @param export A logical vector that allows writing the final gene list in
#'   excel file. Default is TRUE.
#' @param quiet if `TRUE`, suppresses printed output
#' @param filename Name of the exported file (if `export=TRUE`)
#' @importFrom stats quantile var fitted.values pchisq p.adjust median aggregate
#' @importFrom graphics plot axis abline points lines
#' @importFrom statmod glmgam.fit
#' @note This function should be used only if the dataset has ERCC.
#' @return The DISCBIO-class object input with the noiseF slot filled.
#' @examples
#' sc <- DISCBIO(valuesG1msTest) # changes signature of data
#' sd_filtered <- NoiseFiltering(sc, export = FALSE)
#' str(sd_filtered)
#'
setGeneric(
  name = "NoiseFiltering",
  def = function(object, percentile = 0.8, CV = 0.3, geneCol = "yellow",
                 FgeneCol = "black", erccCol = "blue", Val = TRUE, plot = TRUE,
                 export = FALSE, quiet = FALSE,
                 filename = "Noise_filtering_genes_test") {
    standardGeneric("NoiseFiltering")
  }
)

#' @export
#' @rdname NoiseFiltering
setMethod(
  f = "NoiseFiltering",
  signature = "DISCBIO",
  definition = function(
    object, percentile, CV, geneCol, FgeneCol, erccCol, Val, plot, export,
    quiet, filename
  ) {
    if (!is.numeric(percentile)) {
      stop("percentile has to be a positive number")
    } else if (percentile <= 0) {
      stop("percentile has to be a positive number")
    }
    if (!is.numeric(CV)) {
      stop("CV has to be a positive number")
    } else if (CV <= 0) {
      stop("CV has to be a positive number")
    }

    # Split data into sub tables based on the factor data geneTypes
    GeneList <- rownames(object@expdata)
    data <- object@expdataAll
    shortNames <- substr(rownames(data), 1, 4)
    geneTypes <-
      factor(c(ENSG = "ENSG", ERCC = "ERCC")[shortNames])

    # calculate normalisation for counts\n",
    countsG1ms <- data[which(geneTypes == "ENSG"), ]
    countsERCC <- data[which(geneTypes == "ERCC"), ]

    estimateSizeFactorsForMatrix <- function(counts, locfunc = median) {
      loggeomeans <- rowMeans(log(counts))
      apply(counts, 2, function(cnts) {
        exp(locfunc((
          log(cnts) - loggeomeans
        )[is.finite(loggeomeans)]))
      })
    }

    sfERCC <- estimateSizeFactorsForMatrix(countsERCC)
    sfG1ms <- estimateSizeFactorsForMatrix(countsG1ms)

    # Divide columns by size factors to normalize counts
    nCountsERCC <- t(t(countsERCC) / sfERCC)
    nCountsG1ms <- t(t(countsG1ms) / sfG1ms)

    # perform fit, define sample moments per gene
    meansG1ms <- rowMeans(nCountsG1ms)
    varsG1ms <- apply(nCountsG1ms, 1, var)
    cv2G1ms <- varsG1ms / meansG1ms^2
    meansERCC <- rowMeans(nCountsERCC)
    varsERCC <- apply(nCountsERCC, 1, var)
    cv2ERCC <- varsERCC / meansERCC^2
    minMeanForFit <- unname(
      quantile(meansERCC[which(cv2ERCC > CV)], percentile)
    )

    if (!quiet) {
      message(
        "Cut-off value for the ERCCs= ",
        round(minMeanForFit, digits = 2),
        "\n"
      )
    }

    # Perform the fit of technical noise strength on average count.
    # We regress cv2HeLa on 1/meansForHeLa. We use the glmgam.fit function
    # from the statmod package to perform the regression as a GLM fit of
    # the gamma family with log link. The 'cbind' construct serves to
    # produce a model matrix with an intercept.
    useForFit <- meansERCC >= minMeanForFit
    fit <- glmgam.fit(
      cbind(a0 = 1, a1tilde = 1 / meansERCC[useForFit]),
      cv2ERCC[useForFit]
    )

    if (!quiet) {
      message("Coefficients of the fit:")
      print(fit$coefficients)
    }

    # To get the actual noise coefficients, we need to subtract Xi
    xi <- mean(1 / sfERCC)
    a0 <- unname(fit$coefficients["a0"])
    a1 <- unname(fit$coefficients["a1tilde"] - xi)

    # how much variance does the fit explain?
    residual <-
      var(log(fitted.values(fit)) - log(cv2ERCC[useForFit]))
    total <- var(log(cv2ERCC[useForFit]))

    if (!quiet) {
      message(
        "Explained variances of log CV^2 values= ",
        c(round(1 - residual / total, digits = 2)),
        "\n"
      )
    }

    ## Pick out genes above noise line

    # test which entries are above the line
    idx_test <- cv2G1ms > (xi + a1) / meansG1ms + a0

    # pick out genes that fulfill statement
    genes_test <- GeneList[idx_test] # pick out genes
    genes_test <- genes_test[!is.na(genes_test)] # remove na entries
    meansG1ms_test <-
      meansG1ms[idx_test] # take out mean values for fulfilled genes
    meansG1ms_test <-
      meansG1ms_test[!is.na(meansG1ms_test)] # remove na entries
    cv2G1ms_test <-
      cv2G1ms[idx_test] # take out cv2 values for fulfilled genes
    cv2G1ms_test <-
      cv2G1ms_test[!is.na(cv2G1ms_test)] # remove na entries
    genes_test <- genes_test[which(!vapply(genes_test, is.null, FALSE))]
    genes_test <- vapply(genes_test, paste0, collapse = "", "")

    if (!quiet) {
      message(
        "Number of genes that passed the filtering = ",
        length(genes_test),
        "\n"
      )
    }

    if (export) {
      write.csv(genes_test, file = paste0(filename, ".csv"))
      if (!quiet) {
        message(
          "The filtered gene list was saved as ",
          paste0(filename, ".csv")
        )
      }
    }

    if (plot) {
      plot(
        NULL,
        xaxt = "n",
        yaxt = "n",
        log = "xy",
        xlim = c(1e-1, 3e5),
        ylim = c(.005, 100),
        main = "Gene filtration by accounting for technical noise",
        xlab = "Average normalized read count",
        ylab = "Squared coefficient of variation (CV^2)"
      )
      axis(
        1,
        10^(-1:5),
        c(
          "0.1",
          "1",
          "10",
          "100",
          "1000",
          expression(10^4),
          expression(10^5)
        )
      )
      axis(2, 10^(-2:2), c("0.01", "0.1", "1", "10", "100"), las = 2)
      abline(
        h = 10^(-2:1),
        v = 10^(-1:5),
        col = "#D0D0D0",
        lwd = 2
      )
      # Plot the genes, use a different color if they are highly variable
      points(
        meansG1ms,
        cv2G1ms,
        pch = 20,
        cex = .2,
        col = geneCol
      )

      # highlight gene list from test
      points(
        meansG1ms_test,
        cv2G1ms_test,
        pch = 20,
        cex = .22,
        col = FgeneCol
      )

      # Add the technical noise fit, as before
      xg <- 10^seq(-2, 6, length.out = 1000)
      lines(xg, (xi + a1) / xg + a0, col = "red", lwd = 5)

      # Add the normalised ERCC points
      if (Val) {
        points(
          meansERCC[useForFit],
          cv2ERCC[useForFit],
          pch = 20,
          cex = 1.5,
          col = erccCol
        ) # Showing only the valied ERCCs
      } else {
        points(
          meansERCC,
          cv2ERCC,
          pch = 20,
          cex = 2,
          col = erccCol
        ) # Showing all the valied ERCCs
      }
      add_legend(
        "topleft",
        legend = c(
          "Noise curve",
          "ERCC spike-ins",
          "Genes above the noise line"
        ),
        pch = c(15, 20, 20),
        col = c("red", erccCol, FgeneCol),
        horiz = TRUE,
        bty = "n",
        cex = 0.85
      )
    }
    object@noiseF <- genes_test
    return(object)
  }
)
