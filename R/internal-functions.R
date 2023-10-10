#' @importFrom fpc clusterboot cluster.stats calinhara dudahart2
#' @importFrom withr with_par with_options local_options
clustfun <- function(
    x,
    clustnr = 20,
    bootnr = 50,
    metric = "pearson",
    do.gap = TRUE,
    SE.method = "Tibs2001SEmax",
    SE.factor = .25,
    B.gap = 50,
    cln = 0,
    rseed = rseed,
    quiet = FALSE) {
  if (clustnr < 2) stop("Choose clustnr > 1")
  di <- dist.gen(t(x), method = metric)
  if (do.gap | cln > 0) {
    gpr <- NULL
    if (do.gap) {
      set.seed(rseed)
      gpr <- clusGap(
        as.matrix(di),
        FUNcluster = kmeans,
        K.max = clustnr,
        B = B.gap,
        verbose = !quiet
      )
      if (cln == 0) {
        cln <- maxSE(
          gpr$Tab[, 3],
          gpr$Tab[, 4],
          method = SE.method,
          SE.factor
        )
      }
    }
    if (cln <= 1) {
      clb <- list(
        result   = list(partition = rep(1, dim(x)[2])),
        bootmean = 1
      )
      names(clb$result$partition) <- names(x)
      return(list(x = x, clb = clb, gpr = gpr, di = di))
    }
    clb <- clusterboot(
      di,
      B             = bootnr,
      distances     = FALSE,
      bootmethod    = "boot",
      clustermethod = KmeansCBI,
      krange        = cln,
      scaling       = FALSE,
      multipleboot  = FALSE,
      bscompare     = TRUE,
      seed          = rseed,
      count         = !quiet
    )
    return(list(x = x, clb = clb, gpr = gpr, di = di))
  }
}

Kmeansruns <- function(
    data,
    krange = 2:10,
    criterion = "ch",
    iter.max = 100,
    runs = 100,
    scaledata = FALSE,
    alpha = 0.001,
    critout = FALSE,
    plot = FALSE,
    method = "euclidean",
    ...) {
  data <- as.matrix(data)
  if (criterion == "asw") sdata <- dist(data)
  if (scaledata) data <- scale(data)
  cluster1 <- 1 %in% krange
  crit <- numeric(max(krange))
  km <- list()
  for (k in krange) {
    if (k > 1) {
      minSS <- Inf
      kmopt <- NULL
      for (i in 1:runs) {
        opar <- withr::local_options(show.error.messages = FALSE)
        repeat {
          kmm <- try(kmeans(data, k))
          if (!is(kmm, "try-error")) break
        }
        opar <- withr::local_options(show.error.messages = TRUE)
        on.exit(withr::local_options(opar))
        swss <- sum(kmm$withinss)
        if (swss < minSS) {
          kmopt <- kmm
          minSS <- swss
        }
        if (plot) {
          opar <- withr::local_par(ask = TRUE)
          on.exit(withr::local_par(opar))
          pairs(data, col = kmm$cluster, main = swss)
        }
      }
      km[[k]] <- kmopt
      crit[k] <- switch(criterion,
        asw = cluster.stats(sdata, km[[k]]$cluster)$avg.silwidth,
        ch  = calinhara(data, km[[k]]$cluster)
      )
      if (critout) {
        message(k, " clusters ", crit[k], "\n")
      }
    }
  }
  if (cluster1) {
    cluster1 <- dudahart2(data, km[[2]]$cluster, alpha = alpha)$cluster1
  }
  k.best <- which.max(crit)
  if (cluster1) k.best <- 1
  km[[k.best]]$crit <- crit
  km[[k.best]]$bestk <- k.best
  out <- km[[k.best]]
  return(out)
}

KmeansCBI <- function(
    data,
    krange,
    k = NULL,
    scaling = FALSE,
    runs = 1,
    criterion = "ch",
    method = "euclidean",
    ...) {
  if (!is.null(k)) krange <- k
  if (!identical(scaling, FALSE)) {
    sdata <- scale(data, center = TRUE, scale = scaling)
  } else {
    sdata <- data
  }
  c1 <- Kmeansruns(
    sdata,
    krange,
    runs = runs,
    criterion = criterion,
    method = method,
    ...
  )
  partition <- c1$cluster
  cl <- list()
  nc <- krange
  for (i in 1:nc) cl[[i]] <- partition == i
  out <- list(
    result = c1,
    nc = nc,
    clusterlist = cl,
    partition = partition,
    clustermethod = "kmeans"
  )
  return(out)
}

dist.gen <- function(x, method = "euclidean") {
  if (method %in% c("spearman", "pearson", "kendall")) {
    as.dist(1 - cor(t(x), method = method))
  } else {
    dist(x, method = method)
  }
}

binompval <- function(p, N, n) {
  pval <- pbinom(n, round(N, 0), p, lower.tail = TRUE)
  filter <- !is.na(pval) & pval > 0.5
  pval[filter] <- 1 - pval[filter]
  return(pval)
}

add_legend <- function(...) {
  opar <- withr::local_par(
    fig = c(0, 1, 0, 1),
    oma = c(0, 0, 0, 0),
    mar = c(0, 0, 0, 0),
    new = TRUE
  )
  on.exit(withr::local_par(opar))
  plot(
    x    = 0,
    y    = 0,
    type = "n",
    bty  = "n",
    xaxt = "n",
    yaxt =  "n"
  )
  legend(...)
}

downsample <- function(x, n, dsn) {
  x <- round(x[, apply(x, 2, sum, na.rm = TRUE) >= n], 0)
  nn <- min(apply(x, 2, sum))
  for (j in 1:dsn) {
    z <- data.frame(GENEID = rownames(x))
    rownames(z) <- rownames(x)
    initv <- rep(0, nrow(z))
    for (i in seq_len(ncol(x))) {
      y <- aggregate(
        rep(1, nn), list(sample(
          rep(rownames(x), x[, i]), nn
        )),
        sum
      )
      na <- names(x)[i]
      names(y) <- c("GENEID", na)
      rownames(y) <- y$GENEID
      z[, na] <- initv
      k <- intersect(rownames(z), y$GENEID)
      z[k, na] <- y[k, na]
      z[is.na(z[, na]), na] <- 0
    }
    rownames(z) <- as.vector(z$GENEID)
    ds <- if (j == 1) z[, -1] else ds + z[, -1]
  }
  ds <- ds / dsn + .1
  return(ds)
}

eval.pred <- function(pred.class, true.class, class1, performance) {
  for (index in seq_len(length(pred.class))) {
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

#' @title Reformat Siggenes Table
#' @description Reformats the Siggenes table output from the SAMR package
#' @param table output from `samr::samr.compute.siggenes.table`
#' @seealso replaceDecimals
#' @author Waldir Leoncio
reformatSiggenes <- function(table) {
  if (is.null(table)) {
    return(table)
  }
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

#' @title Replace Decimals
#' @description Replaces decimals separators between comma and periods on a
#' character vector
#' @note This function was especially designed to be used with retormatSiggenes
#' @param x vector of characters
#' @param from decimal separator on input file
#' @param to decimal separator for output file
#' @seealso reformatSiggenes
replaceDecimals <- function(x, from = ",", to = ".") {
  x <- gsub(",", ".", x)
  return(x)
}

#' @title Prepare Example Dataset
#' @description Internal function that prepares a pre-treated dataset for use in
#' several examples
#' @param dataset Dataset used for transformation
#' @param save save results?
#' @details This function serves the purpose of treating datasets such as
#' valuesG1msReduced to reduce examples of other functions by bypassing some
#' analysis steps covered in the vignettes.
#' @return Two rda files, ones for K-means clustering and another for
#' Model-based clustering.
#' @author Waldir Leoncio
prepExampleDataset <- function(dataset, save = TRUE) {
  # ==========================================================================
  # Initial data treatment
  # ==========================================================================
  message("Treating dataset")
  sc <- DISCBIO(dataset)
  sc <- NoiseFiltering(
    sc,
    percentile = 0.9, CV = 0.2, export = FALSE, plot = FALSE, quiet = TRUE
  )
  sc <- Normalizedata(sc)
  sc <- FinalPreprocessing(sc, export = FALSE, quiet = TRUE)
  # ==========================================================================
  # Clustering
  # ==========================================================================
  message("K-means clustering")
  sc_k <- Clustexp(sc, cln = 3, quiet = TRUE)
  sc_k <- comptSNE(sc_k, quiet = TRUE)
  valuesG1msReduced_treated_K <- sc_k
  message("Model-based clustering")
  sc_mb <- Exprmclust(sc, quiet = TRUE)
  sc_mb <- comptSNE(sc_mb, rseed = 15555, quiet = TRUE)
  valuesG1msReduced_treated_MB <- sc_mb
  # ==========================================================================
  # Output
  # ==========================================================================
  message("Saving datasets")
  if (save) {
    save(
      valuesG1msReduced_treated_K,
      file = file.path("data", "valuesG1msReduced_treated_K.rda")
    )
    save(
      valuesG1msReduced_treated_MB,
      file = file.path("data", "valuesG1msReduced_treated_MB.rda")
    )
  } else {
    message("Not saving dataset because (save == FALSE)")
  }
}

#' @title Retries a URL
#' @description Retries a URL
#' @param data A gene list
#' @param species The taxonomy name/id. Default is "9606" for Homo sapiens
#' @param outputFormat format of the output. Can be "highres_image", "tsv",
#' "json", "tsv-no-header", "xml"
#' @param maxRetries maximum number of attempts to connect to the STRING api.
#' @param successCode Status code number that represents success
#' @return either the output of httr::GET or an error message
#' @importFrom httr GET status_code
#' @author Waldir Leoncio
retrieveURL <- function(
    data, species, outputFormat, maxRetries = 3, successCode = 200) {
  # ======================================================== #
  # Setting up retrieval                                     #
  # ======================================================== #
  string_api_url <- "https://string-db.org/api/"
  method <- "network"
  url <- paste0(
    string_api_url, outputFormat, "/", method, "?identifiers=",
    paste(as.character(data), collapse = "%0d"), "&species=",
    species
  )

  # ======================================================== #
  # Retrieving URL                                           #
  # ======================================================== #
  message("Retrieving URL. Please wait...")
  repos <- GET(url)
  failedGET <- status_code(repos) != successCode
  r <- 1
  while (failedGET & (r <= maxRetries)) {
    message("Failed retrieval. Retry ", r, " out of ", maxRetries, ".")
    repos <- GET(url)
    failedGET <- status_code(repos) != successCode
    r <- r + 1
  }

  # ======================================================== #
  # Final output                                             #
  # ======================================================== #
  if (failedGET) {
    stop(
      "Unable to retrieve URL. Please check the parameters ",
      "passed to the Networking() function, increase the ",
      "'maxRetries' parameter or try again later."
    )
  } else {
    message("Successful retrieval.")
    return(repos)
  }
}
