# This script contains customized versions of functions found in the samr
# package. This is necessary because samr seems to have been abandoned, so an
# upstream collaboration doesn't seem possible at the time of writing.
# ATTENTION: The source code in this file is licensed under LGPL-3.

# ==============================================================================
# Constants
# ==============================================================================
samr.const.twoclass.unpaired.response <- "Two class unpaired"
samr.const.twoclass.paired.response <- "Two class paired"
samr.const.oneclass.response <- "One class"
samr.const.quantitative.response <- "Quantitative"
samr.const.multiclass.response <- "Multiclass"
samr.const.twoclass.unpaired.timecourse.response <-
  "Two class unpaired timecourse"
samr.const.twoclass.paired.timecourse.response <- "Two class paired timecourse"
samr.const.oneclass.timecourse.response <- "One class timecourse"
samr.const.survival.response <- "Survival"
samr.const.patterndiscovery.response <- "Pattern discovery"

# ==============================================================================
# Functions
# ==============================================================================

#' @title Significance analysis of microarrays
#' @description This function is an adaptation of `samr::samr`
#' @param data Data object with components x- p by n matrix of features, one
#' observation per column (missing values allowed); y- n-vector of outcome
#' measurements; censoring.status- n-vector of censoring censoring.status
#' (1= died or event occurred, 0=survived, or event was censored), needed for a
#' censored survival outcome
#' @param resp.type Problem type: "Quantitative" for a continuous parameter
#' (Available for both array and sequencing data); "Two class unpaired" (for
#' both array and sequencing data); "Survival" for censored survival outcome
#' (for both array and sequencingdata); "Multiclass": more than 2 groups (for
#' both array and sequencing data); "One class" for a single group (only for
#' array data); "Two class paired" for two classes with paired observations
#' (for both array and sequencing data); "Two class unpaired timecourse" (only
#' for array data), "One class time course" (only for array data),
#' "Two class.paired timecourse" (only for array data), or "Pattern discovery"
#' (only for array data)
#' @param assay.type Assay type: "array" for microarray data, "seq" for counts
#' from sequencing
#' @param s0 Exchangeability factor for denominator of test statistic; Default
#' is automatic choice. Only used for array data.
#' @param s0.perc Percentile of standard deviation values to use for s0; default
#' is automatic choice; -1 means s0=0 (different from s0.perc=0, meaning
#' s0=zeroeth percentile of standard deviation values= min of sd values.
#' Only used for array data.
#' @param nperms Number of permutations used to estimate false discovery rates
#' @param center.arrays Should the data for each sample (array) be median
#' centered at the outset? Default =FALSE. Only used for array data.
#' @param testStatistic Test statistic to use in two class unpaired case.Either
#' "standard" (t-statistic) or ,"wilcoxon" (Two-sample wilcoxon or Mann-Whitney
#' test). Only used for array data.
#' @param time.summary.type Summary measure for each time course: "slope", or
#' "signed.area"). Only used for array data.
#' @param regression.method Regression method for quantitative case: "standard",
#' (linear least squares) or "ranks" (linear least squares on ranked data).
#' Only used for array data.
#' @param return.x Should the matrix of feature values be returned? Only useful
#' for time course data, where x contains summaries of the features over time.
#' Otherwise x is the same as the input data
#' @param knn.neighbors Number of nearest neighbors to use for imputation of
#' missing features values. Only used for array data.
#' @param random.seed Optional initial seed for random number generator
#' (integer)
#' @param nresamp For assay.type="seq", number of resamples used to construct
#' test statistic. Default 20. Only used for sequencing data.
#' @param nresamp.perm For assay.type="seq", number of resamples used to
#' construct test statistic for permutations. Default is equal to nresamp and it
#' must be at most nresamp. Only used for sequencing data.
#' @param xl.mode Used by Excel interface
#' @param xl.time Used by Excel interface
#' @param xl.prevfit Used by Excel interface
#' @importFrom impute impute.knn
sammy <- function(
    data, resp.type = c(
      "Quantitative", "Two class unpaired",
      "Survival", "Multiclass", "One class", "Two class paired",
      "Two class unpaired timecourse", "One class timecourse",
      "Two class paired timecourse", "Pattern discovery"
    ), assay.type = c(
      "array",
      "seq"
    ), s0 = NULL, s0.perc = NULL, nperms = 100, center.arrays = FALSE,
    testStatistic = c("standard", "wilcoxon"), time.summary.type = c(
      "slope",
      "signed.area"
    ), regression.method = c("standard", "ranks"),
    return.x = FALSE, knn.neighbors = 10, random.seed = NULL,
    nresamp = 20, nresamp.perm = NULL, xl.mode = c(
      "regular",
      "firsttime", "next20", "lasttime"
    ), xl.time = NULL, xl.prevfit = NULL) {
  this.call <- match.call()
  resp.type.arg <- match.arg(resp.type)
  assay.type <- match.arg(assay.type)
  xl.mode <- match.arg(xl.mode)
  set.seed(random.seed)
  if (is.null(nresamp.perm)) nresamp.perm <- nresamp
  nresamp.perm <- min(nresamp, nresamp.perm)
  if (xl.mode == "regular" | xl.mode == "firsttime") {
    x <- NULL
    xresamp <- NULL
    ttstar0 <- NULL
    evo <- NULL
    ystar <- NULL
    sdstar.keep <- NULL
    censoring.status <- NULL
    pi0 <- NULL
    stand.contrasts <- NULL
    stand.contrasts.star <- NULL
    stand.contrasts.95 <- NULL
    foldchange <- NULL
    foldchange.star <- NULL
    perms <- NULL
    permsy <- NULL
    eigengene <- NULL
    eigengene.number <- NULL
    testStatistic <- match.arg(testStatistic)
    time.summary.type <- match.arg(time.summary.type)
    regression.method <- match.arg(regression.method)
    x <- data$x
    y <- data$y
    argy <- y
    if (!is.null(data$eigengene.number)) {
      eigengene.number <- data$eigengene.number
    }
    if (sum(is.na(x)) > 0) {
      x <- impute.knn(x, k = knn.neighbors)
      if (!is.matrix(x)) {
        x <- x$data
      }
    }
    are.blocks.specified <- FALSE
    cond <- resp.type %in% c(
      "One class", "Two class unpaired timecourse",
      "One class unpaired timecourse", "Two class paired timecourse",
      "Pattern discovery"
    )
    if (assay.type == "seq" & cond) {
      stop(paste("Resp.type=", resp.type, " not allowed when assay.type='seq'"))
    }
    if (assay.type == "seq" & min(x) < 0) {
      stop(paste("Negative values not allowed when assay.type='seq'"))
    }
    if (assay.type == "seq" & (sum(x %% 1 != 0) != 0)) {
      stop("Non-integer values not alled when assay.type='seq'")
    }
    if (assay.type == "seq" & center.arrays) {
      stop(paste("Centering  not allowed when assay.type='seq'"))
    }
    if (assay.type == "seq" & regression.method == "ranks") {
      stop(paste("regression.method==ranks not allowed when assay.type='seq'"))
    }
    if (center.arrays) {
      x <- scale(x, center = apply(x, 2, median), scale = FALSE)
    }
    depth <- rep(NA, ncol(x))
    if (assay.type == "seq") {
      message("Estimating sequencing depths...")
      depth <- samr.estimate.depth(x)
      message("Resampling to get new data matrices...")
      xresamp <- resa(x, depth, nresamp = nresamp)
    }
    if (resp.type == samr.const.twoclass.unpaired.response) {
      if (substring(y[1], 2, 6) == "Block" | substring(
        y[1],
        2, 6
      ) == "block") {
        junk <- parse.block.labels.for.2classes(y)
        y <- junk$y
        blocky <- junk$blocky
        are.blocks.specified <- TRUE
      }
    }
    if (resp.type == samr.const.twoclass.unpaired.response |
      resp.type == samr.const.twoclass.paired.response |
      resp.type == samr.const.oneclass.response |
      resp.type == samr.const.quantitative.response |
      resp.type == samr.const.multiclass.response
    ) {
      y <- as.numeric(y)
    }
    sd.internal <- NULL
    if (resp.type == samr.const.twoclass.unpaired.timecourse.response |
      resp.type == samr.const.twoclass.paired.timecourse.response |
      resp.type == samr.const.oneclass.timecourse.response) {
      junk <- parse.time.labels.and.summarize.data(
        x, y,
        resp.type, time.summary.type
      )
      y <- junk$y
      x <- junk$x
      sd.internal <- sqrt(rowMeans(junk$sd^2))
      if (min(table(y)) == 1) {
        warning(
          "Only one timecourse in one or more classes;\n",
          "SAM plot and FDRs will be unreliable;",
          "only gene scores are informative"
        )
      }
    }
    if (resp.type == samr.const.twoclass.unpaired.timecourse.response) {
      resp.type <- samr.const.twoclass.unpaired.response
    }
    if (resp.type == samr.const.twoclass.paired.timecourse.response) {
      resp.type <- samr.const.twoclass.paired.response
    }
    if (resp.type == samr.const.oneclass.timecourse.response) {
      resp.type <- samr.const.oneclass.response
    }
    stand.contrasts <- NULL
    stand.contrasts.95 <- NULL
    if (resp.type == samr.const.survival.response) {
      censoring.status <- data$censoring.status
    }
    check.format(
      y,
      resp.type = resp.type, censoring.status = censoring.status
    )
    if (resp.type == samr.const.quantitative.response & regression.method ==
      "ranks") {
      y <- rank(y)
      x <- t(apply(x, 1, rank))
    }
    n <- nrow(x)
    sd <- NULL
    numer <- NULL
    if (resp.type == samr.const.twoclass.unpaired.response &
      testStatistic == "standard" & assay.type == "array") {
      init.fit <- ttest.func(x, y, sd = sd.internal)
      numer <- init.fit$numer
      sd <- init.fit$sd
    }
    if (resp.type == samr.const.twoclass.unpaired.response &
      testStatistic == "wilcoxon" & assay.type == "array") {
      init.fit <- wilcoxon.func(x, y)
      numer <- init.fit$numer
      sd <- init.fit$sd
    }
    if (resp.type == samr.const.oneclass.response & assay.type ==
      "array") {
      init.fit <- onesample.ttest.func(x, y, sd = sd.internal)
      numer <- init.fit$numer
      sd <- init.fit$sd
    }
    if (resp.type == samr.const.twoclass.paired.response &
      assay.type == "array") {
      init.fit <- paired.ttest.func(x, y, sd = sd.internal)
      numer <- init.fit$numer
      sd <- init.fit$sd
    }
    if (resp.type == samr.const.survival.response & assay.type ==
      "array") {
      init.fit <- cox.func(x, y, censoring.status)
      numer <- init.fit$numer
      sd <- init.fit$sd
    }
    if (resp.type == samr.const.multiclass.response & assay.type ==
      "array") {
      init.fit <- multiclass.func(x, y)
      numer <- init.fit$numer
      sd <- init.fit$sd
    }
    if (resp.type == samr.const.quantitative.response & assay.type ==
      "array") {
      init.fit <- quantitative.func(x, y)
      numer <- init.fit$numer
      sd <- init.fit$sd
    }
    if (resp.type == samr.const.patterndiscovery.response &
      assay.type == "array") {
      init.fit <- patterndiscovery.func(x)
      numer <- init.fit$numer
      sd <- init.fit$sd
    }
    if ((resp.type == samr.const.quantitative.response &
      (testStatistic == "wilcoxon" | regression.method ==
        "ranks" & assay.type == "array") | resp.type ==
      samr.const.patterndiscovery.response) | resp.type ==
      samr.const.twoclass.unpaired.response & assay.type ==
      "array" & testStatistic == "wilcoxon" | (nrow(x) <
      500) & is.null(s0) & is.null(s0.perc)) {
      s0 <- quantile(sd, 0.05)
      s0.perc <- 0.05
    }
    if (is.null(s0) & assay.type == "array") {
      if (!is.null(s0.perc)) {
        if ((s0.perc != -1 & s0.perc < 0) | s0.perc >
          100) {
          stop(
            "Illegal value for s0.perc: must be between 0 and 100, ",
            "or equal\nto (-1) (meaning that s0 should be set to zero)"
          )
        }
        if (s0.perc == -1) {
          s0 <- 0
        }
        if (s0.perc >= 0) {
          s0 <- quantile(init.fit$sd, s0.perc / 100)
        }
      }
      if (is.null(s0.perc)) {
        s0 <- est.s0(init.fit$tt, init.fit$sd)$s0.hat
        s0.perc <- 100 * sum(init.fit$sd < s0) / length(init.fit$sd)
      }
    }
    if (assay.type == "seq") {
      s0 <- 0
      s0.perc <- 0
    }
    if (resp.type == samr.const.twoclass.unpaired.response &
      testStatistic == "standard" & assay.type == "array") {
      tt <- ttest.func(x, y, s0 = s0, sd = sd.internal)$tt
    }
    if (resp.type == samr.const.twoclass.unpaired.response &
      testStatistic == "wilcoxon" & assay.type == "array") {
      tt <- wilcoxon.func(x, y, s0 = s0)$tt
    }
    if (resp.type == samr.const.oneclass.response & assay.type ==
      "array") {
      tt <- onesample.ttest.func(x, y, s0 = s0, sd = sd.internal)$tt
    }
    if (resp.type == samr.const.twoclass.paired.response &
      assay.type == "array") {
      tt <- paired.ttest.func(x, y, s0 = s0, sd = sd.internal)$tt
    }
    if (resp.type == samr.const.survival.response & assay.type ==
      "array") {
      tt <- cox.func(x, y, censoring.status, s0 = s0)$tt
    }
    if (resp.type == samr.const.multiclass.response & assay.type ==
      "array") {
      junk2 <- multiclass.func(x, y, s0 = s0)
      tt <- junk2$tt
      stand.contrasts <- junk2$stand.contrasts
    }
    if (resp.type == samr.const.quantitative.response & assay.type ==
      "array") {
      tt <- quantitative.func(x, y, s0 = s0)$tt
    }
    if (resp.type == samr.const.patterndiscovery.response &
      assay.type == "array") {
      junk <- patterndiscovery.func(
        x, s0 = s0, eigengene.number = eigengene.number
      )
      tt <- junk$tt
      eigengene <- junk$eigengene
    }
    if (resp.type == samr.const.twoclass.unpaired.response &
      assay.type == "seq") {
      junk <- wilcoxon.unpaired.seq.func(xresamp, y)
      tt <- junk$tt
      numer <- junk$numer
      sd <- junk$sd
    }
    if (resp.type == samr.const.twoclass.paired.response &
      assay.type == "seq") {
      junk <- wilcoxon.paired.seq.func(xresamp, y)
      tt <- junk$tt
      numer <- junk$numer
      sd <- junk$sd
    }
    if (resp.type == samr.const.quantitative.response & assay.type ==
      "seq") {
      junk <- quantitative.seq.func(xresamp, y)
      tt <- junk$tt
      numer <- junk$numer
      sd <- junk$sd
    }
    if (resp.type == samr.const.survival.response & assay.type ==
      "seq") {
      junk <- cox.seq.func(xresamp, y, censoring.status)
      tt <- junk$tt
      numer <- junk$numer
      sd <- junk$sd
    }
    if (resp.type == samr.const.multiclass.response & assay.type ==
      "seq") {
      junk2 <- multiclass.seq.func(xresamp, y)
      tt <- junk2$tt
      numer <- junk2$numer
      sd <- junk2$sd
      stand.contrasts <- junk2$stand.contrasts
    }
    if (
      resp.type == samr.const.quantitative.response |
        resp.type == samr.const.multiclass.response |
        resp.type == samr.const.survival.response
    ) {
      junk <- getperms(y, nperms)
      perms <- junk$perms
      all.perms.flag <- junk$all.perms.flag
      nperms.act <- junk$nperms.act
    }
    if (resp.type == samr.const.twoclass.unpaired.response) {
      if (are.blocks.specified) {
        junk <- compute.block.perms(y, blocky, nperms)
        permsy <- matrix(junk$permsy, ncol = length(y))
        all.perms.flag <- junk$all.perms.flag
        nperms.act <- junk$nperms.act
      } else {
        junk <- getperms(y, nperms)
        permsy <- matrix(y[junk$perms], ncol = length(y))
        all.perms.flag <- junk$all.perms.flag
        nperms.act <- junk$nperms.act
      }
    }
    if (resp.type == samr.const.oneclass.response) {
      if ((length(y) * log(2)) < log(nperms)) {
        allii <- 0:((2^length(y)) - 1)
        nperms.act <- 2^length(y)
        all.perms.flag <- 1
      } else {
        nperms.act <- nperms
        all.perms.flag <- 0
      }
      permsy <- matrix(NA, nrow = nperms.act, ncol = length(y))
      if (all.perms.flag == 1) {
        k <- 0
        for (i in allii) {
          junk <- integer.base.b(i, b = 2)
          if (length(junk) < length(y)) {
            junk <- c(
              rep(0, length(y) - length(junk)),
              junk
            )
          }
          k <- k + 1
          permsy[k, ] <- y * (2 * junk - 1)
        }
      } else {
        for (i in 1:nperms.act) {
          permsy[i, ] <- sample(c(-1, 1),
            size = length(y),
            replace = TRUE
          )
        }
      }
    }
    if (resp.type == samr.const.twoclass.paired.response) {
      junk <- compute.block.perms(y, abs(y), nperms)
      permsy <- junk$permsy
      all.perms.flag <- junk$all.perms.flag
      nperms.act <- junk$nperms.act
    }
    if (resp.type == samr.const.patterndiscovery.response) {
      nperms.act <- nperms
      perms <- NULL
      permsy <- NULL
      all.perms.flag <- FALSE
    }
    sdstar.keep <- NULL
    if (assay.type != "seq") {
      sdstar.keep <- matrix(0, ncol = nperms.act, nrow = nrow(x))
    }
    ttstar <- matrix(0, nrow = nrow(x), ncol = nperms.act)
    foldchange.star <- NULL
    if (resp.type == samr.const.twoclass.unpaired.response |
      resp.type == samr.const.twoclass.paired.response) {
      foldchange.star <- matrix(0, nrow = nrow(x), ncol = nperms.act)
    }
    if (resp.type == samr.const.multiclass.response) {
      stand.contrasts.star <- array(NA, c(
        nrow(x), length(table(y)),
        nperms.act
      ))
    }
  }
  if (xl.mode == "next20" | xl.mode == "lasttime") {
    evo <- xl.prevfit$evo
    tt <- xl.prevfit$tt
    numer <- xl.prevfit$numer
    eigengene <- xl.prevfit$eigengene
    eigengene.number <- xl.prevfit$eigengene.number
    sd <- xl.prevfit$sd - xl.prevfit$s0
    sd.internal <- xl.prevfit$sd.internal
    ttstar <- xl.prevfit$ttstar
    ttstar0 <- xl.prevfit$ttstar0
    n <- xl.prevfit$n
    pi0 <- xl.prevfit$pi0
    foldchange <- xl.prevfit$foldchange
    y <- xl.prevfit$y
    x <- xl.prevfit$x
    xresamp <- xl.prevfit$xresamp
    censoring.status <- xl.prevfit$censoring.status
    argy <- xl.prevfit$argy
    testStatistic <- xl.prevfit$testStatistic
    foldchange.star <- xl.prevfit$foldchange.star
    s0 <- xl.prevfit$s0
    s0.perc <- xl.prevfit$s0.perc
    resp.type <- xl.prevfit$resp.type
    resp.type.arg <- xl.prevfit$resp.type.arg
    assay.type <- xl.prevfit$assay.type
    sdstar.keep <- xl.prevfit$sdstar.keep
    resp.type <- xl.prevfit$resp.type
    stand.contrasts <- xl.prevfit$stand.contrasts
    stand.contrasts.star <- xl.prevfit$stand.contrasts.star
    stand.contrasts.95 <- xl.prevfit$stand.contrasts.95
    perms <- xl.prevfit$perms
    permsy <- xl.prevfit$permsy
    nperms <- xl.prevfit$nperms
    nperms.act <- xl.prevfit$nperms.act
    all.perms.flag <- xl.prevfit$all.perms.flag
    depth <- xl.prevfit$depth
    nresamp <- xl.prevfit$nresamp
    nresamp.perm <- xl.prevfit$nresamp.perm
  }
  if (xl.mode == "regular") {
    first <- 1
    last <- nperms.act
  }
  if (xl.mode == "firsttime") {
    first <- 1
    last <- 1
  }
  if (xl.mode == "next20") {
    first <- xl.time
    last <- min(xl.time + 19, nperms.act - 1)
  }
  if (xl.mode == "lasttime") {
    first <- nperms.act
    last <- nperms.act
  }
  for (b in first:last) {
    message(c("perm = ", b))
    if (assay.type == "array") {
      xstar <- x
    }
    if (assay.type == "seq") {
      xstar <- xresamp[, , 1:nresamp.perm]
    }
    if (resp.type == samr.const.oneclass.response) {
      ystar <- permsy[b, ]
      if (testStatistic == "standard") {
        ttstar[, b] <- onesample.ttest.func(xstar, ystar,
          s0 = s0, sd = sd.internal
        )$tt
      }
    }
    if (resp.type == samr.const.twoclass.paired.response) {
      ystar <- permsy[b, ]
      if (assay.type == "array") {
        ttstar[, b] <- paired.ttest.func(xstar, ystar,
          s0 = s0, sd = sd.internal
        )$tt
        foldchange.star[, b] <- foldchange.paired(
          xstar,
          ystar, data$logged2
        )
      }
      if (assay.type == "seq") {
        ttstar[, b] <- wilcoxon.paired.seq.func(
          xstar,
          ystar
        )$tt
        foldchange.star[, b] <- foldchange.seq.twoclass.paired(
          x,
          ystar, depth
        )
      }
    }
    if (resp.type == samr.const.twoclass.unpaired.response) {
      ystar <- permsy[b, ]
      if (assay.type == "array") {
        if (testStatistic == "standard") {
          junk <- ttest.func(xstar, ystar, s0 = s0, sd = sd.internal)
        }
        if (testStatistic == "wilcoxon") {
          junk <- wilcoxon.func(xstar, ystar, s0 = s0)
        }
        ttstar[, b] <- junk$tt
        sdstar.keep[, b] <- junk$sd
        foldchange.star[, b] <- foldchange.twoclass(
          xstar,
          ystar, data$logged2
        )
      }
      if (assay.type == "seq") {
        ttstar[, b] <- wilcoxon.unpaired.seq.func(
          xstar,
          ystar
        )$tt
        foldchange.star[, b] <- foldchange.seq.twoclass.unpaired(
          x,
          ystar, depth
        )
      }
    }
    if (resp.type == samr.const.survival.response) {
      o <- perms[b, ]
      if (assay.type == "array") {
        ttstar[, b] <- cox.func(
          xstar, y[o],
          censoring.status = censoring.status[o], s0 = s0
        )$tt
      }
      if (assay.type == "seq") {
        ttstar[, b] <- cox.seq.func(
          xstar, y[o],
          censoring.status = censoring.status[o]
        )$tt
      }
    }
    if (resp.type == samr.const.multiclass.response) {
      ystar <- y[perms[b, ]]
      if (assay.type == "array") {
        junk <- multiclass.func(xstar, ystar, s0 = s0)
        ttstar[, b] <- junk$tt
        sdstar.keep[, b] <- junk$sd
        stand.contrasts.star[, , b] <- junk$stand.contrasts
      }
      if (assay.type == "seq") {
        junk <- multiclass.seq.func(xstar, ystar)
        ttstar[, b] <- junk$tt
        stand.contrasts.star[, , b] <- junk$stand.contrasts
      }
    }
    if (resp.type == samr.const.quantitative.response) {
      ystar <- y[perms[b, ]]
      if (assay.type == "array") {
        junk <- quantitative.func(xstar, ystar, s0 = s0)
        ttstar[, b] <- junk$tt
        sdstar.keep[, b] <- junk$sd
      }
      if (assay.type == "seq") {
        junk <- quantitative.seq.func(xstar, ystar)
        ttstar[, b] <- junk$tt
      }
    }
    if (resp.type == samr.const.patterndiscovery.response) {
      xstar <- permute.rows(x)
      junk <- patterndiscovery.func(
        xstar,
        s0 = s0, eigengene.number = eigengene.number
      )
      ttstar[, b] <- junk$tt
      sdstar.keep[, b] <- junk$sd
    }
  }
  if (xl.mode == "regular" | xl.mode == "lasttime") {
    ttstar0 <- ttstar
    for (j in seq_len(ncol(ttstar))) {
      ttstar[, j] <- -1 * sort(-1 * ttstar[, j])
    }
    for (i in seq_len(nrow(ttstar))) {
      ttstar[i, ] <- sort(ttstar[i, ])
    }
    evo <- apply(ttstar, 1, mean)
    evo <- evo[seq(length(evo), 1)]
    pi0 <- 1
    if (resp.type != samr.const.multiclass.response) {
      qq <- quantile(ttstar, c(0.25, 0.75))
    }
    if (resp.type == samr.const.multiclass.response) {
      qq <- quantile(ttstar, c(0, 0.5))
    }
    pi0 <- sum(tt > qq[1] & tt < qq[2]) / (0.5 * length(tt))
    foldchange <- NULL
    if (resp.type == samr.const.twoclass.unpaired.response &
      assay.type == "array") {
      foldchange <- foldchange.twoclass(x, y, data$logged2)
    }
    if (resp.type == samr.const.twoclass.paired.response &
      assay.type == "array") {
      foldchange <- foldchange.paired(x, y, data$logged2)
    }
    if (resp.type == samr.const.oneclass.response & assay.type ==
      "array") {
    }
    stand.contrasts.95 <- NULL
    if (resp.type == samr.const.multiclass.response) {
      stand.contrasts.95 <- quantile(
        stand.contrasts.star,
        c(0.025, 0.975)
      )
    }
    if (resp.type == samr.const.twoclass.unpaired.response &
      assay.type == "seq") {
      foldchange <- foldchange.seq.twoclass.unpaired(
        x,
        y, depth
      )
    }
    if (resp.type == samr.const.twoclass.paired.response &
      assay.type == "seq") {
      foldchange <- foldchange.seq.twoclass.paired(
        x, y,
        depth
      )
    }
    if (return.x == FALSE) {
      x <- NULL
    }
  }
  return(
    list(
      n = n, x = x, xresamp = xresamp, y = y, argy = argy,
      censoring.status = censoring.status, testStatistic = testStatistic,
      nperms = nperms, nperms.act = nperms.act, tt = tt, numer = numer,
      sd = sd + s0, sd.internal = sd.internal, s0 = s0, s0.perc = s0.perc,
      evo = evo, perms = perms, permsy = permsy, nresamp = nresamp,
      nresamp.perm = nresamp.perm, all.perms.flag = all.perms.flag,
      ttstar = ttstar, ttstar0 = ttstar0, eigengene = eigengene,
      eigengene.number = eigengene.number, pi0 = pi0, foldchange = foldchange,
      foldchange.star = foldchange.star, sdstar.keep = sdstar.keep,
      resp.type = resp.type, resp.type.arg = resp.type.arg,
      assay.type = assay.type, stand.contrasts = stand.contrasts,
      stand.contrasts.star = stand.contrasts.star,
      stand.contrasts.95 = stand.contrasts.95,
      depth = depth, call = this.call
    )
  )
}

#' @title Estimate sequencing depths
#' @param x data matrix. nrow=#gene, ncol=#sample
#' @return depth: estimated sequencing depth. a vector with len sample.
samr.estimate.depth <- function(x) {
  iter <- 5
  cmeans <- colSums(x) / sum(x)
  for (i in 1:iter) {
    n0 <- rowSums(x) %*% t(cmeans)
    prop <- rowSums((x - n0)^2 / (n0 + 1e-08))
    qs <- quantile(prop, c(0.25, 0.75))
    keep <- (prop >= qs[1]) & (prop <= qs[2])
    cmeans <- colMeans(x[keep, ])
    cmeans <- cmeans / sum(cmeans)
  }
  depth <- cmeans / mean(cmeans)
  return(depth)
}

#' @title Resampling
#' @param x data matrix. nrow=#gene, ncol=#sample
#' @param d estimated sequencing depth
#' @param nresamp number of resamplings
#' @return xresamp: an rank array with dim #gene*#sample*nresamp
#' @description Corresponds to `samr::resample`
#' @importFrom stats rpois runif
resa <- function(x, d, nresamp = 20) {
  ng <- nrow(x)
  ns <- ncol(x)
  dbar <- exp(mean(log(d)))
  xresamp <- array(0, dim = c(ng, ns, nresamp))
  for (k in 1:nresamp) {
    for (j in 1:ns) {
      xresamp[, j, k] <- rpois(n = ng, lambda = (dbar / d[j]) *
        x[, j]) + runif(ng) * 0.1
    }
  }
  for (k in 1:nresamp) {
    xresamp[, , k] <- t(rankcols(t(xresamp[, , k])))
  }
  return(xresamp)
}

#' @title Rank columns
#' @description Ranks the elements within each col of the matrix x and returns
#' these ranks in a matrix
#' @note this function is equivalent to `samr::rankcol`, but uses `apply` to
#' rank the colums instead of a compiled Fortran function which was causing our
#' DEGanalysis functions to freeze in large datasets.
#' @param x x
rankcols <- function(x) {
  # ranks the elements within each col of the matrix x
  # and returns these ranks in a matrix
  n <- nrow(x)
  p <- ncol(x)
  mode(n) <- "integer"
  mode(p) <- "integer"
  mode(x) <- "double"
  xr <- apply(x, 2, rank)
  return(xr)
}

#' @title Check format
#' @param y y
#' @param resp.type resp type
#' @param censoring.status censoring status
check.format <- function(y, resp.type, censoring.status = NULL) {
  # here i do some format checks for the input data$y
  # note that checks for time course data are done in the
  #   parse function for time course;
  #  we then check the output from the parser in this function
  if (resp.type == samr.const.twoclass.unpaired.response |
    resp.type == samr.const.twoclass.unpaired.timecourse.response) {
    if (sum(y == 1) + sum(y == 2) != length(y)) {
      stop(paste(
        "Error in input response data: response type ",
        resp.type, " specified; values must be 1 or 2"
      ))
    }
  }
  if (resp.type == samr.const.twoclass.paired.response | resp.type ==
    samr.const.twoclass.paired.timecourse.response) {
    if (sum(y) != 0) {
      stop(paste(
        "Error in input response data: response type ",
        resp.type, " specified; values must be -1, 1, -2, 2, etc"
      ))
    }
    if (sum(table(y[y > 0]) != abs(table(y[y < 0])))) {
      stop(paste(
        "Error in input response data:  response type ",
        resp.type, " specified; values must be -1, 1, -2, 2, etc"
      ))
    }
  }
  if (resp.type == samr.const.oneclass.response | resp.type ==
    samr.const.oneclass.timecourse.response) {
    if (sum(y == 1) != length(y)) {
      stop(paste(
        "Error in input response data: response type ",
        resp.type, " specified;  values must all be 1"
      ))
    }
  }
  if (resp.type == samr.const.multiclass.response) {
    tt <- table(y)
    nc <- length(tt)
    if (sum(y <= nc & y > 0) < length(y)) {
      stop(
        "Error in input response data: response type ", resp.type,
        " specified; values must be 1,2, ... number of classes"
      )
    }
    for (k in 1:nc) {
      if (sum(y == k) < 2) {
        stop(paste(
          "Error in input response data: response type ",
          resp.type, " specified; there must be >1 sample per class"
        ))
      }
    }
  }
  if (resp.type == samr.const.quantitative.response) {
    if (!is.numeric(y)) {
      stop(paste(
        "Error in input response data: response type",
        resp.type, " specified; values must be numeric"
      ))
    }
  }
  if (resp.type == samr.const.survival.response) {
    if (is.null(censoring.status)) {
      stop(paste(
        "Error in input response data: response type ",
        resp.type, " specified; error in censoring indicator"
      ))
    }
    if (!is.numeric(y) | sum(y < 0) > 0) {
      stop(
        "Error in input response data:  response type ", resp.type,
        " specified; survival times  must be numeric and nonnegative"
      )
      if (sum(censoring.status == 0) + sum(censoring.status ==
        1) != length(censoring.status)) {
        stop(
          "Error in input response data: response type ", resp.type,
          " specified; censoring indicators must be ",
          "0 (censored) or 1 (failed)"
        )
      }
    }
    if (sum(censoring.status == 1) < 1) {
      stop(paste(
        "Error in input response data:   response type ",
        resp.type, " specified; there are no uncensored observations"
      ))
    }
  }
  return()
}

#' @title Twoclass Wilcoxon statistics
#' @param xresamp an rank array with dim #gene*#sample*nresamp
#' @param y outcome vector of values 1 and 2
#' @return the statistic.
wilcoxon.unpaired.seq.func <- function(xresamp, y) {
  tt <- rep(0, dim(xresamp)[1])
  for (i in seq_len(dim(xresamp)[3])) {
    tt <- tt + rowSums(xresamp[, y == 2, i]) - sum(y == 2) *
      (length(y) + 1) / 2
  }
  tt <- tt / dim(xresamp)[3]
  return(list(tt = tt, numer = tt, sd = rep(1, length(tt))))
}
wilcoxon.paired.seq.func <- function(xresamp, y) {
  tt <- rep(0, dim(xresamp)[1])
  for (i in seq_len(dim(xresamp)[3])) {
    tt <- tt + rowSums(xresamp[, y > 0, i]) - sum(y > 0) *
      (length(y) + 1) / 2
  }
  tt <- tt / dim(xresamp)[3]
  return(list(tt = tt, numer = tt, sd = rep(1, length(tt))))
}
getperms <- function(y, nperms) {
  total.perms <- factorial(length(y))
  if (total.perms <= nperms) {
    perms <- permute(seq_len(length(y)))
    all.perms.flag <- 1
    nperms.act <- total.perms
  }
  if (total.perms > nperms) {
    perms <- matrix(NA, nrow = nperms, ncol = length(y))
    for (i in 1:nperms) {
      perms[i, ] <- sample(seq_len(length(y)), size = length(y))
    }
    all.perms.flag <- 0
    nperms.act <- nperms
  }
  return(list(
    perms = perms, all.perms.flag = all.perms.flag,
    nperms.act = nperms.act
  ))
}
foldchange.twoclass <- function(x, y, logged2) {
  m1 <- rowMeans(x[, y == 1, drop = FALSE])
  m2 <- rowMeans(x[, y == 2, drop = FALSE])
  if (!logged2) {
    fc <- m2 / m1
  }
  if (logged2) {
    fc <- 2^{
      m2 - m1
    }
  }
  return(fc)
}
#' @title Foldchange of twoclass unpaired sequencing data
#' @param x x
#' @param y y
#' @param depth depth
foldchange.seq.twoclass.unpaired <- function(x, y, depth) {
  x.norm <- scale(x, center = FALSE, scale = depth) + 1e-08
  fc <- apply(x.norm[, y == 2], 1, median) /
    apply(x.norm[, y ==
      1], 1, median)
  return(fc)
}
foldchange.seq.twoclass.paired <- function(x, y, depth) {
  nc <- ncol(x) / 2
  o1 <- o2 <- rep(0, nc)
  for (j in 1:nc) {
    o1[j] <- which(y == -j)
    o2[j] <- which(y == j)
  }
  x.norm <- scale(x, center = FALSE, scale = depth) + 1e-08
  d <- x.norm[, o2, drop = FALSE] / x.norm[, o1, drop = FALSE]
  fc <- lapply(d, 1, function(x) median(x, na.rm = TRUE))
  return(fc)
}
permute <- function(elem) {
  # generates all perms of the vector elem
  if (!missing(elem)) {
    if (length(elem) == 2) {
      return(matrix(c(elem, elem[2], elem[1]), nrow = 2))
    }
    last.matrix <- permute(elem[-1])
    dim.last <- dim(last.matrix)
    new.matrix <- matrix(0, nrow = dim.last[1] * (dim.last[2] +
      1), ncol = dim.last[2] + 1)
    for (row in 1:(dim.last[1])) {
      for (col in 1:(dim.last[2] + 1)) {
        new.matrix[row + (col - 1) * dim.last[1], ] <-
          insert.value(last.matrix[row, ], elem[1], col)
      }
    }
    return(new.matrix)
  } else {
    message("Usage: permute(elem)\n\twhere elem is a vector")
  }
}
insert.value <- function(vec, newval, pos) {
  if (pos == 1) {
    return(c(newval, vec))
  }
  lvec <- length(vec)
  if (pos > lvec) {
    return(c(vec, newval))
  }
  return(c(vec[1:pos - 1], newval, vec[pos:lvec]))
}
parse.block.labels.for.2classes <- function(y) {
  # this only works for 2 class case- having form jBlockn,
  #   where j=1 or 2
  n <- length(y)
  y.act <- rep(NA, n)
  blocky <- rep(NA, n)
  for (i in 1:n) {
    blocky[i] <- as.numeric(substring(y[i], 7, nchar(y[i])))
    y.act[i] <- as.numeric(substring(y[i], 1, 1))
  }
  return(list(y.act = y.act, blocky = blocky))
}
parse.time.labels.and.summarize.data <- function(
    x,
    y, resp.type, time.summary.type) {
  # parse time labels, and summarize time data for each
  #   person, via a slope or area
  # does some error checking too
  n <- length(y)
  last5char <- rep(NA, n)
  last3char <- rep(NA, n)
  for (i in 1:n) {
    last3char[i] <- substring(y[i], nchar(y[i]) - 2, nchar(y[i]))
    last5char[i] <- substring(y[i], nchar(y[i]) - 4, nchar(y[i]))
  }
  if (sum(last3char == "End") != sum(last5char == "Start")) {
    stop("Error in format of time course data: a Start or End tag is missing")
  }
  y.act <- rep(NA, n)
  timey <- rep(NA, n)
  person.id <- rep(NA, n)
  k <- 1
  end.flag <- FALSE
  person.id[1] <- 1
  if (substring(y[1], nchar(y[1]) - 4, nchar(y[1])) != "Start") {
    stop(
      "Error in format of time course data: ",
      "first cell should have a Start tag"
    )
  }
  for (i in 1:n) {
    message(i)
    j <- 1
    while (substring(y[i], j, j) != "T") {
      j <- j + 1
    }
    end.of.y <- j - 1
    y.act[i] <- as.numeric(substring(y[i], 1, end.of.y))
    timey[i] <- substring(y[i], end.of.y + 5, nchar(y[i]))
    if (nchar(timey[i]) > 3 & substring(timey[i], nchar(timey[i]) -
      2, nchar(timey[i])) == "End") {
      end.flag <- TRUE
      timey[i] <- substring(timey[i], 1, nchar(timey[i]) -
        3)
    }
    if (nchar(timey[i]) > 3 & substring(timey[i], nchar(timey[i]) -
      4, nchar(timey[i])) == "Start") {
      timey[i] <- substring(timey[i], 1, nchar(timey[i]) -
        5)
    }
    if (i < n & !end.flag) {
      person.id[i + 1] <- k
    }
    if (i < n & end.flag) {
      k <- k + 1
      person.id[i + 1] <- k
    }
    end.flag <- FALSE
  }
  timey <- as.numeric(timey)
  # do a check that the format was correct
  tt <- table(person.id, y.act)
  junk <- function(x) {
    sum(x != 0)
  }
  if (sum(apply(tt, 1, junk) != 1) > 0) {
    num <- seq_len(nrow(tt))[apply(tt, 1, junk) > 1]
    stop(paste(
      "Error in format of  time course data, timecourse #",
      as.character(num)
    ))
  }
  npeople <- length(unique(person.id))
  newx <- matrix(NA, nrow = nrow(x), ncol = npeople)
  sd <- matrix(NA, nrow = nrow(x), ncol = npeople)
  for (j in 1:npeople) {
    jj <- person.id == j
    tim <- timey[jj]
    xc <- t(scale(t(x[, jj, drop = FALSE]), center = TRUE, scale = FALSE))
    if (time.summary.type == "slope") {
      junk <- quantitative.func(xc, tim - mean(tim))
      newx[, j] <- junk$numer
      sd[, j] <- junk$sd
    }
    if (time.summary.type == "signed.area") {
      junk <- timearea.func(x[, jj, drop = FALSE], tim)
      newx[, j] <- junk$numer
      sd[, j] <- junk$sd
    }
  }
  y.unique <- y.act[!duplicated(person.id)]
  return(list(y = y.unique, x = newx, sd = sd))
}
ttest.func <- function(x, y, s0 = 0, sd = NULL) {
  n1 <- sum(y == 1)
  n2 <- sum(y == 2)
  m1 <- rowMeans(x[, y == 1, drop = FALSE])
  m2 <- rowMeans(x[, y == 2, drop = FALSE])
  if (is.null(sd)) {
    sd <- sqrt(((n2 - 1) * varr(x[, y == 2], meanx = m2) +
      (n1 - 1) * varr(x[, y == 1], meanx = m1)) * (1 / n1 +
      1 / n2) / (n1 + n2 - 2))
  }
  numer <- m2 - m1
  dif.obs <- (numer) / (sd + s0)
  return(list(tt = dif.obs, numer = numer, sd = sd))
}

wilcoxon.func <- function(x, y, s0 = 0) {
  n1 <- sum(y == 1)
  n2 <- sum(y == 2)
  p <- nrow(x)
  r2 <- rowSums(t(apply(x, 1, rank))[, y == 2, drop = FALSE])
  numer <- r2 - (n2 / 2) * (n2 + 1) - (n1 * n2) / 2
  sd <- sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
  tt <- (numer) / (sd + s0)
  return(list(tt = tt, numer = numer, sd = rep(sd, p)))
}

onesample.ttest.func <- function(x, y, s0 = 0, sd = NULL) {
  n <- length(y)
  x <- x * matrix(y, nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
  m <- rowMeans(x)
  if (is.null(sd)) {
    sd <- sqrt(varr(x, meanx = m) / n)
  }
  dif.obs <- m / (sd + s0)
  return(list(tt = dif.obs, numer = m, sd = sd))
}

patterndiscovery.func <- function(x, s0 = 0, eigengene.number = 1) {
  a <- mysvd(x, n.components = eigengene.number)
  v <- a$v[, eigengene.number]
  # here we try to guess the most interpretable orientation
  #   for the eigengene
  om <- abs(a$u[, eigengene.number]) > quantile(
    abs(a$u[, eigengene.number]),
    0.95
  )
  if (median(a$u[om, eigengene.number]) < 0) {
    v <- -1 * v
  }
  aa <- quantitative.func(x, v, s0 = s0)
  eigengene <- cbind(seq_len(nrow(a$v)), v)
  dimnames(eigengene) <- list(NULL, c("sample number", "value"))
  return(
    list(tt = aa$tt, numer = aa$numer, sd = aa$sd, eigengene = eigengene)
  )
}

paired.ttest.func <- function(x, y, s0 = 0, sd = NULL) {
  nc <- ncol(x) / 2
  o <- 1:nc
  o1 <- rep(0, ncol(x) / 2)
  o2 <- o1
  for (j in 1:nc) {
    o1[j] <- seq_len(ncol(x))[y == -o[j]]
  }
  for (j in 1:nc) {
    o2[j] <- seq_len(ncol(x))[y == o[j]]
  }
  d <- x[, o2, drop = FALSE] - x[, o1, drop = FALSE]
  if (is.matrix(d)) {
    m <- rowMeans(d)
  }
  if (!is.matrix(d)) {
    m <- mean(d)
  }
  if (is.null(sd)) {
    if (is.matrix(d)) {
      sd <- sqrt(varr(d, meanx = m) / nc)
    }
    if (!is.matrix(d)) {
      sd <- sqrt(var(d) / nc)
    }
  }
  dif.obs <- m / (sd + s0)
  return(list(tt = dif.obs, numer = m, sd = sd))
}

cox.func <- function(x, y, censoring.status, s0 = 0) {
  # find the index matrix
  Dn <- sum(censoring.status == 1)
  Dset <- seq_len(ncol(x))[censoring.status == 1] # the set of observed
  ind <- matrix(0, ncol(x), Dn)
  # get the matrix
  for (i in 1:Dn) {
    ind[y > y[Dset[i]] - 1e-08, i] <- 1 / sum(y > y[Dset[i]] -
      1e-08)
  }
  ind.sums <- rowSums(ind)
  x.ind <- x %*% ind
  # get the derivatives
  numer <- x %*% (censoring.status - ind.sums)
  sd <- sqrt((x * x) %*% ind.sums - rowSums(x.ind * x.ind))
  tt <- numer / (sd + s0)
  return(list(tt = tt, numer = numer, sd = sd))
}

multiclass.func <- function(x, y, s0 = 0) {
  ## assumes y is coded 1,2...
  nn <- table(y)
  m <- matrix(0, nrow = nrow(x), ncol = length(nn))
  v <- m
  for (j in seq_len(length(nn))) {
    m[, j] <- rowMeans(x[, y == j])
    v[, j] <- (nn[j] - 1) * varr(x[, y == j], meanx = m[
      ,
      j
    ])
  }
  mbar <- rowMeans(x)
  mm <- m - matrix(mbar, nrow = length(mbar), ncol = length(nn))
  fac <- (sum(nn) / prod(nn))
  scor <- sqrt(fac * (apply(matrix(nn,
    nrow = nrow(m), ncol = ncol(m),
    byrow = TRUE
  ) * mm * mm, 1, sum)))
  sd <- sqrt(rowSums(v) * (1 / sum(nn - 1)) * sum(1 / nn))
  tt <- scor / (sd + s0)
  mm.stand <- t(scale(t(mm), center = FALSE, scale = sd))
  return(list(tt = tt, numer = scor, sd = sd, stand.contrasts = mm.stand))
}

est.s0 <- function(tt, sd, s0.perc = seq(0, 1, by = 0.05)) {
  ## estimate s0 (exchangeability) factor for denominator.
  ## returns the actual estimate s0 (not a percentile)
  br <- unique(quantile(sd, seq(0, 1, len = 101)))
  nbr <- length(br)
  a <- cut(sd, br, labels = FALSE)
  a[is.na(a)] <- 1
  cv.sd <- rep(0, length(s0.perc))
  for (j in seq_len(length(s0.perc))) {
    w <- quantile(sd, s0.perc[j])
    w[j == 1] <- 0
    tt2 <- tt * sd / (sd + w)
    tt2[tt2 == Inf] <- NA
    sds <- rep(0, nbr - 1)
    for (i in 1:(nbr - 1)) {
      sds[i] <- stats::mad(tt2[a == i], na.rm = TRUE)
    }
    cv.sd[j] <- sqrt(var(sds)) / mean(sds)
  }
  o <- seq_len(length(s0.perc))[cv.sd == min(cv.sd)]
  # we don;t allow taking s0.hat to be 0th percentile when
  #   min sd is 0
  s0.hat <- quantile(sd[sd != 0], s0.perc[o])
  return(list(s0.perc = s0.perc, cv.sd = cv.sd, s0.hat = s0.hat))
}

permute.rows <- function(x) {
  dd <- dim(x)
  n <- dd[1]
  p <- dd[2]
  mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
  matrix(t(x)[order(mm)], n, p, byrow = TRUE)
}

foldchange.paired <- function(x, y, logged2) {
  nc <- ncol(x) / 2
  o <- 1:nc
  o1 <- rep(0, ncol(x) / 2)
  o2 <- o1
  for (j in 1:nc) {
    o1[j] <- seq_len(ncol(x))[y == -o[j]]
  }
  for (j in 1:nc) {
    o2[j] <- seq_len(ncol(x))[y == o[j]]
  }
  if (!logged2) {
    d <- x[, o2, drop = FALSE] / x[, o1, drop = FALSE]
  }
  if (logged2) {
    d <- x[, o2, drop = FALSE] - x[, o1, drop = FALSE]
  }
  if (!logged2) {
    fc <- rowMeans(d)
  }
  if (logged2) {
    fc <- 2^rowMeans(d)
  }
  return(fc)
}
foldchange.seq.twoclass.unpaired <- function(x, y, depth) {
  x.norm <- scale(x, center = FALSE, scale = depth) + 1e-08
  fc <- apply(x.norm[, y == 2], 1, median) /
    apply(x.norm[, y == 1], 1, median)
  return(fc)
}
integer.base.b <- function(x, b = 2) {
  xi <- as.integer(x)
  if (xi == 0) {
    return(0)
  }
  if (any(is.na(xi) | ((x - xi) != 0))) {
    print(list(ERROR = "x not integer", x = x))
  }
  N <- length(x)
  xMax <- max(x)
  ndigits <- (floor(logb(xMax, base = 2)) + 1)
  Base.b <- array(NA, dim = c(N, ndigits))
  for (i in 1:ndigits) {
    Base.b[, ndigits - i + 1] <- (x %% b)
    x <- (x %/% b)
  }
  if (N == 1) {
    Base.b[1, ]
  } else {
    Base.b
  }
}
varr <- function(x, meanx = NULL) {
  n <- ncol(x)
  p <- nrow(x)
  Y <- matrix(1, nrow = n, ncol = 1)
  if (is.null(meanx)) {
    meanx <- rowMeans(x)
  }
  ans <- rep(1, p)
  xdif <- x - meanx %*% t(Y)
  ans <- (xdif^2) %*% rep(1 / (n - 1), n)
  ans <- drop(ans)
  return(ans)
}
quantitative.func <- function(x, y, s0 = 0) {
  # regression of x on y
  my <- mean(y)
  yy <- y - my
  temp <- x %*% yy
  mx <- rowMeans(x)
  syy <- sum(yy^2)
  scor <- temp / syy
  b0hat <- mx - scor * my
  ym <- matrix(y, nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
  xhat <- matrix(b0hat, nrow = nrow(x), ncol = ncol(x)) + ym *
    matrix(scor, nrow = nrow(x), ncol = ncol(x))
  sigma <- sqrt(rowSums((x - xhat)^2) / (ncol(xhat) - 2))
  sd <- sigma / sqrt(syy)
  tt <- scor / (sd + s0)
  return(list(tt = tt, numer = scor, sd = sd))
}
timearea.func <- function(x, y, s0 = 0) {
  n <- ncol(x)
  xx <- 0.5 * (x[, 2:n] + x[, 1:(n - 1)]) * matrix(diff(y),
    nrow = nrow(x), ncol = n - 1, byrow = TRUE
  )
  numer <- rowMeans(xx)
  sd <- sqrt(varr(xx, meanx = numer) / n)
  tt <- numer / sqrt(sd + s0)
  return(list(tt = tt, numer = numer, sd = sd))
}
cox.seq.func <- function(xresamp, y, censoring.status) {
  # get the dimensions
  ns <- ncol(xresamp)
  # prepare for the calculation
  # find the index matrix
  Dn <- sum(censoring.status == 1)
  Dset <- seq_len(ns)[censoring.status == 1] # the set of died
  ind <- matrix(0, ns, Dn)
  # get the matrix
  for (i in 1:Dn) {
    ind[y >= y[Dset[i]] - 1e-08, i] <- 1 / sum(y >= y[Dset[i]] -
      1e-08)
  }
  ind.sums <- rowSums(ind)
  # calculate the score statistic
  tt <- apply(xresamp, 3, function(x, cen.ind, ind.para, ind.sums.para) {
    dev1 <- x %*% cen.ind
    x.ind <- x %*% ind.para
    dev2 <- (x * x) %*% ind.sums.para - rowSums(x.ind * x.ind)
    dev1 / (sqrt(dev2) + 1e-08)
  }, (censoring.status - ind.sums), ind, ind.sums)
  tt <- rowMeans(tt)
  return(list(tt = tt, numer = tt, sd = rep(1, length(tt))))
}
compute.block.perms <- function(y, blocky, nperms) {
  # y are the data (eg class label 1 vs 2; or -1,1, -2,2 for
  #   paired data)
  # blocky are the block labels (abs(y) for paired daatr)
  ny <- length(y)
  nblocks <- length(unique(blocky))
  tab <- table(blocky)
  total.nperms <- prod(factorial(tab))
  # block.perms is a list of all possible permutations
  block.perms <- vector("list", nblocks)
  # first enumerate all perms, when possible
  if (total.nperms <= nperms) {
    all.perms.flag <- 1
    nperms.act <- total.nperms
    for (i in 1:nblocks) {
      block.perms[[i]] <- permute(y[blocky == i])
    }
    kk <- 0:(factorial(max(tab))^nblocks - 1)
    # the rows of the matrix outerm runs through the 'outer
    #   product'
    # first we assume that all blocks have max(tab) members;
    #   then we remove rows of outerm that
    #  are illegal (ie when a block has fewer members)
    outerm <- matrix(0, nrow = length(kk), ncol = nblocks)
    for (i in seq_len(length(kk))) {
      kkkk <- integer.base.b(kk[i], b = factorial(max(tab)))
      if (length(kkkk) > nblocks) {
        kkkk <- kkkk[(length(kkkk) - nblocks + 1):length(kkkk)]
      }
      outerm[i, (nblocks - length(kkkk) + 1):nblocks] <- kkkk
    }
    outerm <- outerm + 1
    # now remove rows that are illegal perms
    ind <- rep(TRUE, nrow(outerm))
    for (j in seq_len(ncol(outerm))) {
      ind <- ind & outerm[, j] <= factorial(tab[j])
    }
    outerm <- outerm[ind, , drop = FALSE]
    # finally, construct permutation matrix from outer product
    permsy <- matrix(NA, nrow = total.nperms, ncol = ny)
    for (i in 1:total.nperms) {
      junk <- NULL
      for (j in 1:nblocks) {
        junk <- c(junk, block.perms[[j]][outerm[i, j], ])
      }
      permsy[i, ] <- junk
    }
  }
  # next handle case when there are too many perms to enumerate
  if (total.nperms > nperms) {
    all.perms.flag <- 0
    nperms.act <- nperms
    permsy <- NULL
    block.perms <- vector("list", nblocks)
    for (j in 1:nblocks) {
      block.perms[[j]] <- sample.perms(y[blocky == j], nperms = nperms)
    }
    for (j in 1:nblocks) {
      permsy <- cbind(permsy, block.perms[[j]])
    }
  }
  return(list(
    permsy = permsy, all.perms.flag = all.perms.flag,
    nperms.act = nperms.act
  ))
}
sample.perms <- function(elem, nperms) {
  # randomly generates  nperms of the vector elem
  res <- permute.rows(matrix(elem,
    nrow = nperms, ncol = length(elem),
    byrow = TRUE
  ))
  return(res)
}
mysvd <- function(x, n.components = NULL) {
  # finds PCs of matrix x
  p <- nrow(x)
  n <- ncol(x)
  # center the observations (rows)
  feature.means <- rowMeans(x)
  x <- t(scale(t(x), center = feature.means, scale = FALSE))
  if (is.null(n.components)) {
    n.components <- min(n, p)
  }
  if (p > n) {
    a <- eigen(t(x) %*% x)
    v <- a$vec[, 1:n.components, drop = FALSE]
    d <- sqrt(a$val[1:n.components, drop = FALSE])
    u <- scale(x %*% v, center = FALSE, scale = d)
    return(list(u = u, d = d, v = v))
  } else {
    junk <- svd(x, LINPACK = TRUE)
    nc <- min(ncol(junk$u), n.components)
    return(list(u = junk$u[, 1:nc], d = junk$d[1:nc], v = junk$v[
      ,
      1:nc
    ]))
  }
}
quantitative.seq.func <- function(xresamp, y) {
  tt <- rep(0, dim(xresamp)[1])
  for (i in seq_len(dim(xresamp)[3])) {
    y.ranked <- rank(y, ties.method = "random") - (dim(xresamp)[2] +
      1) / 2
    tt <- tt + (xresamp[, , i] - (dim(xresamp)[2] + 1) / 2) %*%
      y.ranked
  }
  ns <- dim(xresamp)[2]
  tt <- tt / (dim(xresamp)[3] * (ns^3 - ns) / 12)
  return(list(tt = as.vector(tt), numer = as.vector(tt), sd = rep(
    1,
    length(tt)
  )))
}
multiclass.seq.func <- function(xresamp, y) {
  # number of classes and number of samples in each class
  K <- max(y)
  n.each <- rep(0, K)
  for (k in 1:K) {
    n.each[k] <- sum(y == k)
  }
  # the statistic
  tt <- temp <- rep(0, dim(xresamp)[1])
  stand.contrasts <- matrix(0, dim(xresamp)[1], K)

  for (i in seq_len(dim(xresamp)[3])) {
    for (k in 1:K) {
      temp <- rowSums(xresamp[, y == k, i])
      tt <- tt + temp^2 / n.each[k]
      stand.contrasts[, k] <- stand.contrasts[, k] + temp
    }
  }
  # finalize
  nresamp <- dim(xresamp)[3]
  ns <- dim(xresamp)[2]
  tt <- tt / nresamp * 12 / ns / (ns + 1) - 3 * (ns + 1)
  stand.contrasts <- stand.contrasts / nresamp
  stand.contrasts <- scale(stand.contrasts,
    center = n.each * (ns + 1) / 2,
    scale = sqrt(n.each * (ns - n.each) * (ns + 1) / 12)
  )
  return(list(
    tt = tt, numer = tt, sd = rep(1, length(tt)),
    stand.contrasts = stand.contrasts
  ))
}

# ==============================================================================
# samr.compute.delta.table
# ==============================================================================
## Jun added starts
samr.compute.delta.table <- function(
    samr.obj, min.foldchange = 0,
    dels = NULL, nvals = 50) {
  res <- NULL
  if (samr.obj$assay.type == "array") {
    res <- samr.compute.delta.table.array(
      samr.obj, min.foldchange,
      dels, nvals
    )
  } else if (samr.obj$assay.type == "seq") {
    res <- samr.compute.delta.table.seq(
      samr.obj, min.foldchange,
      dels
    )
  }
  return(res)
}
## Jun added ends

## Jun added the first row below, and commented the row
#   after it
samr.compute.delta.table.array <- function(
    samr.obj,
    min.foldchange = 0, dels = NULL, nvals = 50) {
  # samr.compute.delta.table <- function(samr.obj,
  #   min.foldchange=0, dels=NULL, nvals=50) {
  # computes delta table, starting with samr object 'a', for
  #   nvals values of delta
  lmax <- sqrt(max(abs(sort(samr.obj$tt) - samr.obj$evo)))
  if (is.null(dels)) {
    dels <- (seq(0, lmax, length = nvals)^2)
  }
  res1 <- NULL
  foldchange.cond.up <- matrix(TRUE,
    nrow = nrow(samr.obj$ttstar),
    ncol = ncol(samr.obj$ttstar)
  )
  foldchange.cond.lo <- matrix(TRUE,
    nrow = nrow(samr.obj$ttstar),
    ncol = ncol(samr.obj$ttstar)
  )
  if (!is.null(samr.obj$foldchange[1]) & (min.foldchange >
    0)) {
    foldchange.cond.up <- samr.obj$foldchange.star >= min.foldchange
    foldchange.cond.lo <- samr.obj$foldchange.star <= 1 / min.foldchange
  }
  cutup <- rep(NA, length(dels))
  cutlow <- rep(NA, length(dels))
  g2 <- rep(NA, length(dels))
  errup <- matrix(NA, ncol = length(dels), nrow = ncol(samr.obj$ttstar0))
  errlow <- matrix(NA, ncol = length(dels), nrow = ncol(samr.obj$ttstar0))
  cat("", fill = TRUE)
  cat("Computing delta table", fill = TRUE)
  for (ii in seq_len(length(dels))) {
    cat(ii, fill = TRUE)
    ttt <- detec.slab(samr.obj, dels[ii], min.foldchange)
    cutup[ii] <- 1e+10
    if (length(ttt$pup > 0)) {
      cutup[ii] <- min(samr.obj$tt[ttt$pup])
    }
    cutlow[ii] <- -1e+10
    if (length(ttt$plow) > 0) {
      cutlow[ii] <- max(samr.obj$tt[ttt$plow])
    }
    g2[ii] <- sumlengths(ttt)
    errup[, ii] <- colSums(samr.obj$ttstar0 > cutup[ii] &
      foldchange.cond.up)
    errlow[, ii] <- colSums(samr.obj$ttstar0 < cutlow[ii] &
      foldchange.cond.lo)
  }
  gmed <- apply(errup + errlow, 2, median)
  g90 <- apply(errup + errlow, 2, quantile, 0.9)
  res1 <- cbind(
    samr.obj$pi0 * gmed, samr.obj$pi0 * g90, g2,
    samr.obj$pi0 * gmed / g2, samr.obj$pi0 * g90 / g2, cutlow,
    cutup
  )
  res1 <- cbind(dels, res1)
  dimnames(res1) <- list(NULL, c(
    "delta", "# med false pos",
    "90th perc false pos", "# called", "median FDR", "90th perc FDR",
    "cutlo", "cuthi"
  ))
  return(res1)
}

#######################################################################
# \tcompute the delta table for sequencing data
#######################################################################
samr.compute.delta.table.seq <- function(
    samr.obj,
    min.foldchange = 0, dels = NULL) {
  res1 <- NULL
  flag <- TRUE
  ## check whether any gene satisfies the foldchange
  #   restrictions
  if ((samr.obj$resp.type == samr.const.twoclass.unpaired.response |
    samr.obj$resp.type == samr.const.twoclass.paired.response) &
    (min.foldchange > 0)) {
    sat.up <- (samr.obj$foldchange >= min.foldchange) & (samr.obj$evo >
      0)
    sat.dn <- (samr.obj$foldchange <= 1 / min.foldchange) &
      (samr.obj$evo < 0)
    if (sum(sat.up) + sum(sat.dn) == 0) {
      flag <- FALSE
    }
  }
  if (flag) {
    if (is.null(dels)) {
      dels <- generate.dels(samr.obj, min.foldchange = min.foldchange)
    }
    cat("Number of thresholds chosen (all possible thresholds) =",
      length(dels),
      fill = TRUE
    )
    if (length(dels) > 0) {
      ## sort delta to make the fast calculation right
      dels <- sort(dels)
      ## get the upper and lower cutoffs
      cat("Getting all the cutoffs for the thresholds...\n")
      slabs <- samr.seq.detec.slabs(samr.obj, dels, min.foldchange)
      cutup <- slabs$cutup
      cutlow <- slabs$cutlow
      g2 <- slabs$g2
      ## get the number of errors under the null hypothesis
      cat("Getting number of false positives in the permutation...\n")
      errnum <- samr.seq.null.err(
        samr.obj, min.foldchange,
        cutup, cutlow
      )
      res1 <- NULL
      gmed <- apply(errnum, 2, median)
      g90 <- apply(errnum, 2, quantile, 0.9)
      res1 <- cbind(samr.obj$pi0 * gmed, samr.obj$pi0 *
        g90, g2, samr.obj$pi0 * gmed / g2, samr.obj$pi0 *
        g90 / g2, cutlow, cutup)
      res1 <- cbind(dels, res1)
      dimnames(res1) <- list(NULL, c(
        "delta", "# med false pos",
        "90th perc false pos", "# called", "median FDR",
        "90th perc FDR", "cutlo", "cuthi"
      ))
    }
  }
  return(res1)
}

# ==============================================================================
# samr.plot
# ==============================================================================
samr.plot <- function(samr.obj, del = NULL, min.foldchange = 0) {
  ## make observed-expected plot
  ## takes foldchange into account too
  if (is.null(del)) {
    del <- sqrt(max(abs(sort(samr.obj$tt) - samr.obj$evo)))
  }
  b <- detec.slab(samr.obj, del, min.foldchange)
  foldchange.cond.up <- rep(TRUE, length(samr.obj$evo))
  foldchange.cond.lo <- rep(TRUE, length(samr.obj$evo))
  if (!is.null(samr.obj$foldchange[1]) & (min.foldchange >
    0)) {
    foldchange.cond.up <- samr.obj$foldchange >= min.foldchange
    foldchange.cond.lo <- samr.obj$foldchange <= 1 / min.foldchange
  }
  col <- rep(1, length(samr.obj$evo))
  col[b$plow] <- 3
  col[b$pup] <- 2
  if (!is.null(samr.obj$foldchange[1]) & (min.foldchange >
    0)) {
    col[!foldchange.cond.lo & !foldchange.cond.up] <- 1
  }
  col.ordered <- col[order(samr.obj$tt)]
  ylims <- range(samr.obj$tt)
  xlims <- range(samr.obj$evo)
  plot(samr.obj$evo, sort(samr.obj$tt),
    xlab = "expected score",
    ylab = "observed score", ylim = ylims, xlim = xlims,
    type = "n"
  )
  points(samr.obj$evo, sort(samr.obj$tt), col = col.ordered)
  abline(0, 1)
  abline(del, 1, lty = 2)
  abline(-del, 1, lty = 2)
}

# ==============================================================================
# samr.compute.siggenes.table
# ==============================================================================
samr.compute.siggenes.table <- function(
    samr.obj, del,
    data, delta.table, min.foldchange = 0, all.genes = FALSE,
    compute.localfdr = FALSE) {
  ## computes significant genes table, starting with samr
  #   object 'a' and 'delta.table'
  ##  for a  **single** value del
  ## if all.genes is true, all genes are printed (and value
  #   of del is ignored)
  if (is.null(data$geneid)) {
    data$geneid <- paste("g", seq_len(nrow(data$x)), sep = "")
  }
  if (is.null(data$genenames)) {
    data$genenames <- paste("g", seq_len(nrow(data$x)), sep = "")
  }
  if (!all.genes) {
    sig <- detec.slab(samr.obj, del, min.foldchange)
  }
  if (all.genes) {
    p <- length(samr.obj$tt)
    pup <- (1:p)[samr.obj$tt >= 0]
    plo <- (1:p)[samr.obj$tt < 0]
    sig <- list(pup = pup, plo = plo)
  }
  if (compute.localfdr) {
    aa <- localfdr(samr.obj, min.foldchange)
    if (length(sig$pup) > 0) {
      fdr.up <- predictlocalfdr(aa$smooth.object, samr.obj$tt[sig$pup])
    }
    if (length(sig$plo) > 0) {
      fdr.lo <- predictlocalfdr(aa$smooth.object, samr.obj$tt[sig$plo])
    }
  }
  qvalues <- NULL
  if (length(sig$pup) > 0 | length(sig$plo) > 0) {
    qvalues <- qvalue.func(samr.obj, sig, delta.table)
  }
  res.up <- NULL
  res.lo <- NULL
  done <- FALSE

  # two class unpaired or paired  (foldchange is reported)
  if ((samr.obj$resp.type == samr.const.twoclass.unpaired.response |
    samr.obj$resp.type == samr.const.twoclass.paired.response)) {
    if (!is.null(sig$pup)) {
      res.up <- cbind(
        sig$pup + 1, data$genenames[sig$pup],
        data$geneid[sig$pup], samr.obj$tt[sig$pup], samr.obj$numer[sig$pup],
        samr.obj$sd[sig$pup], samr.obj$foldchange[sig$pup],
        qvalues$qvalue.up
      )
      if (compute.localfdr) {
        res.up <- cbind(res.up, fdr.up)
      }
      temp.names <- list(NULL, c(
        "Row", "Gene ID", "Gene Name",
        "Score(d)", "Numerator(r)", "Denominator(s+s0)",
        "Fold Change", "q-value(%)"
      ))
      if (compute.localfdr) {
        temp.names[[2]] <- c(temp.names[[2]], "localfdr(%)")
      }
      dimnames(res.up) <- temp.names
    }
    if (!is.null(sig$plo)) {
      res.lo <- cbind(
        sig$plo + 1, data$genenames[sig$plo],
        data$geneid[sig$plo], samr.obj$tt[sig$plo], samr.obj$numer[sig$plo],
        samr.obj$sd[sig$plo], samr.obj$foldchange[sig$plo],
        qvalues$qvalue.lo
      )
      if (compute.localfdr) {
        res.lo <- cbind(res.lo, fdr.lo)
      }
      temp.names <- list(NULL, c(
        "Row", "Gene ID", "Gene Name",
        "Score(d)", "Numerator(r)", "Denominator(s+s0)",
        "Fold Change", "q-value(%)"
      ))
      if (compute.localfdr) {
        temp.names[[2]] <- c(temp.names[[2]], "localfdr(%)")
      }
      dimnames(res.lo) <- temp.names
    }
    done <- TRUE
  }

  # multiclass
  if (samr.obj$resp.type == samr.const.multiclass.response) {
    if (!is.null(sig$pup)) {
      res.up <- cbind(
        sig$pup + 1, data$genenames[sig$pup],
        data$geneid[sig$pup], samr.obj$tt[sig$pup], samr.obj$numer[sig$pup],
        samr.obj$sd[sig$pup], samr.obj$stand.contrasts[sig$pup, ],
        qvalues$qvalue.up
      )

      if (compute.localfdr) {
        res.up <- cbind(res.up, fdr.up)
      }

      collabs.contrast <- paste(
        "contrast-", as.character(seq_len(ncol(samr.obj$stand.contrasts))),
        sep = ""
      )
      temp.names <- list(NULL, c(
        "Row", "Gene ID", "Gene Name",
        "Score(d)", "Numerator(r)", "Denominator(s+s0)",
        collabs.contrast, "q-value(%)"
      ))

      if (compute.localfdr) {
        temp.names[[2]] <- c(temp.names[[2]], "localfdr(%)")
      }
      dimnames(res.up) <- temp.names
    }
    res.lo <- NULL
    done <- TRUE
  }

  # all other cases
  if (!done) {
    if (!is.null(sig$pup)) {
      res.up <- cbind(
        sig$pup + 1, data$genenames[sig$pup],
        data$geneid[sig$pup], samr.obj$tt[sig$pup], samr.obj$numer[sig$pup],
        samr.obj$sd[sig$pup], samr.obj$foldchange[sig$pup],
        qvalues$qvalue.up
      )
      if (compute.localfdr) {
        res.up <- cbind(res.up, fdr.up)
      }
      temp.names <- list(NULL, c(
        "Row", "Gene ID", "Gene Name",
        "Score(d)", "Numerator(r)", "Denominator(s+s0)",
        "q-value(%)"
      ))
      if (compute.localfdr) {
        temp.names[[2]] <- c(temp.names[[2]], "localfdr(%)")
      }
      dimnames(res.up) <- temp.names
    }
    if (!is.null(sig$plo)) {
      res.lo <- cbind(
        sig$plo + 1, data$genenames[sig$plo],
        data$geneid[sig$plo], samr.obj$tt[sig$plo], samr.obj$numer[sig$plo],
        samr.obj$sd[sig$plo], samr.obj$foldchange[sig$plo],
        qvalues$qvalue.lo
      )
      if (compute.localfdr) {
        res.lo <- cbind(res.lo, fdr.lo)
      }
      temp.names <- list(NULL, c(
        "Row", "Gene ID", "Gene Name",
        "Score(d)", "Numerator(r)", "Denominator(s+s0)",
        "q-value(%)"
      ))
      if (compute.localfdr) {
        temp.names[[2]] <- c(temp.names[[2]], "localfdr(%)")
      }
      dimnames(res.lo) <- temp.names
    }
    done <- TRUE
  }
  if (!is.null(res.up)) {
    o1 <- order(-samr.obj$tt[sig$pup])
    res.up <- res.up[o1, , drop = FALSE]
  }
  if (!is.null(res.lo)) {
    o2 <- order(samr.obj$tt[sig$plo])
    res.lo <- res.lo[o2, , drop = FALSE]
  }
  color.ind.for.multi <- NULL
  if (
    samr.obj$resp.type == samr.const.multiclass.response & !is.null(sig$pup)
  ) {
    condition_1 <-
      samr.obj$stand.contrasts[sig$pup, ] > samr.obj$stand.contrasts.95[2]
    condition_2 <-
      samr.obj$stand.contrasts[sig$pup, ] < samr.obj$stand.contrasts.95[1]
    color.ind.for.multi <- 1 * condition_1 + (-1) * condition_2
  }
  ngenes.up <- nrow(res.up)
  if (is.null(ngenes.up)) {
    ngenes.up <- 0
  }
  ngenes.lo <- nrow(res.lo)
  if (is.null(ngenes.lo)) {
    ngenes.lo <- 0
  }
  return(
    list(
      genes.up = res.up, genes.lo = res.lo,
      color.ind.for.multi = color.ind.for.multi, ngenes.up = ngenes.up,
      ngenes.lo = ngenes.lo
    )
  )
}
generate.dels <- function(samr.obj, min.foldchange = 0) {
  dels <- NULL
  ## initialize calculation
  tag <- order(samr.obj$tt)
  if ((samr.obj$resp.type == samr.const.twoclass.unpaired.response |
    samr.obj$resp.type == samr.const.twoclass.paired.response) &
    (min.foldchange > 0)) {
    res.mat <- data.frame(
      tt = samr.obj$tt[tag], fc = samr.obj$foldchange[tag],
      evo = samr.obj$evo, dif = samr.obj$tt[tag] - samr.obj$evo
    )
    res.up <- res.mat[res.mat$evo > 0, ]
    res.lo <- res.mat[res.mat$evo < 0, ]
    res.up <- res.up[res.up$fc >= min.foldchange, ]
    res.lo <- res.lo[res.lo$fc <= 1 / min.foldchange, ]
  } else {
    res.mat <- data.frame(
      tt = samr.obj$tt[tag], evo = samr.obj$evo,
      dif = samr.obj$tt[tag] - samr.obj$evo
    )
    res.up <- res.mat[res.mat$evo > 0, ]
    res.lo <- res.mat[res.mat$evo < 0, ]
  }
  ## for the upper part
  up.vec <- rep(NA, nrow(res.up))
  if (nrow(res.up) > 0) {
    st <- 1e-08
    i.cur <- 1
    for (i in seq_len(nrow(res.up))) {
      if (res.up$dif[i] > st) {
        st <- res.up$dif[i]
        up.vec[i.cur] <- st
        i.cur <- i.cur + 1
      }
    }
  }
  ## for the lower part
  lo.vec <- rep(NA, nrow(res.lo))
  if (nrow(res.lo) > 0) {
    st <- -1e-08
    i.cur <- 1
    for (i in seq(nrow(res.lo), 1)) {
      if (res.lo$dif[i] < st) {
        st <- res.lo$dif[i]
        lo.vec[i.cur] <- st
        i.cur <- i.cur + 1
      }
    }
  }
  ## combine them
  vec <- c(up.vec, -lo.vec)
  vec <- vec[!is.na(vec)]
  vec <- vec - 1e-08
  dels <- sort(unique(vec))
  return(dels)
}
samr.seq.detec.slabs <- function(samr.obj, dels, min.foldchange) {
  ## initialize calculation
  tag <- order(samr.obj$tt)
  if ((samr.obj$resp.type == samr.const.twoclass.unpaired.response |
    samr.obj$resp.type == samr.const.twoclass.paired.response) &
    (min.foldchange > 0)) {
    res.mat <- data.frame(
      tt = samr.obj$tt[tag], fc = samr.obj$foldchange[tag],
      evo = samr.obj$evo, dif = samr.obj$tt[tag] - samr.obj$evo
    )
    res.up <- res.mat[res.mat$evo > 0, ]
    res.lo <- res.mat[res.mat$evo < 0, ]
    res.up <- res.up[res.up$fc >= min.foldchange, ]
    res.lo <- res.lo[res.lo$fc <= 1 / min.foldchange, ]
  } else {
    res.mat <- data.frame(
      tt = samr.obj$tt[tag], evo = samr.obj$evo,
      dif = samr.obj$tt[tag] - samr.obj$evo
    )
    res.up <- res.mat[res.mat$evo > 0, ]
    res.lo <- res.mat[res.mat$evo < 0, ]
  }
  ## begin calculating
  cutup <- rep(1e+10, length(dels))
  cutlow <- rep(-1e+10, length(dels))
  g2.up <- g2.lo <- rep(0, length(dels))
  if (nrow(res.up) > 0) {
    res.up <- data.frame(
      dif = res.up$dif, tt = res.up$tt,
      num = seq(nrow(res.up), 1)
    )
    ## get the upper part
    j <- 1
    ii <- 1
    while (j <= nrow(res.up) & ii <= length(dels)) {
      if (res.up$dif[j] > dels[ii]) {
        cutup[ii] <- res.up$tt[j]
        g2.up[ii] <- res.up$num[j]
        ii <- ii + 1
      } else {
        j <- j + 1
      }
    }
  }
  if (nrow(res.lo) > 0) {
    res.lo <- data.frame(
      dif = res.lo$dif, tt = res.lo$tt,
      num = seq_len(nrow(res.lo))
    )
    ## get the lower part
    j <- nrow(res.lo)
    ii <- 1
    while (j >= 1 & ii <= length(dels)) {
      if (res.lo$dif[j] < -dels[ii]) {
        cutlow[ii] <- res.lo$tt[j]
        g2.lo[ii] <- res.lo$num[j]
        ii <- ii + 1
      } else {
        j <- j - 1
      }
    }
  }
  g2 <- g2.up + g2.lo
  return(list(cutup = cutup, cutlow = cutlow, g2 = g2))
}
sumlengths <- function(aa) {
  length(aa$pl) + length(aa$pu)
}

samr.seq.null.err <- function(
    samr.obj, min.foldchange,
    cutup, cutlow) {
  errup <- matrix(NA, ncol = length(cutup), nrow = ncol(samr.obj$ttstar0))
  errlow <- matrix(NA, ncol = length(cutlow), nrow = ncol(samr.obj$ttstar0))
  cutup.rank <- rank(cutup, ties.method = "min")
  cutlow.rank <- rank(-cutlow, ties.method = "min")
  for (jj in seq_len(ncol(samr.obj$ttstar0))) {
    keep.up <- keep.dn <- samr.obj$ttstar0[, jj]
    if ((samr.obj$resp.type == samr.const.twoclass.unpaired.response |
      samr.obj$resp.type == samr.const.twoclass.paired.response) &
      (min.foldchange > 0)) {
      keep.up <- keep.up[samr.obj$foldchange.star[, jj] >=
        min.foldchange]
      keep.dn <- keep.dn[samr.obj$foldchange.star[, jj] <=
        1 / min.foldchange]
    }
    errup[jj, ] <- length(keep.up) - (rank(c(cutup, keep.up),
      ties.method = "min"
    )[seq_len(length(cutup))] - cutup.rank)
    errlow[jj, ] <- length(keep.dn) - (rank(c(-cutlow, -keep.dn),
      ties.method = "min"
    )[seq_len(length(cutlow))] - cutlow.rank)
  }
  errnum <- errup + errlow
  return(errnum)
}
detec.slab <- function(samr.obj, del, min.foldchange) {
  ## find genes above and below the slab of half-width del
  # this calculation is tricky- for consistency, the slab
  #   condition picks
  # all genes that are beyond the first departure from the
  #   slab
  # then the fold change condition is applied (if applicable)
  n <- length(samr.obj$tt)
  tt <- samr.obj$tt
  evo <- samr.obj$evo
  tag <- order(tt)
  pup <- NULL
  foldchange.cond.up <- rep(TRUE, length(evo))
  foldchange.cond.lo <- rep(TRUE, length(evo))
  if (!is.null(samr.obj$foldchange[1]) & (min.foldchange >
    0)) {
    foldchange.cond.up <- samr.obj$foldchange >= min.foldchange
    foldchange.cond.lo <- samr.obj$foldchange <= 1 / min.foldchange
  }
  o1 <- (1:n)[(tt[tag] - evo > del) & evo > 0]
  if (length(o1) > 0) {
    o1 <- o1[1]
    o11 <- o1:n
    o111 <- rep(FALSE, n)
    o111[tag][o11] <- TRUE
    pup <- (1:n)[o111 & foldchange.cond.up]
  }
  plow <- NULL
  o2 <- (1:n)[(evo - tt[tag] > del) & evo < 0]
  if (length(o2) > 0) {
    o2 <- o2[length(o2)]
    o22 <- 1:o2
    o222 <- rep(FALSE, n)
    o222[tag][o22] <- TRUE
    plow <- (1:n)[o222 & foldchange.cond.lo]
  }
  return(list(plow = plow, pup = pup))
}

#' @importFrom stats smooth.spline
localfdr <- function(
    samr.obj, min.foldchange, perc = 0.01,
    df = 10) {
  ## estimates compute.localfdr at score 'd', using SAM
  #   object 'samr.obj'
  ## 'd' can be a vector of d scores
  ## returns estimate of symmetric fdr  as a percentage
  # this version uses a 1% symmetric window, and does not
  #   estimate fdr in
  # windows  having fewer than 100 genes
  ## to use: first run samr and then pass the resulting fit
  #   object to
  ## localfdr
  ## NOTE: at most 20 of the perms are used to estimate the
  #   fdr (for speed sake)
  # I try two window shapes: symmetric and an assymetric one
  # currently I use the symmetric window to estimate the
  #   compute.localfdr
  ngenes <- length(samr.obj$tt)
  mingenes <- 50
  # perc is increased, in order to get at least mingenes in a
  #   window
  perc <- max(perc, mingenes / length(samr.obj$tt))
  nperms.to.use <- min(20, ncol(samr.obj$ttstar))
  nperms <- ncol(samr.obj$ttstar)
  d <- seq(sort(samr.obj$tt)[1], sort(samr.obj$tt)[ngenes],
    length = 100
  )
  ndscore <- length(d)
  dvector <- rep(NA, ndscore)
  ind.foldchange <- rep(TRUE, length(samr.obj$tt))
  if (!is.null(samr.obj$foldchange[1]) & min.foldchange > 0) {
    ind.foldchange <- (samr.obj$foldchange >= min.foldchange) |
      (samr.obj$foldchange <= min.foldchange)
  }
  fdr.temp <- function(temp, dlow, dup, pi0, ind.foldchange) {
    return(sum(pi0 * (temp >= dlow & temp <= dup & ind.foldchange)))
  }
  for (i in 1:ndscore) {
    pi0 <- samr.obj$pi0
    r <- sum(samr.obj$tt < d[i])
    r22 <- round(max(r - length(samr.obj$tt) * perc / 2, 1))
    dlow.sym <- sort(samr.obj$tt)[r22]
    r22 <- min(r + length(samr.obj$tt) * perc / 2, length(samr.obj$tt))
    dup.sym <- sort(samr.obj$tt)[r22]
    oo <- samr.obj$tt >= dlow.sym & samr.obj$tt <= dup.sym &
      ind.foldchange
    if (!is.null(samr.obj$foldchange[1]) & min.foldchange >
      0) {
      temp <- as.vector(samr.obj$foldchange.star[, 1:nperms.to.use])
      ind.foldchange <- (temp >= min.foldchange) | (temp <=
        min.foldchange)
    }
    temp <- samr.obj$ttstar0[, sample(1:nperms, size = nperms.to.use)]
    fdr.sym <- median(apply(
      temp, 2, fdr.temp, dlow.sym,
      dup.sym, pi0, ind.foldchange
    ))
    fdr.sym <- 100 * fdr.sym / sum(oo)
    dlow.sym <- dlow.sym
    dup.sym <- dup.sym
    dvector[i] <- fdr.sym
  }
  om <- !is.na(dvector) & (dvector != Inf)
  aa <- smooth.spline(d[om], dvector[om], df = df)
  return(list(smooth.object = aa, perc = perc, df = df))
}

predictlocalfdr <- function(smooth.object, d) {
  yhat <- predict(smooth.object, d)$y
  yhat <- pmin(yhat, 100)
  yhat <- pmax(yhat, 0)
  return(yhat)
}

qvalue.func <- function(samr.obj, sig, delta.table) {
  # returns q-value as a percentage (out of 100)
  LARGE <- 1e+10
  qvalue.up <- rep(NA, length(sig$pup))
  o1 <- sig$pup
  cutup <- delta.table[, 8]
  FDR <- delta.table[, 5]
  ii <- 0
  for (i in o1) {
    o <- abs(cutup - samr.obj$tt[i])
    o[is.na(o)] <- LARGE
    oo <- seq_len(length(o))[o == min(o)]
    oo <- oo[length(oo)]
    ii <- ii + 1
    qvalue.up[ii] <- FDR[oo]
  }
  qvalue.lo <- rep(NA, length(sig$plo))
  o2 <- sig$plo
  cutlo <- delta.table[, 7]
  ii <- 0
  for (i in o2) {
    o <- abs(cutlo - samr.obj$tt[i])
    o[is.na(o)] <- LARGE
    oo <- seq_len(length(o))[o == min(o)]
    oo <- oo[length(oo)]
    ii <- ii + 1
    qvalue.lo[ii] <- FDR[oo]
  }
  # any qvalues that are missing, are set to 1 (the highest
  #   value)
  qvalue.lo[is.na(qvalue.lo)] <- 1
  qvalue.up[is.na(qvalue.up)] <- 1
  # ensure that each qvalue vector is monotone non-increasing
  o1 <- order(samr.obj$tt[sig$plo])
  qv1 <- qvalue.lo[o1]
  qv11 <- qv1
  if (length(qv1) > 1) {
    for (i in 2:length(qv1)) {
      if (qv11[i] < qv11[i - 1]) {
        qv11[i] <- qv11[i - 1]
      }
    }
    qv111 <- qv11
    qv111[o1] <- qv11
  } else {
    qv111 <- qv1
  }
  o2 <- order(samr.obj$tt[sig$pup])
  qv2 <- qvalue.up[o2]
  qv22 <- qv2
  if (length(qv2) > 1) {
    for (i in 2:length(qv2)) {
      if (qv22[i] > qv22[i - 1]) {
        qv22[i] <- qv22[i - 1]
      }
    }
    qv222 <- qv22
    qv222[o2] <- qv22
  } else {
    qv222 <- qv2
  }
  return(list(qvalue.lo = 100 * qv111, qvalue.up = 100 * qv222))
}
