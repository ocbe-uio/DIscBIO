# This script contains customized versions of functions found in the samr package. This is necessary because samr seems to have been abandoned, so an upstream collaboration doesn't seem possible at the time of writing.

# ==============================================================================
# Constants
# ==============================================================================
samr.const.twoclass.unpaired.response <- "Two class unpaired"
samr.const.twoclass.paired.response <- "Two class paired"
samr.const.oneclass.response <- "One class"
samr.const.quantitative.response <- "Quantitative"
samr.const.multiclass.response <- "Multiclass"
samr.const.twoclass.unpaired.timecourse.response <- "Two class unpaired timecourse"
samr.const.twoclass.paired.timecourse.response <- "Two class paired timecourse"
samr.const.oneclass.timecourse.response <- "One class timecourse"
samr.const.survival.response <- "Survival"
samr.const.patterndiscovery.response <- "Pattern discovery"

# ==============================================================================
# Functions
# ==============================================================================

#' @title Significance analysis of microarrays
#' @description This function is an adaptation of `samr::samr`
sammy <- function (data, resp.type = c("Quantitative", "Two class unpaired",
    "Survival", "Multiclass", "One class", "Two class paired",
    "Two class unpaired timecourse", "One class timecourse",
    "Two class paired timecourse", "Pattern discovery"), assay.type = c("array",
    "seq"), s0 = NULL, s0.perc = NULL, nperms = 100, center.arrays = FALSE,
    testStatistic = c("standard", "wilcoxon"), time.summary.type = c("slope",
        "signed.area"), regression.method = c("standard", "ranks"),
    return.x = FALSE, knn.neighbors = 10, random.seed = NULL,
    nresamp = 20, nresamp.perm = NULL, xl.mode = c("regular",
        "firsttime", "next20", "lasttime"), xl.time = NULL, xl.prevfit = NULL)
{
    this.call = match.call()
    resp.type.arg = match.arg(resp.type)
    assay.type = match.arg(assay.type)
    xl.mode = match.arg(xl.mode)
    if (!is.null(random.seed)) {
        set.seed(random.seed)
    }
    if (is.null(nresamp.perm)) {
        nresamp.perm = nresamp
    }
    nresamp.perm = min(nresamp, nresamp.perm)
    if (xl.mode == "regular" | xl.mode == "firsttime") {
        x = NULL
        xresamp = NULL
        ttstar0 = NULL
        evo = NULL
        ystar = NULL
        sdstar.keep = NULL
        censoring.status = NULL
        sdstar = NULL
        pi0 = NULL
        stand.contrasts = NULL
        stand.contrasts.star = NULL
        stand.contrasts.95 = NULL
        foldchange = NULL
        foldchange.star = NULL
        perms = NULL
        permsy = NULL
        eigengene = NULL
        eigengene.number = NULL
        testStatistic <- match.arg(testStatistic)
        time.summary.type <- match.arg(time.summary.type)
        regression.method <- match.arg(regression.method)
        x = data$x
        y = data$y
        argy = y
        if (!is.null(data$eigengene.number)) {
            eigengene.number = data$eigengene.number
        }
        if (sum(is.na(x)) > 0) {
            x = impute.knn(x, k = knn.neighbors)
            if (!is.matrix(x)) {
                x = x$data
            }
        }
        are.blocks.specified = FALSE
        cond = (resp.type == "One class") | (resp.type == "Two class unpaired timecourse") |
            (resp.type == "One class unpaired timecourse") |
            (resp.type == "Two class paired timecourse") | (resp.type ==
            "Pattern discovery")
        if (assay.type == "seq" & cond) {
            stop(paste("Resp.type=", resp.type, " not allowed when assay.type='seq'"))
        }
        if (assay.type == "seq" & min(x) < 0) {
            stop(paste("Negative values not allowed when assay.type='seq'"))
        }
        if (assay.type == "seq" & (sum(x%%1 != 0) != 0)) {
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
        depth = scaling.factors = rep(NA, ncol(x))
        scaling.factors = (prod(depth)^(1/length(depth)))/depth
        if (assay.type == "seq") {
            cat("Estimating sequencing depths...", fill = T)
            depth = samr.estimate.depth(x)
            cat("Resampling to get new data matrices...", fill = T)
            xresamp = resa(x, depth, nresamp = nresamp)
        }
        scaling.factors = (prod(depth)^(1/length(depth)))/depth
        if (resp.type == samr.const.twoclass.unpaired.response) {
            if (substring(y[1], 2, 6) == "Block" | substring(y[1],
                2, 6) == "block") {
                junk = parse.block.labels.for.2classes(y)
                y = junk$y
                blocky = junk$blocky
                are.blocks.specified = TRUE
            }
        }
        if (resp.type == samr.const.twoclass.unpaired.response |
            resp.type == samr.const.twoclass.paired.response |
            resp.type == samr.const.oneclass.response | resp.type ==
            samr.const.quantitative.response | resp.type == samr.const.multiclass.response) {
            y = as.numeric(y)
        }
        sd.internal = NULL
        if (resp.type == samr.const.twoclass.unpaired.timecourse.response |
            resp.type == samr.const.twoclass.paired.timecourse.response |
            resp.type == samr.const.oneclass.timecourse.response) {
            junk = parse.time.labels.and.summarize.data(x, y,
                resp.type, time.summary.type)
            y = junk$y
            x = junk$x
            sd.internal = sqrt(rowMeans(junk$sd^2))
            if (min(table(y)) == 1) {
                cat("", fill = T)
                cat("Warning: only one timecourse in one or more classes;\nSAM plot and FDRs will be unreliable; only gene scores are informative",
                  fill = T)
            }
        }
        if (resp.type == samr.const.twoclass.unpaired.timecourse.response) {
            resp.type = samr.const.twoclass.unpaired.response
        }
        if (resp.type == samr.const.twoclass.paired.timecourse.response) {
            resp.type = samr.const.twoclass.paired.response
        }
        if (resp.type == samr.const.oneclass.timecourse.response) {
            resp.type = samr.const.oneclass.response
        }
        stand.contrasts = NULL
        stand.contrasts.95 = NULL
        if (resp.type == samr.const.survival.response) {
            censoring.status = data$censoring.status
        }
        check.format(y, resp.type = resp.type, censoring.status = censoring.status)
        if (resp.type == samr.const.quantitative.response & regression.method ==
            "ranks") {
            y = rank(y)
            x = t(apply(x, 1, rank))
        }
        n <- nrow(x)
        ny <- length(y)
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
            s0 = quantile(sd, 0.05)
            s0.perc = 0.05
        }
        if (is.null(s0) & assay.type == "array") {
            if (!is.null(s0.perc)) {
                if ((s0.perc != -1 & s0.perc < 0) | s0.perc >
                  100) {
                  stop("Illegal value for s0.perc: must be between 0 and 100, or equal\nto (-1) (meaning that s0 should be set to zero)")
                }
                if (s0.perc == -1) {
                  s0 = 0
                }
                if (s0.perc >= 0) {
                  s0 <- quantile(init.fit$sd, s0.perc/100)
                }
            }
            if (is.null(s0.perc)) {
                s0 = est.s0(init.fit$tt, init.fit$sd)$s0.hat
                s0.perc = 100 * sum(init.fit$sd < s0)/length(init.fit$sd)
            }
        }
        if (assay.type == "seq") {
            s0 = 0
            s0.perc = 0
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
            tt = junk2$tt
            stand.contrasts = junk2$stand.contrasts
        }
        if (resp.type == samr.const.quantitative.response & assay.type ==
            "array") {
            tt <- quantitative.func(x, y, s0 = s0)$tt
        }
        if (resp.type == samr.const.patterndiscovery.response &
            assay.type == "array") {
            junk <- patterndiscovery.func(x, s0 = s0, eigengene.number = eigengene.number)
            tt <- junk$tt
            eigengene = junk$eigengene
        }
        if (resp.type == samr.const.twoclass.unpaired.response &
            assay.type == "seq") {
            junk = wilcoxon.unpaired.seq.func(xresamp, y)
            tt = junk$tt
            numer = junk$numer
            sd = junk$sd
        }
        if (resp.type == samr.const.twoclass.paired.response &
            assay.type == "seq") {
            junk <- wilcoxon.paired.seq.func(xresamp, y)
            tt = junk$tt
            numer = junk$numer
            sd = junk$sd
        }
        if (resp.type == samr.const.quantitative.response & assay.type ==
            "seq") {
            junk <- quantitative.seq.func(xresamp, y)
            tt = junk$tt
            numer = junk$numer
            sd = junk$sd
        }
        if (resp.type == samr.const.survival.response & assay.type ==
            "seq") {
            junk <- cox.seq.func(xresamp, y, censoring.status)
            tt = junk$tt
            numer = junk$numer
            sd = junk$sd
        }
        if (resp.type == samr.const.multiclass.response & assay.type ==
            "seq") {
            junk2 <- multiclass.seq.func(xresamp, y)
            tt = junk2$tt
            numer = junk2$numer
            sd = junk2$sd
            stand.contrasts = junk2$stand.contrasts
        }
        if (resp.type == samr.const.quantitative.response | resp.type ==
            samr.const.multiclass.response | resp.type == samr.const.survival.response) {
            junk <- getperms(y, nperms)
            perms = junk$perms
            all.perms.flag = junk$all.perms.flag
            nperms.act = junk$nperms.act
        }
        if (resp.type == samr.const.twoclass.unpaired.response) {
            if (are.blocks.specified) {
                junk = compute.block.perms(y, blocky, nperms)
                permsy = matrix(junk$permsy, ncol = length(y))
                all.perms.flag = junk$all.perms.flag
                nperms.act = junk$nperms.act
            }
            else {
                junk <- getperms(y, nperms)
                permsy = matrix(y[junk$perms], ncol = length(y))
                all.perms.flag = junk$all.perms.flag
                nperms.act = junk$nperms.act
            }
        }
        if (resp.type == samr.const.oneclass.response) {
            if ((length(y) * log(2)) < log(nperms)) {
                allii = 0:((2^length(y)) - 1)
                nperms.act = 2^length(y)
                all.perms.flag = 1
            }
            else {
                nperms.act = nperms
                all.perms.flag = 0
            }
            permsy = matrix(NA, nrow = nperms.act, ncol = length(y))
            if (all.perms.flag == 1) {
                k = 0
                for (i in allii) {
                  junk = integer.base.b(i, b = 2)
                  if (length(junk) < length(y)) {
                    junk = c(rep(0, length(y) - length(junk)),
                      junk)
                  }
                  k = k + 1
                  permsy[k, ] = y * (2 * junk - 1)
                }
            }
            else {
                for (i in 1:nperms.act) {
                  permsy[i, ] = sample(c(-1, 1), size = length(y),
                    replace = TRUE)
                }
            }
        }
        if (resp.type == samr.const.twoclass.paired.response) {
            junk = compute.block.perms(y, abs(y), nperms)
            permsy = junk$permsy
            all.perms.flag = junk$all.perms.flag
            nperms.act = junk$nperms.act
        }
        if (resp.type == samr.const.patterndiscovery.response) {
            nperms.act = nperms
            perms = NULL
            permsy = NULL
            all.perms.flag = FALSE
        }
        sdstar.keep <- NULL
        if (assay.type != "seq") {
            sdstar.keep <- matrix(0, ncol = nperms.act, nrow = nrow(x))
        }
        ttstar <- matrix(0, nrow = nrow(x), ncol = nperms.act)
        foldchange.star = NULL
        if (resp.type == samr.const.twoclass.unpaired.response |
            resp.type == samr.const.twoclass.paired.response) {
            foldchange.star <- matrix(0, nrow = nrow(x), ncol = nperms.act)
        }
        if (resp.type == samr.const.multiclass.response) {
            stand.contrasts.star = array(NA, c(nrow(x), length(table(y)),
                nperms.act))
        }
    }
    if (xl.mode == "next20" | xl.mode == "lasttime") {
        evo = xl.prevfit$evo
        tt = xl.prevfit$tt
        numer = xl.prevfit$numer
        eigengene = xl.prevfit$eigengene
        eigengene.number = xl.prevfit$eigengene.number
        sd = xl.prevfit$sd - xl.prevfit$s0
        sd.internal = xl.prevfit$sd.internal
        ttstar = xl.prevfit$ttstar
        ttstar0 = xl.prevfit$ttstar0
        n = xl.prevfit$n
        pi0 = xl.prevfit$pi0
        foldchange = xl.prevfit$foldchange
        y = xl.prevfit$y
        x = xl.prevfit$x
        xresamp = xl.prevfit$xresamp
        censoring.status = xl.prevfit$censoring.status
        argy = xl.prevfit$argy
        testStatistic = xl.prevfit$testStatistic
        foldchange.star = xl.prevfit$foldchange.star
        s0 = xl.prevfit$s0
        s0.perc = xl.prevfit$s0.perc
        resp.type = xl.prevfit$resp.type
        resp.type.arg = xl.prevfit$resp.type.arg
        assay.type = xl.prevfit$assay.type
        sdstar.keep = xl.prevfit$sdstar.keep
        resp.type = xl.prevfit$resp.type
        stand.contrasts = xl.prevfit$stand.contrasts
        stand.contrasts.star = xl.prevfit$stand.contrasts.star
        stand.contrasts.95 = xl.prevfit$stand.contrasts.95
        perms = xl.prevfit$perms
        permsy = xl.prevfit$permsy
        nperms = xl.prevfit$nperms
        nperms.act = xl.prevfit$nperms.act
        all.perms.flag = xl.prevfit$all.perms.flag
        depth = xl.prevfit$depth
        scaling.factors = xl.prevfit$scaling.factors
        nresamp = xl.prevfit$nresamp
        nresamp.perm = xl.prevfit$nresamp.perm
    }
    if (xl.mode == "regular") {
        first = 1
        last = nperms.act
    }
    if (xl.mode == "firsttime") {
        first = 1
        last = 1
    }
    if (xl.mode == "next20") {
        first = xl.time
        last = min(xl.time + 19, nperms.act - 1)
    }
    if (xl.mode == "lasttime") {
        first = nperms.act
        last = nperms.act
    }
    for (b in first:last) {
        cat(c("perm=", b), fill = TRUE)
        if (assay.type == "array") {
            xstar <- x
        }
        if (assay.type == "seq") {
            xstar <- xresamp[, , 1:nresamp.perm]
        }
        if (resp.type == samr.const.oneclass.response) {
            ystar = permsy[b, ]
            if (testStatistic == "standard") {
                ttstar[, b] <- onesample.ttest.func(xstar, ystar,
                  s0 = s0, sd = sd.internal)$tt
            }
        }
        if (resp.type == samr.const.twoclass.paired.response) {
            ystar = permsy[b, ]
            if (assay.type == "array") {
                ttstar[, b] <- paired.ttest.func(xstar, ystar,
                  s0 = s0, sd = sd.internal)$tt
                foldchange.star[, b] = foldchange.paired(xstar,
                  ystar, data$logged2)
            }
            if (assay.type == "seq") {
                ttstar[, b] <- wilcoxon.paired.seq.func(xstar,
                  ystar)$tt
                foldchange.star[, b] <- foldchange.seq.twoclass.paired(x,
                  ystar, depth)
            }
        }
        if (resp.type == samr.const.twoclass.unpaired.response) {
            ystar = permsy[b, ]
            if (assay.type == "array") {
                if (testStatistic == "standard") {
                  junk <- ttest.func(xstar, ystar, s0 = s0, sd = sd.internal)
                }
                if (testStatistic == "wilcoxon") {
                  junk <- wilcoxon.func(xstar, ystar, s0 = s0)
                }
                ttstar[, b] <- junk$tt
                sdstar.keep[, b] <- junk$sd
                foldchange.star[, b] = foldchange.twoclass(xstar,
                  ystar, data$logged2)
            }
            if (assay.type == "seq") {
                ttstar[, b] <- wilcoxon.unpaired.seq.func(xstar,
                  ystar)$tt
                foldchange.star[, b] <- foldchange.seq.twoclass.unpaired(x,
                  ystar, depth)
            }
        }
        if (resp.type == samr.const.survival.response) {
            o <- perms[b, ]
            if (assay.type == "array") {
                ttstar[, b] <- cox.func(xstar, y[o], censoring.status = censoring.status[o],
                  s0 = s0)$tt
            }
            if (assay.type == "seq") {
                ttstar[, b] <- cox.seq.func(xstar, y[o], censoring.status = censoring.status[o])$tt
            }
        }
        if (resp.type == samr.const.multiclass.response) {
            ystar = y[perms[b, ]]
            if (assay.type == "array") {
                junk <- multiclass.func(xstar, ystar, s0 = s0)
                ttstar[, b] <- junk$tt
                sdstar.keep[, b] <- junk$sd
                stand.contrasts.star[, , b] = junk$stand.contrasts
            }
            if (assay.type == "seq") {
                junk <- multiclass.seq.func(xstar, ystar)
                ttstar[, b] <- junk$tt
                stand.contrasts.star[, , b] <- junk$stand.contrasts
            }
        }
        if (resp.type == samr.const.quantitative.response) {
            ystar = y[perms[b, ]]
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
            xstar = permute.rows(x)
            junk <- patterndiscovery.func(xstar, s0 = s0, eigengene.number = eigengene.number)
            ttstar[, b] <- junk$tt
            sdstar.keep[, b] <- junk$sd
        }
    }
    if (xl.mode == "regular" | xl.mode == "lasttime") {
        ttstar0 <- ttstar
        for (j in 1:ncol(ttstar)) {
            ttstar[, j] <- -1 * sort(-1 * ttstar[, j])
        }
        for (i in 1:nrow(ttstar)) {
            ttstar[i, ] <- sort(ttstar[i, ])
        }
        evo <- apply(ttstar, 1, mean)
        evo <- evo[length(evo):1]
        sdstar <- sdstar.keep
        pi0 = 1
        if (resp.type != samr.const.multiclass.response) {
            qq <- quantile(ttstar, c(0.25, 0.75))
        }
        if (resp.type == samr.const.multiclass.response) {
            qq <- quantile(ttstar, c(0, 0.5))
        }
        pi0 <- sum(tt > qq[1] & tt < qq[2])/(0.5 * length(tt))
        foldchange = NULL
        if (resp.type == samr.const.twoclass.unpaired.response &
            assay.type == "array") {
            foldchange = foldchange.twoclass(x, y, data$logged2)
        }
        if (resp.type == samr.const.twoclass.paired.response &
            assay.type == "array") {
            foldchange = foldchange.paired(x, y, data$logged2)
        }
        if (resp.type == samr.const.oneclass.response & assay.type ==
            "array") {
        }
        stand.contrasts.95 = NULL
        if (resp.type == samr.const.multiclass.response) {
            stand.contrasts.95 = quantile(stand.contrasts.star,
                c(0.025, 0.975))
        }
        if (resp.type == samr.const.twoclass.unpaired.response &
            assay.type == "seq") {
            foldchange <- foldchange.seq.twoclass.unpaired(x,
                y, depth)
        }
        if (resp.type == samr.const.twoclass.paired.response &
            assay.type == "seq") {
            foldchange <- foldchange.seq.twoclass.paired(x, y,
                depth)
        }
        if (return.x == FALSE) {
            x = NULL
        }
    }
    return(list(n = n, x = x, xresamp = xresamp, y = y, argy = argy,
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
        stand.contrasts.star = stand.contrasts.star, stand.contrasts.95 = stand.contrasts.95,
        depth = depth, call = this.call))
}

#' @title Estimate sequencing depths
#' @param x data matrix. nrow=#gene, ncol=#sample
#' @return depth: estimated sequencing depth. a vector with len sample.
samr.estimate.depth <- function(x) {
	iter <- 5
	cmeans <- colSums(x)/sum(x)
	for (i in 1:iter) {
		n0 <- rowSums(x) %*% t(cmeans)
		prop <- rowSums((x - n0)^2/(n0 + 1e-08))
		qs <- quantile(prop, c(0.25, 0.75))
		keep <- (prop >= qs[1]) & (prop <= qs[2])
		cmeans <- colMeans(x[keep, ])
		cmeans <- cmeans/sum(cmeans)
	}
	depth <- cmeans/mean(cmeans)
	return(depth)
}

#' @title Resampling
#' @param x data matrix. nrow=#gene, ncol=#sample
#' @param d estimated sequencing depth
#' @param nresamp number of resamplings
#' @return xresamp: an rank array with dim #gene*#sample*nresamp
#' @description Corresponds to `samr::resample`
resa <- function(x, d, nresamp = 20) {
	ng <- nrow(x)
	ns <- ncol(x)
	dbar <- exp(mean(log(d)))
	xresamp <- array(0, dim = c(ng, ns, nresamp))
	for (k in 1:nresamp) {
		for (j in 1:ns) {
			xresamp[, j, k] <- rpois(n = ng, lambda = (dbar/d[j]) *
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
rankcols <- function(x) {
	# ranks the elements within each col of the matrix x
	# and returns these ranks in a matrix
	n = nrow(x)
	p = ncol(x)
	mode(n) = "integer"
	mode(p) = "integer"
	mode(x) = "double"
	xr <- apply(x, 2, rank)
	return(xr)
}

#' @title Check format
check.format = function(y, resp.type, censoring.status = NULL) {
	# here i do some format checks for the input data$y
	# note that checks for time course data are done in the
	#   parse function for time course;
	#  we then check the output from the parser in this function
	if (resp.type == samr.const.twoclass.unpaired.response |
		resp.type == samr.const.twoclass.unpaired.timecourse.response) {
		if (sum(y == 1) + sum(y == 2) != length(y)) {
			stop(paste("Error in input response data: response type ",
				resp.type, " specified; values must be 1 or 2"))
		}
	}
	if (resp.type == samr.const.twoclass.paired.response | resp.type ==
		samr.const.twoclass.paired.timecourse.response) {
		if (sum(y) != 0) {
			stop(paste("Error in input response data: response type ",
				resp.type, " specified; values must be -1, 1, -2, 2, etc"))
		}
		if (sum(table(y[y > 0]) != abs(table(y[y < 0])))) {
			stop(paste("Error in input response data:  response type ",
				resp.type, " specified; values must be -1, 1, -2, 2, etc"))
		}
	}
	if (resp.type == samr.const.oneclass.response | resp.type ==
		samr.const.oneclass.timecourse.response) {
		if (sum(y == 1) != length(y)) {
			stop(paste("Error in input response data: response type ",
				resp.type, " specified;  values must all be 1"))
		}
	}
	if (resp.type == samr.const.multiclass.response) {
		tt = table(y)
		nc = length(tt)
		if (sum(y <= nc & y > 0) < length(y)) {
			stop(paste("Error in input response data: response type ",
				resp.type, " specified; values must be 1,2, ... number of classes"))
		}
		for (k in 1:nc) {
			if (sum(y == k) < 2) {
				stop(paste("Error in input response data: response type ",
				  resp.type, " specified; there must be >1 sample per class"))
			}
		}
	}
	if (resp.type == samr.const.quantitative.response) {
		if (!is.numeric(y)) {
			stop(paste("Error in input response data: response type",
				resp.type, " specified; values must be numeric"))
		}
	}
	if (resp.type == samr.const.survival.response) {
		if (is.null(censoring.status)) {
			stop(paste("Error in input response data: response type ",
				resp.type, " specified; error in censoring indicator"))
		}
		if (!is.numeric(y) | sum(y < 0) > 0) {
			stop(paste("Error in input response data:  response type ",
				resp.type, " specified; survival times  must be numeric and nonnegative"))
			if (sum(censoring.status == 0) + sum(censoring.status ==
				1) != length(censoring.status)) {
				stop(paste("Error in input response data: response type ",
				  resp.type, " specified; censoring indicators must be 0 (censored) or 1 (failed)"))
			}
		}
		if (sum(censoring.status == 1) < 1) {
			stop(paste("Error in input response data:   response type ",
				resp.type, " specified; there are no uncensored observations"))
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
	for (i in 1:dim(xresamp)[3]) {
		tt <- tt + rowSums(xresamp[, y == 2, i]) - sum(y == 2) *
			(length(y) + 1)/2
	}
	tt <- tt/dim(xresamp)[3]
	return(list(tt = tt, numer = tt, sd = rep(1, length(tt))))
}

getperms = function(y, nperms) {
	total.perms = factorial(length(y))
	if (total.perms <= nperms) {
		perms = permute(1:length(y))
		all.perms.flag = 1
		nperms.act = total.perms
	}
	if (total.perms > nperms) {
		perms = matrix(NA, nrow = nperms, ncol = length(y))
		for (i in 1:nperms) {
			perms[i, ] = sample(1:length(y), size = length(y))
		}
		all.perms.flag = 0
		nperms.act = nperms
	}
	return(list(perms = perms, all.perms.flag = all.perms.flag,
		nperms.act = nperms.act))
}

#' @title Foldchange of twoclass unpaired sequencing data
foldchange.seq.twoclass.unpaired <- function(x, y, depth)
{
	x.norm <- scale(x, center = F, scale = depth) + 1e-08
	fc <- rowMedians(x.norm[, y == 2])/rowMedians(x.norm[, y ==
		1])
	return(fc)
}
