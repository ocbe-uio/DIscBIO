#' @title Determining differentially expressed genes (DEGs) between two
#'   particular clusters.
#' @description This function defines DEGs between particular clusters generated
#'   by either K-means or model based clustering.
#' @param object \code{DISCBIO} class object.
#' @param Clustering Clustering has to be one of the following:
#'   ["K-means","MB"]. Default is "K-means"
#' @param K A numeric value of the number of clusters.
#' @param fdr A numeric value of the false discovery rate. Default is 0.05.
#' @param name A string vector showing the name to be used to save the resulted
#'   tables.
#' @param First A string vector showing the first target cluster.  Default is
#'   "CL1"
#' @param Second A string vector showing the second target cluster.  Default is
#'   "CL2"
#' @param export A logical vector that allows writing the final gene list in
#'   excel file. Default is TRUE.
#' @param quiet if `TRUE`, suppresses intermediate text output
#' @param plot if `TRUE`, plots are generated
#' @param filename_deg Name of the exported DEG table
#' @param filename_sigdeg Name of the exported sigDEG table
#' @param ... additional parameters to be passed to samr()
#' @importFrom graphics title
#' @importFrom utils write.csv capture.output
#' @importFrom AnnotationDbi keys
#' @return A list containing two tables.
setGeneric(
    "DEGanalysis2clust",
    function(
        object, K, Clustering = "K-means", fdr = 0.05, name = "Name",
        First = "CL1", Second = "CL2",  export = FALSE, quiet = FALSE,
        plot = TRUE, filename_deg = "DEGsTable", filename_sigdeg = "sigDEG",
        ...
    )
    standardGeneric("DEGanalysis2clust")
)

#' @export
#' @rdname DEGanalysis2clust
setMethod(
    "DEGanalysis2clust",
    signature = "DISCBIO",
    definition = function(
        object, K, Clustering, fdr, name, First, Second, export, quiet, plot,
        filename_deg, filename_sigdeg, ...)
    {
        if (!(Clustering %in% c("K-means", "MB"))) {
            stop("Clustering has to be either K-means or MB")
        }
        gene_list <- object@FinalGeneList
        gene_names <- rownames(object@expdata)
        idx_genes <- is.element(gene_names, gene_list)
        gene_names2 <- gene_names[idx_genes]
        dataset <- object@expdata[gene_names2, ]
        Nam <- colnames(dataset)
        if (Clustering == "K-means") {
            Cluster_ID = object@cpart
            if (length(object@cpart) < 1)
                stop("run Clustexp before running DEGanalysis2clust")
        }

        if (Clustering == "MB") {
            Cluster_ID = object@MBclusters$clusterid
            if (length(object@MBclusters$clusterid) < 1)
                stop("run ExprmclustMB before running DEGanalysis2clust")
        }
        num <- c(1:K)
        num1 <- paste("CL", num, sep = "")
        for (n in num) {
            Nam <- ifelse((Cluster_ID == n), num1[n], Nam)
        }
        colnames(dataset) <- Nam
        sg1 <- dataset[, which(colnames(dataset) == First)]
        sg2 <- dataset[, which(colnames(dataset) == Second)]
        sg <- cbind(sg1, sg2)

        sg3 <- factor(
            gsub(paste0("(", First, "|", Second, ").*"), "\\1", colnames(sg)),
            levels = c(paste0(First), paste0(Second))
        )
        sg3 <- sg3[!is.na(sg3)]

        colnames(sg) <- sg3
        len <- c(
            length(sg[, which(colnames(sg) == First)]),
            length(sg[, which(colnames(sg) == Second)])
        )
        y <- c(rep(1:2, len))
        L <- as.matrix(sg)
        gname <- rownames(sg)
        x <- L
        data = list(x = x, y = y, geneid = gname)
        if (quiet) {
            invisible(capture.output({
                samr.obj <- sammy(
                    data,
                    resp.type = "Two class unpaired",
                    assay.type = "seq",
                    testStatistic = "wilcoxon",
                    random.seed = 15,
                    ...
                )
                delta.table <- samr.compute.delta.table(samr.obj)
            }))
        } else {
            samr.obj <- sammy(
                data,
                resp.type = "Two class unpaired",
                assay.type = "seq",
                testStatistic = "wilcoxon",
                random.seed = 15,
                ...
            )
            delta.table <- samr.compute.delta.table(samr.obj)
        }
        DEGsTable <- data.frame()
        DEGsE <- c()
        DEGsS <- c()
        wm <- which.min(delta.table[, 5])
        if (delta.table[wm, 5] <= fdr) {
            w <- which(delta.table[, 5] <= fdr)
            if (is.null (w)) stop("No suitable deltas. Try a lower FDR value.")
            delta <- delta.table[w[1], 1] - 0.001
            if (plot) {
                samr.plot(samr.obj, delta)
                title(paste("DEGs in the", Second, "in", First, "VS", Second))
            }
            siggenes.table <- samr.compute.siggenes.table(
                samr.obj, delta, data, delta.table
            )
            # ------------------------------------------------------------------
            # Reformat siggenes.table as data.frame
            # ------------------------------------------------------------------
            siggenes.table$genes.lo <- reformatSiggenes(siggenes.table$genes.lo)
            siggenes.table$genes.up <- reformatSiggenes(siggenes.table$genes.up)

            FDRl <- as.numeric(siggenes.table$genes.lo[, 8]) / 100
            FDRu <- as.numeric(siggenes.table$genes.up[, 8]) / 100

            siggenes.table$genes.lo[, 8] <- FDRl
            siggenes.table$genes.up[, 8] <- FDRu

            DEGsTable[1, 1] <- paste0(First, " VS ", Second)
            DEGsTable[1, 2] <- Second
            DEGsTable[1, 3] <- length(FDRu)
            DEGsTable[1, 4] <- paste0(
                "Up-regulated-", name, Second, "in", First, "VS", Second,
                ".csv"
            )
            DEGsTable[1, 5] <- length(FDRl)
            DEGsTable[1, 6] <- paste0(
                "Low-regulated-", name, Second, "in", First, "VS", Second,
                ".csv"
            )
            DEGsTable[2, 1] <- paste0(First, " VS ", Second)
            DEGsTable[2, 2] <- First
            DEGsTable[2, 3] <- length(FDRu)
            DEGsTable[2, 4] <- paste0(
                "Low-regulated-", name, First, "in", First, "VS", Second,
                ".csv"
            )
            DEGsTable[2, 5] <- length(FDRl)
            DEGsTable[2, 6] <- paste0(
                "Up-regulated-", name, First, "in", First, "VS", Second,
                ".csv"
            )
            FinalDEGsL <- data.frame()
            if (length(FDRl) > 0) {
                genes <- siggenes.table$genes.lo[, 3]
                if (quiet) {
                    suppressMessages(
                        geneList <- AnnotationDbi::select(
                            org.Hs.eg.db,
                            keys = keys(org.Hs.eg.db),
                            columns = c("SYMBOL", "ENSEMBL")
                        )
                    )
                    GL <- c(1, "MTRNR2", "ENSG00000210082")
                    GL1 <- c(1, "MTRNR1", "ENSG00000211459")
                    geneList <- rbind(geneList, GL, GL1)
                } else {
                    geneList <- AnnotationDbi::select(
                        org.Hs.eg.db,
                        keys = keys(org.Hs.eg.db),
                        columns = c("SYMBOL", "ENSEMBL")
                    )
                    GL <- c(1, "MTRNR2", "ENSG00000210082")
                    GL1 <- c(1, "MTRNR1", "ENSG00000211459")
                    geneList <- rbind(geneList, GL, GL1)
                }
                FinalDEGsL <- cbind(genes, siggenes.table$genes.lo)
                gene_list <- geneList[, 3]
                idx_genes <- is.element(gene_list, genes)
                genes2 <- geneList[idx_genes, ]
                if (!is.null(FinalDEGsL)) {
                    FinalDEGsL <- merge(
                        FinalDEGsL,
                        genes2,
                        by.x = "genes",
                        by.y = "ENSEMBL",
                        all.x = TRUE
                    )
                    FinalDEGsL[, 3] <- FinalDEGsL[, 11]
                    FinalDEGsL <- FinalDEGsL[, c(-1, -10, -11)]
                    FinalDEGsL <- FinalDEGsL[order(FinalDEGsL[, 8]), ]
                    FinalDEGsL[is.na(FinalDEGsL[, 2]), c(2, 3)] <-
                        FinalDEGsL[is.na(FinalDEGsL[, 2]), 3]
                }
                if (export) {
                    message("The results of DEGs are saved in your directory")
                    message(
                        "Low-regulated genes in the ", Second, " in ",
                        First, " VS ", Second, "\n"
                    )
                    write.csv(
                        FinalDEGsL,
                        file = paste0(
                            "Low-regulated-", name, Second, "in", First,
                            "VS", Second, ".csv"
                        )
                    )
                    write.csv(
                        FinalDEGsL,
                        file = paste0(
                            "Up-regulated-", name, First, "in", First, "VS",
                            Second, ".csv"
                        )
                    )
                }
                DEGsS <- c(DEGsS, FinalDEGsL[, 2])
                DEGsE <- c(DEGsE, as.character(FinalDEGsL[, 3]))
            }
            FinalDEGsU <- data.frame()
            if (length(FDRu) > 0) {
                genes <- siggenes.table$genes.up[, 3]
                if (quiet) {
                    suppressMessages(
                        geneList <- AnnotationDbi::select(
                            org.Hs.eg.db,
                            keys = keys(org.Hs.eg.db),
                            columns = c("SYMBOL", "ENSEMBL")
                        )
                    )
                    GL <- c(1, "MTRNR2", "ENSG00000210082")
                    geneList <- rbind(geneList, GL)
                } else {
                    geneList <- AnnotationDbi::select(
                        org.Hs.eg.db,
                        keys = keys(org.Hs.eg.db),
                        columns = c("SYMBOL", "ENSEMBL")
                    )
                    GL <- c(1, "MTRNR2", "ENSG00000210082")
                    geneList <- rbind(geneList, GL)
                }
                FinalDEGsU <- cbind(genes, siggenes.table$genes.up)
                gene_list <- geneList[, 3]
                idx_genes <- is.element(gene_list, genes)
                genes2 <- geneList[idx_genes, ]
                if (!is.null(FinalDEGsU)) {
                    FinalDEGsU <- merge(
                        FinalDEGsU,
                        genes2,
                        by.x = "genes",
                        by.y = "ENSEMBL",
                        all.x = TRUE
                    )
                    FinalDEGsU[, 3] <- FinalDEGsU[, 11]
                    FinalDEGsU <- FinalDEGsU[, c(-1, -10, -11)]
                    FinalDEGsU <- FinalDEGsU[order(FinalDEGsU[, 8]), ]
                    FinalDEGsU[is.na(FinalDEGsU[, 2]), c(2, 3)] <-
                        FinalDEGsU[is.na(FinalDEGsU[, 2]), 3]
                }
                if (export) {
                    message("The results of DEGs are saved in your directory")
                    message(
                        "Up-regulated genes in the ", Second, " in ", First,
                        " VS ", Second, "\n"
                    )
                    write.csv(
                        FinalDEGsU,
                        file = paste0(
                            "Up-regulated-", name, Second, "in", First, "VS",
                            Second, ".csv"
                        )
                    )
                    write.csv(
                        FinalDEGsU,
                        file = paste0(
                            "Low-regulated-", name, First, "in", First, "VS",
                            Second, ".csv"
                        )
                    )
                }
                DEGsS <- c(DEGsS, FinalDEGsU[, 2])
                DEGsE <- c(DEGsE, as.character(FinalDEGsU[, 3]))
            }
        } else{
            DEGsTable[1, 1] <- paste0(First, " VS ", Second)
            DEGsTable[1, 2] <- Second
            DEGsTable[1, 3] <- NA
            DEGsTable[1, 4] <- NA
            DEGsTable[1, 5] <- NA
            DEGsTable[1, 6] <- NA

            DEGsTable[2, 1] <- paste0(First, " VS ", Second)
            DEGsTable[2, 2] <- First
            DEGsTable[2, 3] <- NA
            DEGsTable[2, 4] <- NA
            DEGsTable[2, 5] <- NA
            DEGsTable[2, 6] <- NA
        }
        colnames(DEGsTable) <- c(
            "Comparisons",
            "Target cluster",
            "Gene number",
            "File name",
            "Gene number",
            "File name"
        )
        if (!quiet) print(DEGsTable)
        sigDEG <- cbind(DEGsE, DEGsS)
        if (export) {
            write.csv(DEGsTable, file = paste0(filename_deg, ".csv"))
            write.csv(sigDEG, file = paste0(filename_sigdeg, ".csv"))
        }
        return(
            list(
                sigDEG = sigDEG,
                DEGsTable = DEGsTable,
                FinalDEGsL = FinalDEGsL,
                FinalDEGsU = FinalDEGsU

            )
        )
    }
)