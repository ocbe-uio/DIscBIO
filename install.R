if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(
    c(
        "pheatmap**currentlyMASS", "cluster", "mclust", "flexmix",
        "lattice", "fpc", "amap", "RColorBrewer", "locfit", "TSCAN",
        "genefilter", "statmod", "ggplot2", "gplots", "DESeq2",
        "matrixStats", "robustbase", "philentropy", "igraph", "boot",
        "biomaRt", "tidyr", "calibrate", "partykit", "RWeka", "rpart",
        "rpart.plot", "imager", "png", "NetIndices", "httr", "jsonlite",
        "tidyverse", "samr", "tidyverse","org.Hs.eg.db", "AnnotationDbi",
        "enrichR"
    )
)

library("devtools")
devtools::install(".", build_vignettes=TRUE)
