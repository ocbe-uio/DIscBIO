# This test makes sure the package works with respect to the interactive notebook

context("Reproducing Jupyter demo")

test_that("Loading datasets generate the expected output", {
    expect_equal(dim(valuesG1ms), c(59838, 94))
    expect_equal(dim(MLSrawWithoutERCC), c(59746, 94))
})

test_that("Data pre-processing results are reproduced", {
    # Code from notebook =======================================================
    percentile <- 0.8
    CV <- 0.3 
    Object <- valuesG1ms
    gene_names <- GeneNames

    # Remove .x from gene ID's, since they can't be handled later on
    gene_names <- as.list(sub("*\\..*", "", unlist(gene_names)))

    gene_names2 <- head(gene_names, -92)

    geneCol  <- "yellow"
    FgeneCol <- "black"
    erccCol  <- "blue"
    noiseF   <- NoiseFiltering(
        Object, percentile, CV, gene_names2, geneCol, FgeneCol, erccCol, Val=T,
        plot = FALSE, export = FALSE, quiet = TRUE
    )       # Val=F  will plot all the ERCC spike-ins

    sc <- PSCANseq(MLSrawWithoutERCC)
    sc <- Normalizedata(sc, mintotal = 1000)
    gene_list <- noiseF
    gene_names <- rownames(sc@ndata)
    idx_genes <- is.element(gene_names, gene_list)
    gene_names2 <- gene_names[idx_genes]
    filteredDataset <- sc@ndata[gene_names2, ]

    # ASK: should these cats be part of some function output?
    # cat("The gene filtering method= Noise filtering","\n","\n") 
    # cat("The Filtered Normalized dataset is called: filteredDataset","\n","\n") 
    # cat("The filtered Normalized dataset contains:","\n","Genes:",length(filteredDataset[,1]),"\n","cells:",length(filteredDataset[1,]),"\n","\n")
    # save(filteredDataset,file="filteredDataset.Rdata")
    sc@fdata <- filteredDataset

    # ############# Generating a filtered dataset with raw data for DEG analysis
    gene_list <- noiseF
    gene_names <- rownames(MLSrawWithoutERCC)
    idx_genes <- is.element(gene_names, gene_list)
    gene_names2 <- gene_names[idx_genes]
    LipoNoisFilteredRawDataset <- MLSrawWithoutERCC[gene_names2, ]
    
    # ASK: should these cats be part of some function output?  
    # cat("The Noise Filtered Raw dataset contains:","\n","Genes:",length(LipoNoisFilteredRawDataset[,1]),"\n","cells:",length(LipoNoisFilteredRawDataset[1,]))
    # # save(LipoNoisFilteredRawDataset,file="LipoNoisFilteredRawDataset.Rdata")

    # Unit tests ===============================================================
    expect_equal(class(noiseF), "character")
    expect_equal(length(noiseF), 5684)
    expect_equal(noiseF[1234], "ENSG00000105852")
    expect_equal(dim(filteredDataset), c(5684, 94))
    expect_equal(dim(LipoNoisFilteredRawDataset), c(5684, 94))
})