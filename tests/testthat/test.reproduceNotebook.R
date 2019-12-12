# This test makes sure the package works with respect to the interactive notebook

context("Reproducing Jupyter demo")

test_that("Loading datasets generate the expected output", {
    expect_equal(length(valuesG1ms[, 1]), 59838)
    expect_equal(length(valuesG1ms[1, ]), 94)
    expect_equal(length(MLSrawWithoutERCC[, 1]), 59746)
    expect_equal(length(MLSrawWithoutERCC[1, ]), 94)
})

test_that("Data pre-processing results are reproduced", {
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
        plot = FALSE, export = FALSE
    )       # Val=F  will plot all the ERCC spike-ins

    # sc <- PSCANseq(MLSrawWithoutERCC)
    # sc<- Normalizedata(sc,mintotal=1000)
    # gene_list<-noiseF
    # gene_names<-rownames(sc@ndata)
    # idx_genes <- is.element(gene_names, gene_list)
    # gene_names2 <- gene_names[idx_genes]
    # filteredDataset<- sc@ndata[gene_names2,]
    # cat("The gene filtering method= Noise filtering","\n","\n") 
    # cat("The Filtered Normalized dataset is called: filteredDataset","\n","\n") 
    # cat("The filtered Normalized dataset contains:","\n","Genes:",length(filteredDataset[,1]),"\n","cells:",length(filteredDataset[1,]),"\n","\n")
    # save(filteredDataset,file="filteredDataset.Rdata")
    # sc@fdata<-filteredDataset

    # ############# Generating a filtered dataset with raw data for DEG analysis
    # gene_list<-noiseF
    # gene_names<-rownames(MLSrawWithoutERCC)
    # idx_genes <- is.element(gene_names, gene_list)
    # gene_names2 <- gene_names[idx_genes]
    # LipoNoisFilteredRawDataset<- MLSrawWithoutERCC[gene_names2,]
    # cat("The Nois Filtered Raw dataset contains:","\n","Genes:",length(LipoNoisFilteredRawDataset[,1]),"\n","cells:",length(LipoNoisFilteredRawDataset[1,]))
    # # save(LipoNoisFilteredRawDataset,file="LipoNoisFilteredRawDataset.Rdata")

    expect_equal(class(noiseF), "character")
    expect_equal(length(noiseF), 5684)
    expect_equal(noiseF[1234], "ENSG00000105852")
})