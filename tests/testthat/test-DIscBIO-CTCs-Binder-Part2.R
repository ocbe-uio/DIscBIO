# ======================================================== #
# Loading and rearranging files                            #
# ======================================================== #

notebook_path <- ifelse(interactive(), "../../notebook", "notebook")
load(file.path(notebook_path, "SC.RData"))
load(file.path(notebook_path, "Ndata.RData"))
load(file.path(notebook_path, "expdata.RData"))
load(file.path(notebook_path, "DATAforDT.RData"))

sc <- SC
sc@ndata <- Ndata
sc@expdata <- expdata

rm(Ndata)
rm(expdata)
rm(SC)

# ======================================================== #
# differential expression analysis                         #
# ======================================================== #

context("Binder tests, part 2: Differential expression analysis")

cdiff <- DEGanalysis2clust(sc, 4, quiet=TRUE, plot=FALSE)

test_that("DEG analysis between 2 clusters", {
	expect_equal(
		object = head(cdiff[[1]])[, 1],
		expected = c(
			'ENSG00000008988', 'ENSG00000010278', 'ENSG00000034510',
			'ENSG00000071082', 'ENSG00000071127', 'ENSG00000075624'
		)
	)
	expect_equal(
		object = head(cdiff[[1]])[, 2],
		expected = c('RPS20', 'CD9', 'TMSB10', 'RPL31', 'WDR1', 'ACTB')
	)
	expect_equivalent(
		object = as.character(head(cdiff[[2]])[1, ]),
		expected = c(
			'CL1 VS CL2', 'CL2', '106', 'Up-regulated-NameCL2inCL1VSCL2.csv',
			'82', 'Low-regulated-NameCL2inCL1VSCL2.csv'
		)
	)
	expect_equal(
		object = as.character(head(cdiff[[2]])[2, ]),
		expected = c(
			'CL1 VS CL2', 'CL1', '106', 'Low-regulated-NameCL1inCL1VSCL2.csv',
			'82', 'Up-regulated-NameCL1inCL1VSCL2.csv'
		)
	)
})

cdiffBinomial <- ClustDiffGenes(sc, 4, quiet=TRUE)

test_that("Cluster differences", {
		expect_equal(
		object = head(cdiffBinomial[[1]])[, 1],
		expected = c(
			'ENSG00000001630', 'ENSG00000002586', 'ENSG00000003402',
			'ENSG00000003436', 'ENSG00000003756', 'ENSG00000004059'
		)
	)
	expect_equal(
		object = head(cdiffBinomial[[1]])[, 2],
		expected = c('CYP51A1', 'CD99', 'CFLAR', 'TFPI', 'RBM5', 'ARF5')
	)
	expect_equivalent(
		object = as.character(head(cdiffBinomial[[2]])[1, ]),
		expected = c(
			'Cluster 1', 'Remaining Clusters', '1052', 'Up-DEG-cluster1.csv',
			'678', 'Down-DEG-cluster1.csv'
		)
	)
	expect_equal(
		object = as.character(head(cdiffBinomial[[2]])[2, ]),
		expected = c(
			'Cluster 2', 'Remaining Clusters', '0', 'Up-DEG-cluster2.csv', '1',
			'Down-DEG-cluster2.csv'
		)
	)
})


# ======================================================== #
# Decision trees                                           #
# ======================================================== #

context("Binder tests, part 2: Decision trees")

j48dt <- J48DT(DATAforDT, plot=FALSE, quiet=TRUE)
rpartDT <- RpartDT(DATAforDT, plot=FALSE, quiet=TRUE)

test_that("J48 trees", {
	expect_true(is(summary(j48dt), "Weka_classifier_evaluation"))
	expect_output(str(rpartDT), "List of 14")
})

# ======================================================== #
# Networking                                               #
# ======================================================== #

context("Binder tests, part 2: Networking")

data <- cdiffBinomial[[1]] [1:200,2] # only the firat 200 genes
ppi <- suppressMessages(PPI(data))
networking <- suppressMessages(NetAnalysis(ppi))

test_that("Networks", {
	expect_equal(dim(ppi), c(902, 13))
	expect_equal(dim(networking), c(174, 3))
})
