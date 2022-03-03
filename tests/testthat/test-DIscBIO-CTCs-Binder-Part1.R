notebook_path <- ifelse(interactive(), "../../notebook", "notebook")
FileName <- "CTCdataset"
load(file.path(notebook_path, paste0(FileName,".rda")))
load(file.path(notebook_path, "SC.RData"))
load(file.path(notebook_path, "Ndata.RData"))
load(file.path(notebook_path, "expdata.RData"))
DataSet <- get(FileName)

context("Binder tests, part 1")

test_that("Loading CTC dataset", {
	expect_equal(length(DataSet[, 1]), 13181)
	expect_equal(length(DataSet[1,]), 1462)
})

sc <- DISCBIO(DataSet)
S1 <- summary(colSums(DataSet, na.rm=TRUE))
S2 <- summary(rowMeans(DataSet,na.rm=TRUE))
minexpr <- S2[3]
minnumber <- round(length(DataSet[1, ]) / 10)

test_that("Handling datasets", {
	expect_true(is(sc, "DISCBIO"))
	expect_output(str(S1), " 'summaryDefault' Named num")
	expect_equal(minnumber, 146)
})

sc <- Normalizedata(sc, minexpr=minexpr, minnumber=minnumber, rseed=17000)
sc <- suppressMessages(FinalPreprocessing(sc, GeneFlitering="ExpF"))
sc <- SC
sc@ndata <- Ndata
sc@expdata <- expdata

# Removing the unneeded objects
rm(Ndata)
rm(expdata)
rm(DataSet)
rm(SC)

sc <- Clustexp(sc, 4, quiet=TRUE, rseed=17000)
outlg <- round(length(sc@fdata[,1]) * 0.05)
Outliers <- FindOutliers(sc, 4, outlg=outlg, plot=FALSE, quiet=TRUE)
jcrd <- Jaccard(sc, K=4, plot=FALSE)
sc <- pseudoTimeOrdering(sc, quiet=TRUE)

test_that("Post-processing", {
	expect_equivalent(jcrd, c(.426, .648, .323, .602))
	expect_output(
		object = str(sc, max.level=1),
		expected = 'Formal class \'DISCBIO\' [package "DIscBIO"] with 21 slots'
	)
})
