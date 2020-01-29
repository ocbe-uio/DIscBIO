# This test makes sure the package works with respect to the interactive notebook
context("Reproducing notebook: loading and pre-processing")

# Loading datasets =============================================================
test_that("Loading datasets generate the expected output", {
    expect_equal(dim(valuesG1ms), c(59838, 94))
})