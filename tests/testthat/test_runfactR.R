data("factRsample")

obj <- prepTranscriptome(factRsample)
test_that("Test prepTranscriptome", {
    expect_false(identical(obj@transcriptome, factRsample@transcriptome))
    expect_equal(ncol(mcols(obj@transcriptome)), 12)
    expect_equal(length(obj@sets), 3)

    expect_equal(nrow(obj@sets$transcript@rowData), 154)
    expect_equal(nrow(obj@sets$gene@rowData), 41)
    expect_equal(nrow(obj@sets$AS@rowData), 133)

    featdata <- featureData(obj[["transcript"]])
    expect_equal(nrow(featdata[featdata$novel == "yes",]), 46)
})


obj <- buildCDS(obj)
test_that("Test buildCDS functionality", {
    expect_error(buildCDS(factRsample))
    expect_equal(length(obj@transcriptome), 2069)
    featdata <- featureData(obj[["transcript"]])
    expect_equal(nrow(featdata[featdata$cds == "yes",]), 83)
})


obj <- predictNMD(obj)
test_that("Test predictNMD functionality", {
    expect_error(predictNMD(factRsample))
    featdata <- featureData(obj[["transcript"]])
    expect_equal(nrow(featdata[featdata$nmd == "yes",]), 17)
})

obj <- getAAsequence(obj)
test_that("Test getAAsequence functionality", {
    expect_error(getAAsequence(factRsample))
    expect_equal(length(obj@domains$sequence), 83)
    expect_equal(as.character(class(obj@domains$sequence)), "AAStringSet")
})


ref <- .getbestref(obj)
test_that("Test .getbestref function", {
    expect_equal(length(ref), 192)
})

obj <- testASNMDevents(obj)
test_that("Test testASNMD functionality", {
    expect_error(testASNMDevents(factRsample))
    featdata <- featureData(obj[["AS"]])
    expect_equal(nrow(featdata[featdata$ASNMDtype == "Repressing",]), 16)
})


test_that("Test runfactR functionality", {
    expect_equal(obj, runfactR(factRsample))
})
