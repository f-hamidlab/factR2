data("factRsample")




obj <- buildCDS(factRsample)
test_that("Test buildCDS functionality", {
    expect_equal(length(obj@transcriptome), 2069)
    featdata <- obj[["transcript"]]
    expect_equal(nrow(featdata[featdata$cds == "yes",]), 83)
})


obj <- predictNMD(obj)
test_that("Test predictNMD functionality", {
    expect_error(predictNMD(factRsample))
    featdata <- obj[["transcript"]]
    expect_equal(nrow(featdata[featdata$nmd == "yes",]), 17)
})

obj <- getAAsequence(obj)
test_that("Test getAAsequence functionality", {
    expect_error(getAAsequence(factRsample))
    expect_equal(length(obj@domains$sequence), 83)
    expect_equal(as.character(class(obj@domains$sequence)), "AAStringSet")
})


ref <- .getbestref(obj@transcriptome, obj[["transcript"]])
test_that("Test .getbestref function", {
    expect_equal(length(ref), 439)
})

obj <- testASNMDevents(obj)
test_that("Test testASNMD functionality", {
    expect_error(testASNMDevents(factRsample))
    featdata <- obj[["AS"]]
    expect_equal(nrow(featdata[featdata$ASNMDtype == "Repressing",]), 16)
})


test_that("Test runfactR functionality", {
    expect_equal(obj, suppressWarnings(runfactR(factRsample)))
})
