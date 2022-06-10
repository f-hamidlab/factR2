data('factRsample')

counts.path <- system.file("extdata/counts.txt", package = "factR2")
meta.path <- system.file("extdata/meta.txt", package = "factR2")

cts <- read.delim(counts.path)
cts <- as.matrix(cts)
meta <- read.delim(meta.path)

obj <- addTxCounts(factRsample, cts, meta, ~samples)
test_that("Test addCountsData function", {
    expect_error(addTxCounts(factRsample))
    expect_error(addTxCounts(factRsample, cts[1:5,], meta, ~samples))
    countSample.nonames <- cts
    rownames(countSample.nonames) <- NULL
    expect_error(addTxCounts(factRsample, countSample.nonames, meta, ~samples))
    expect_equal(obj@sets$transcript@counts, cts)
})


# add design














obj2 <- addSampleData(factRsample, testSamples)
test_that("Test addSampleData", {

    expect_equal(obj2@colData, testSamples)
    expect_error(addSampleData(obj, testSamples[1:3,]))
    samples.nonames <- testSamples
    rownames(samples.nonames) <- NULL
    expect_error(addSampleData(obj, samples.nonames[1:3,]))
    expect_error(addSampleData(obj, samples.nonames[,c("group", "group") ]))
    dat <- addSampleData(obj, samples.nonames)
    expect_equal(dat@colData, testSamples)

})

obj3 <- addTxCounts(factRsample, countsSample, sampleData = testSamples)
test_that("Test addCountsData with samples added", {
    expect_error(addTxCounts(obj2, countsSample[,1:3]))
    obj4 <- addSampleData(obj, testSamples)
    expect_equal(obj3, obj4)
})

obj5 <- addTxCounts(factRsample, countsSample, testSamples, ~samples)
test_that("Test addDesign function", {
    expect_error(addDesign(factRsample, ~samples))
    expect_equal(addDesign(obj3, ~samples), obj5)
})
