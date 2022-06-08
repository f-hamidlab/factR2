data('factRsample')
data('countsSample')

obj <- addCountData(factRsample, countsSample)
test_that("Test addCountsData function", {
    expect_error(addCountData(factRsample))
    expect_error(addCountData(factRsample, countsSample[1:5,]))
    countSample.nonames <- countsSample
    rownames(countSample.nonames) <- NULL
    expect_error(addCountData(factRsample, countSample.nonames))
    expect_equal(obj@sets$transcript@counts, countsSample)
})

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

obj3 <- addCountData(factRsample, countsSample, sampleData = testSamples)
test_that("Test addCountsData with samples added", {
    expect_error(addCountData(obj2, countsSample[,1:3]))
    obj4 <- addSampleData(obj, testSamples)
    expect_equal(obj3, obj4)
})

obj5 <- addCountData(factRsample, countsSample, testSamples, ~samples)
test_that("Test addDesign function", {
    expect_error(addDesign(factRsample, ~samples))
    expect_equal(addDesign(obj3, ~samples), obj5)
})
