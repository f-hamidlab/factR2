#context("Test factR object")

gtf <- system.file("extdata", "sc_merged_sample.gtf.gz", package = "factR")
out <- CreatefactRObject(gtf, gtf, Mmusculus)

test_that("Test input arguments", {
    expect_error(CreatefactRObject())
    expect_error(CreatefactRObject("wrongpath.gtf", "test", "test"))
    expect_error(CreatefactRObject(gtf, "wronggenome", "test"))
    expect_error(CreatefactRObject(gtf, gtf, "test"))
})


test_that("Test correct object structure", {
    expect_equal(as.character(class(out@ranges$custom)), "GRanges")
    expect_equal(as.character(class(out@ranges$ref)), "GRanges")
    expect_equal(length(out@ranges$custom), 500)
    expect_equal(as.character(class(out@genome$seq)), "BSgenome")

})
