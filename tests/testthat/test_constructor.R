#context("Test factR object")

obj <- CreatefactRObject(gtf, "Mouse")

test_that("Test input arguments", {
    expect_error(CreatefactRObject())
    expect_error(CreatefactRObject(gtf))
    expect_error(CreatefactRObject(gtf, use_own_annotation = "black"))
    expect_error(CreatefactRObject(gtf, use_own_genome = "black"))
    expect_error(CreatefactRObject("wrongpath.gtf", "test"))
    expect_error(CreatefactRObject(gtf, "wronggenome"))
})


test_that("Test correct object structure", {
    expect_equal(as.character(class(obj@custom)), "GRanges")
    expect_equal(as.character(class(obj@reference$ranges)), "GRanges")
    expect_equal(length(obj@custom), 500)
    expect_equal(as.character(class(obj@reference$genome)), "BSgenome")
})
