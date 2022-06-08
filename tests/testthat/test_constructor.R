#context("Test factR object")

obj <- createfactRObject(gtf, "vM25")

test_that("Test input arguments", {
    expect_error(createfactRObject())
    expect_error(createfactRObject(gtf))
    expect_error(createfactRObject(gtf, use_own_annotation = "black"))
    expect_error(createfactRObject(gtf, use_own_genome = "black"))
    expect_error(createfactRObject("wrongpath.gtf", "test"))
    expect_error(createfactRObject(gtf, "wronggenome"))
})


test_that("Test correct object structure", {
    expect_equal(as.character(class(obj@transcriptome)), "GRanges")
    expect_equal(as.character(class(obj@reference$ranges)), "GRanges")
    expect_equal(length(obj@transcriptome), 8117)
    expect_equal(as.character(class(obj@reference$genome)), "BSgenome")
})


test_that("Test .prepTranscriptome", {
    expect_false(identical(obj@transcriptome, factRsample@transcriptome))
    expect_equal(ncol(mcols(obj@transcriptome)), 12)
    expect_equal(length(obj@sets), 3)
    
    expect_equal(nrow(obj@sets$transcript@rowData), 154)
    expect_equal(nrow(obj@sets$gene@rowData), 41)
    expect_equal(nrow(obj@sets$AS@rowData), 133)
    
    featdata <- featureData(obj[["transcript"]])
    expect_equal(nrow(featdata[featdata$novel == "yes",]), 46)
})
