data(genomes)

test_that("Test all links", {
    expect_true(all(RCurl::url.exists(genomes$annotation)))
    expect_true(all(RCurl::url.exists(genomes$genome)))
})
