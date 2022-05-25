ah <- suppressMessages(AnnotationHub::AnnotationHub())

test_that("Test all links", {
    annotationURLs <- stringr::str_subset(genomes$annotation, "^http")
    annotationAH <- stringr::str_subset(genomes$annotation, "^AH")
    expect_true(all(RCurl::url.exists(annotationURLs)))
    expect_true(all(annotationAH %in% names(ah)))

    genomeURLs <- stringr::str_subset(genomes$genome.sec, "^http")
    genomeAH <- stringr::str_subset(genomes$genome.sec, "^AH")
    expect_true(all(RCurl::url.exists(genomeURLs)))
    expect_true(all(genomeAH %in% names(ah)))
})
