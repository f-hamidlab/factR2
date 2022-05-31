object <- buildCDS(factRsample)
object <- getAAsequence(object)
out <- .runDomainSearch(object@domains$sequence[1,])
out2 <- .runDomainSearch(object@domains$sequence[1,], db = "pfam")

test_that("Test getAAsequence output", {
    expect_equal(nrow(object@domains$sequence), 50)
    expect_true(stringr::str_detect(object@domains$sequence[1,]$x,
                                    "MASFGWKRRIGEKVSKATSQQFEAEAADEKDAA"))
})


test_that("Test .runDomainSearch output", {
   expect_equal(out$type, c("DOMAIN", "CHAIN"))
   expect_equal(out$description[1], "TPR-like")
   expect_equal(out2$description[1], "Tetratricopeptide repeat")
})
