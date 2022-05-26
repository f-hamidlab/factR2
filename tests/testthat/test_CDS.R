data("factRsample")

test_that("Test BuildCDS functionality", {
    obj <- BuildCDS(factRsample)
    newgenetxs <- obj@custom$genetxs
    expect_equal(length(obj@custom$ranges), 769)
    expect_equal(nrow(newgenetxs[newgenetxs$cds == "yes",]), 50)
})
