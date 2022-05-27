data("factRsample")

test_that("Test BuildCDS functionality", {
    obj <- BuildCDS(factRsample)
    newgenetxs <- obj@txdata
    expect_equal(length(obj@custom), 769)
    expect_equal(nrow(newgenetxs[newgenetxs$cds == "yes",]), 50)
})
