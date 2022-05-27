data("factRsample")

test_that("Test buildCDS functionality", {
    obj <- buildCDS(factRsample)
    newgenetxs <- obj@txdata
    expect_equal(length(obj@custom), 769)
    expect_equal(nrow(newgenetxs[newgenetxs$cds == "yes",]), 50)
})
