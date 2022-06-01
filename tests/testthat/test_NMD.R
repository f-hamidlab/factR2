data("factRsample")

test_that("Test .runAS functionality", {
    expect_error(predictNMD(factRsample))
    obj <- buildCDS(factRsample)
    obj <- predictNMD(obj)
    genetxs <- txData(obj)

    expect_equal(as.logical(genetxs$is_NMD[1:5]), c(T,F,F,T,F))
    expect_equal(nrow(genetxs[genetxs$nmd == "yes",]), 6)
})
