data("factRsample")

test_that("Test .runAS functionality", {
    expect_error(predictNMD(factRsample))
    obj <- buildCDS(factRsample)
    obj <- predictNMD(obj)

    expect_equal(txData(obj)$is_NMD[1:5], c(T,F,T,T,T))
    expect_equal(nrow(obj@nmd[obj@nmd$is_NMD,]), 6)
})
