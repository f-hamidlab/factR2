data("factRsample")

test_that("Test .runAS functionality", {
    expect_error(predictNMD(factRsample))
    obj <- BuildCDS(factRsample)
    obj <- predictNMD(obj)

    expect_equal(as.character(class(obj@nmd))[1], "tbl_df")
    expect_equal(nrow(obj@nmd), 50)
    expect_equal(nrow(obj@nmd[obj@nmd$is_NMD,]), 6)
})
