obj <- runfactR(factRsample)
ref <- .getbestref(obj)
obj <- testASNMDevents(obj)

test_that("Test .getbestref function", {
    expect_equal(length(ref), 192)
})

test_that("Test testASNMDevents function", {
    expect_equal(obj@ASplicings$ASNMDtype[1:5],
                 c("Repressing",rep("Stimulating", 4)))
})
