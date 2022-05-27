data("factRsample")

test_that("Test head/tail", {
    expect_equal(factR2::head(factRsample), utils::head(factRsample@custom))
    expect_equal(factR2::tail(factRsample), utils::tail(factRsample@custom))
})

test_that("Test View", {
    expect_error(View(factRsample, "WRONGENE"))
    expect_warning(View(factRsample, "WRONGENE", "Selenop"))
})

test_that("Test Plot", {
    expect_error(Plot(factRsample, "WRONGENE"))
    expect_warning(Plot(factRsample, "WRONGENE", "Selenop"))
    plot <- Plot(factRsample, "Selenop")
    expect_equal(c("patchwork", "gg","ggplot"), as.character(class(plot)))
})
