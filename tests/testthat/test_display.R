data("factRsample")

test_that("Test head/tail", {
    expect_equal(factR2::head(factRsample), utils::head(factRsample@custom))
    expect_equal(factR2::tail(factRsample), utils::tail(factRsample@custom))
})

test_that("Test View", {
    expect_error(view(factRsample, "WRONGENE"))
    expect_warning(view(factRsample, "WRONGENE", "Selenop"))
})

test_that("Test plotTranscripts", {
    expect_error(plotTranscripts(factRsample, "WRONGENE"))
    expect_warning(plotTranscripts(factRsample, "WRONGENE", "Selenop"))
    plot <- plotTranscripts(factRsample, "Selenop")
    expect_equal(c("patchwork", "gg","ggplot"), as.character(class(plot)))
})

obj <- runfactR(factRsample)
test_that("Test plotDomains", {
    expect_error(plotDomains(obj, "WRONGENE"))
    expect_warning(plotDomains(obj, "WRONGENE", "Selenop"))
    plot <- suppressWarnings(plotDomains(obj, "Ttc33"))
    expect_equal(c("patchwork", "gg","ggplot"), as.character(class(plot)))
})

