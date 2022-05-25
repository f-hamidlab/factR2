test_that("Test object inputs", {
    out <- .smartimport(genomes, "")
    expect_equal(out, genomes)
})

test_that("Test local inputs", {
    out <- .smartimport(gtf, ".gtf")
    expect_equal(length(out), 500)
    expect_equal(ncol(mcols(out)), 9)
})

test_that("Test FTP inputs", {
    out <- .smartimport(genomes[1,]$annotation, ".gtf")
    expect_equal(length(out), 3542517)
    expect_equal(ncol(mcols(out)), 21)

})

test_that("Test AH inputs", {
    out <- .smartimport("AH49010", ".gtf")
    expect_equal(length(out), 3483880)
    expect_equal(ncol(mcols(out)), 18)
})

test_that("Test wrong input", {
    expect_error(.smartimport("wrong", ".gtf"))
})
