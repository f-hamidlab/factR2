# setup GRanges
gr1 <- GenomicRanges::GRanges(
    seqnames = "1", strand = rep("+", 2),
    ranges = IRanges::IRanges(
        start = c(1, 501),
        end = c(100, 600)
    )
)

gr2 <- GenomicRanges::GRanges(
    seqnames = "1", strand = rep("+", 3),
    ranges = IRanges::IRanges(
        start = c(1, 201, 501),
        end = c(100, 300, 600)
    )
)

gr3 <- GenomicRanges::GRanges(
    seqnames = "1", strand = rep("+", 2),
    ranges = IRanges::IRanges(
        start = c(1,  401),
        end = c(100, 600)
    )
)

gr4 <- GenomicRanges::GRanges(
    seqnames = "1", strand = rep("+", 1),
    ranges = IRanges::IRanges(
        start = c(1),
        end = c(600)
    )
)

exon1 <- GenomicRanges::GRanges(
    seqnames = "1", strand = rep("+", 1),
    ranges = IRanges::IRanges(
        start = c(201),
        end = c(300)
    )
)

exon2 <- GenomicRanges::GRanges(
    seqnames = "1", strand = rep("+", 1),
    ranges = IRanges::IRanges(
        start = c(401),
        end = c(500)
    )
)

exon3 <- GenomicRanges::GRanges(
    seqnames = "1", strand = rep("+", 1),
    ranges = IRanges::IRanges(
        start = c(101),
        end = c(500)
    )
)

exonset1 <- c(exon1,exon2,exon3)
exonset2 <- exonset1 %>% as.data.frame() %>%
    dplyr::mutate(seqnames = c("1", "2", "1")) %>%
    dplyr::mutate(strand = c("-", "+", "-")) %>%
    GenomicRanges::makeGRangesFromDataFrame()
exonset3 <- exonset1 %>% as.data.frame() %>%
    dplyr::mutate(seqnames = c("1", "2", "1")) %>%
    dplyr::mutate(strand = c("+", "+", "-")) %>%
    GenomicRanges::makeGRangesFromDataFrame()

grl1 <- GenomicRanges::GRangesList(gr1,gr1,gr1)
names(grl1) <- c("test1", "test2", "tes3")

grl2 <- GenomicRanges::GRangesList(gr2,gr3,gr4)
names(grl2) <- c("test4", "test5", "test6")


test_that("Test .addExonstoTx function", {
    out <- .addExonstoTx(grl1, exonset1)
    expect_equal(ranges(out[[1]]), ranges(gr2))
    expect_equal(ranges(out[[2]]), ranges(gr3))
    expect_equal(ranges(out[[3]]), ranges(gr4))

    # test lengths of x and y
    expect_error(.addExonstoTx(grl2, c(exon1,exon2)))

    # test function when all exons are on diff chr
    expect_error(.addExonstoTx(grl1, exonset2))

    # test function when some exons are on diff chromosomes/strand
    # #### default parameters
    out <- suppressWarnings(.addExonstoTx(grl1, exonset3))
    expect_equal(length(out), 3)
    expect_equal(as.data.frame(out)$start, c(1,201,501,1,501, 1,501))
    expect_equal(as.data.frame(out)$modified, c(T,T,T,F,F,F,F))

    # #### when drop.unmodified = TRUE
    expect_warning(.addExonstoTx(grl1, exonset3, drop.unmodified = TRUE))
    out <- suppressWarnings(.addExonstoTx(grl1, exonset3, drop.unmodified = TRUE))
    expect_equal(length(out), 1)
    expect_equal(as.data.frame(out)$start, as.data.frame(gr2)$start)

    # when external exons are added but allow.external.exons = FALSE
    exon4 <- GenomicRanges::GRanges(
        seqnames = "1", strand = rep("+", 1),
        ranges = IRanges::IRanges(
            start = c(801),
            end = c(1000)
        )
    )
    out <- suppressWarnings(.addExonstoTx(grl1, c(exon4,exon4,exon4),
                                         allow.external.exons = FALSE))
    expect_equal(length(out), 3)
    expect_equal(as.data.frame(out)$start, as.data.frame(grl1)$start)

    # when external exons are added but allow.external.exons = TRUE
    out <- suppressWarnings(.addExonstoTx(grl1, c(exon4,exon4,exon4),
                                         allow.external.exons = TRUE))
    expect_warning(.addExonstoTx(grl1, c(exon4,exon4,exon4),
                                allow.external.exons = TRUE))
    expect_equal(length(out), 3)
    expect_equal(as.data.frame(out)$start, rep(c(1,501,801),3))

})

test_that("Test .removeExonsfromTx function", {
    # test general functionality
    out <- .removeExonsfromTx(grl2, exonset1)
    expect_equal(ranges(out[[1]]), ranges(gr1))
    expect_equal(ranges(out[[2]]), ranges(gr1))
    expect_equal(ranges(out[[3]]), ranges(gr1))

    # test lengths of x and y
    expect_error(.removeExonsfromTx(grl2, c(exon1,exon2)))

    # test function when all exons are on diff chromosomes/strand
    expect_error(.removeExonsfromTx(grl2, exonset2))

    # test function when some exons are on diff chromosomes/strand
    #### default parameters
    out <- suppressWarnings(.removeExonsfromTx(grl2, exonset3))
    expect_equal(length(out), 3)
    expect_equal(as.data.frame(out)$start, c(1,501,1,401,1))

    #### when ignore.strand = TRUE
    out <- suppressWarnings(.removeExonsfromTx(grl2, exonset3, ignore.strand = T))
    expect_equal(length(out), 3)
    expect_equal(as.data.frame(out)$start, c(1,501, 1, 401, 1,501))

    #### when drop.nonoverlap.tx = TRUE
    expect_warning(.removeExonsfromTx(grl2, exonset3, drop.unmodified = TRUE))
    out <- suppressWarnings(.removeExonsfromTx(grl2, exonset3, drop.unmodified = T))
    expect_equal(length(out), 1)
    expect_equal(as.data.frame(out)$start, as.data.frame(gr1)$start)

    #### when drop = TRUE & ignore.strand = TRUE
    out <- suppressWarnings(.removeExonsfromTx(grl2, exonset3,
                                              drop.unmodified = T, ignore.strand = T))
    expect_warning(.removeExonsfromTx(grl2, exonset3,
                                     drop.unmodified = T, ignore.strand = T))
    expect_equal(length(out), 2)
    expect_equal(as.data.frame(out)$start, c(1,501,1,501))
})

