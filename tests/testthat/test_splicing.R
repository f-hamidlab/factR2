#context("Test factR object")

# setup GRanges
gr1 <- GenomicRanges::GRanges(
    seqnames = "1", strand = rep("+", 5),
    ranges = IRanges::IRanges(
        start = c(1, 401, 801, 1001, 1301),
        end = c(100, 500, 900, 1100, 1400)
    ),
    transcript_id = "transcript1",
    gene_id = "geneA",
    gene_name = "gene001"
)
gr2 <- GenomicRanges::GRanges(
    seqnames = "1", strand = rep("+", 6),
    ranges = IRanges::IRanges(
        start = c(201, 401, 601, 801, 981, 1201),
        end = c(300, 500, 700, 920, 1100, 1250)
    ),
    transcript_id = "transcript2",
    gene_id = "geneA",
    gene_name = "gene001"
)
gr3 <- GenomicRanges::GRanges(
    seqnames = "1", strand = rep("+", 2),
    ranges = IRanges::IRanges(
        start = c(350, 801),
        end = c(500, 1150)
    ),
    transcript_id = "transcript3",
    gene_id = "geneA",
    gene_name = "gene002"
)

gr4 <- GenomicRanges::GRanges(
    seqnames = "1", strand = rep("+", 2),
    ranges = IRanges::IRanges(
        start = c(525, 725),
        end = c(550, 750)
    ),
    transcript_id = "transcript4",
    gene_id = "geneA",
    gene_name = "gene001"
)

gr5 <- GenomicRanges::GRanges(
    seqnames = "1", strand = rep("+", 4),
    ranges = IRanges::IRanges(
        start = c(450, 801, 1001, 1301),
        end = c(500, 900, 1100, 1380)
    ),
    transcript_id = "transcript5",
    gene_id = "geneA",
    gene_name = "gene001"
)

gr6 <- GenomicRanges::GRanges(
    seqnames = "1", strand = rep("+", 3),
    ranges = IRanges::IRanges(
        start = c(10, 401, 801),
        end = c(100, 500, 850)
    ),
    transcript_id = "transcript6",
    gene_id = "geneA",
    gene_name = "gene001"
)

grcomb1 <- c(gr1, gr2, gr3)


test_that("Test .runAS functionality", {
    out <- .runAS(grcomb1)
    expect_equal(out$AStype, c("AF", "CE", "AD", "RI", "AA", "AL"))
    expect_equal(GenomicRanges::start(out), c(201, 601, 901, 901, 981, 1201))
    expect_equal(GenomicRanges::width(out), c(100, 100, 20, 100, 20, 50))

    GenomicRanges::strand(grcomb1) <- "-"
    out <- .runAS(grcomb1)
    expect_equal(out$AStype, c("AL", "CE", "AA", "RI", "AD", "AF"))
})

test_that("Test FindAltSplicing generics", {
    out <- FindAltSplicing(factRsample)
    expect_equal(as.character(class(out@ASplicings$ranges)), "GRanges")
    expect_equal(length(out@ASplicings$ranges), 143)
})

