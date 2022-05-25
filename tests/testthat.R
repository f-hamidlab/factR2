library(testthat)
library(factR2)
library(BSgenome.Mmusculus.UCSC.mm10)

data(genomes)
gtf <- system.file("extdata", "sc_merged_sample.gtf.gz", package = "factR")
test_check("factR2")
