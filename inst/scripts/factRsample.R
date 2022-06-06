# creating factRsample
## get custom transcriptome
path.to.gtf <- system.file("extdata/sc_merged_sample.gtf.gz", package = "factR")
gtf <- rtracklayer::import(path.to.gtf)
gtf <- gtf[!grepl("^Gm", gtf$gene_name)]
gtf <- gtf[gtf$gene_id %in% unique(gtf$gene_id)[1:50]]

## get ref transcriptome
ref <- rtracklayer::import("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz")
ref <- ref[as.character(GenomicRanges::seqnames(ref)) %in% "chr15"]
GenomeInfoDb::seqlevels(ref) <- "chr15"
ref <- IRanges::subsetByOverlaps(ref, range(gtf))

## get fasta
fa <- factR::importFASTA("/media/cdn-bc/RAID/Genomes/GRCm38_mm10/Gencode_vM25/GRCm38.primary_assembly.genome.fa.gz")
fa <- fa["chr15"]
fa$chr15 <- fa$chr15[1:12500000]

## make sample data
factRsample <- createfactRObject(gtf, use_own_annotation = ref,
                                 use_own_genome = fa)
usethis::use_data(factRsample, overwrite = TRUE)
