# **factR v.2**

## Functional Annotation of Custom Transcriptomes in R

## General workflow
<p align="center">
  <img src="man/figures/factR2.png" width="450"/>
</p>

**factR2** is a significant upgrade of its predecessor 
[**factR**](https://fursham-h.github.io/factR/) package with user-friendly tools to work 
with custom-assembled transcriptomes. Below are **factR2**'s key features:

* Core features
  1. Matches gene information on custom transcriptomes to reference annotation
  2. Constructs transcript coding (CDS) information using reference-guided approach
  3. Predicts sensitivity of coding transcripts to nonsense-mediated decay (NMD)
  4. Extracts alternative splicing events and annotate its novelity and contribution to NMD outcome
  4. Tests regulatory potential of AS-NMD events
  5. Quantifies evolutionary conservation scores of alternative exons

* Supporting features 
  1. Predicts protein domains on productively spliced transcripts
  2. Plots transcript architectures 
  3. Plots protein domain architectures
  
## How to install
The development version can be installed using devtools:
```r
# install.packages("devtools")
devtools::install_github("f-hamidlab/factR2")
```

## How to use **factR2**

**factR2** requires a custom-assembled transcriptome in GTF format. Below is 
a quick-start on using **factR2** using a sample GTF file:

```r
library(factR2)

# get path to sample GTF file
gtf <- system.file("extdata/sc_merged_sample.gtf.gz", package = "factR")

# create factRObject
factR.object <- createfactRObject(gtf, "vM25")

# interacting with factRObject
## print out metadata of various levels
genes(factR.object)  # gene metadata
txs(factR.object)  # transcript metadata
ase(factR.object)  # alternative splicing metadata

## print out transcript metadata for particular genes
txs(factR.object, "Dab2")
txs(factR.object, "Dab2", "Osmr")  # accepts multiple input genes

## View transcript isoforms
plotTranscripts(factR.object, "Dab2")  #from specific genes
plotTranscripts(factR.object, "ENSMUST00000078019.12") # of specific transcripts
### to rescale introns (useful for transcripts with long introns)
plotTranscripts(factR.object, "ENSMUST00000078019.12", rescale_introns=TRUE) 
plotTranscripts(factR.object, "AS00179")  # of specific AS events


# run main factR pipeline
factR.object <- runfactR(factR.object)

# export alternative splicing metadata
exportTable(factR.object, data = "AS")

```

