# **factR v.2**

## Functional Annotation of Custom Transcriptomes in R

## General workflow
<p align="center">
  <img src="man/figures/factR2.png" width="450"/>
</p>

**factR2** is an R package with user-friendly tools to work
with custom-assembled transcriptomes.  
It is a significant upgrade of its predecessor package, [**factR**](https://fursham-h.github.io/factR/).

Below are **factR2**'s key features:

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

## Quickstart on using factR2
factR2 requires a custom-assembled transcriptome in GTF format. 
Below is a quick-start on using factR2 using your own GTF file to output
(1) an updated GTF file with CDS and alternative splicing information,
(2) metadata at the gene, transcript and alternative splicing levels. 
```r
library(factR2)

# check supported genomes and take note of ID
listSupportedGenomes()

factR.object <- createfactRObject(gtf = "PATH/TO/GTF",
                                  reference = genome_ID)
factR.object <- runfactR(factR.object)

# export annotated GTF and metadata of genes,transcripts and alternative exons
exportAll(runfactR, path = "PATH/TO/OUTPUT")
```



