# **factR v.2**

## Functional Annotation of Custom Transcriptomes in R

## General workflow
<p align="center">
  <img src="man/figures/factR2.png" width="450"/>
</p>


**factR2** represents a substantial improvement over 
[**factR**](https://fursham-h.github.io/factR/), providing users with
more powerful and user-friendly tools to work
with custom-assembled transcriptomes.  

Below are **factR2**'s latest features:

1. Extracts alternative splicing events and annotate its novelity and 
contribution to NMD outcome
2. Tests regulatory potential of AS-NMD events
3. Quantifies evolutionary conservation scores of alternative exons
4. Interactive plot of transcript architectures

As well as the following features from **factR**:

1. Matches gene information on custom transcriptomes to reference annotation
2. Constructs transcript coding (CDS) information using reference-guided approach
3. Predicts sensitivity of coding transcripts to nonsense-mediated decay (NMD)
4. Predicts protein domains on productively spliced transcripts
  
## How to install
The development version can be installed using devtools:
```r
# install.packages("devtools")
devtools::install_github("f-hamidlab/factR2")
```

## Getting started

See our 
[full vignette](https://htmlpreview.github.io/?https://github.com/f-hamidlab/factR2/blob/master/vignettes/factR2.html) 
on how to use **factR2**.
