# factR2 0.99.0

* First release of factR2

# factR2 0.99.1

* Added functionalities to work with expression levels
    * `addTxCounts()` : Imports transcript-level counts and calculates normalised
    gene-level and transcript-level expressions. Will also deduce the splicing
    efficiencies of alternative exons (PSI values) based on transcript-level
    expression.
    * `testGeneCorr()` : Tests the causality of gene expression changes against
    the splicing of AS-NMD exons to provide insights into the regulatory function
    of  alternative exons.
    * `plot2way()` : Plots the values from any 2 features as a scatterplot.
    
# factR2 0.99.2

* `getAScons()` now computes conservation scores of flanking intronic segments.
* `getAScons()` will now import pre-downloaded GSscore databases and prompts user
to download available databases.

    
