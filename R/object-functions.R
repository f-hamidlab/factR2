setMethod("show",
          "factR",
          function(object) {
              ngenes <- length(unique(object@custom$genetxs$gene_id))
              ntxs <- length(unique(object@custom$genetxs$transcript_id))
              nnovel <- sum(object@custom$genetxs$novel == "yes")
              cat(sprintf("factR object [version %s]\n", object@version))
              cat(sprintf("## Total number of genes: %s\n", ngenes))
              cat(sprintf("## Total number of transcripts: %s [%s novel]\n", ntxs, nnovel))
              cat(sprintf("## Reference species: %s\n", object@reference$species))
              cat(sprintf("## Reference database: %s\n", object@reference$db))
          }
)
