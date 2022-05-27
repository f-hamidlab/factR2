#' @include generics.R

setMethod("show",
          "factR",
          function(object) {
              ngenes <- length(unique(object@txdata$gene_id))
              ntxs <- length(unique(object@txdata$transcript_id))
              nnovel <- sum(object@txdata$novel == "yes")
              ncds <- sum(object@txdata$cds == "yes")
              cat(sprintf("factR object [version %s]\n", object@version))
              cat(sprintf("## Total number of genes: %s\n", ngenes))
              cat(sprintf("## Total number of transcripts: %s [%s novel]\n", ntxs, nnovel))
              cat(sprintf("## Total number of CDSs: %s \n", ncds))
              cat(sprintf("## Reference species: %s\n", object@reference$species))
              cat(sprintf("## Reference database: %s\n", object@reference$db))
          }
)



setMethod("head", "factR", function(object, n = 6L) {
    return(object@custom[1:n])
})

setMethod("tail", "factR", function(object, n = 6L) {
    return(utils::tail(object@custom, n))
})

setMethod("View", "factR", function(object, ...,
                                    type = "transcripts",
                                    in_console = FALSE) {

    # check features
    genetxs <- methods::slot(object, name = "txdata")
    if(missing(...)){
        txs <- genetxs$transcript_id
    } else {
        genetxs.features <- genetxs %>%
            dplyr::mutate(tx = transcript_id) %>%
            tidyr::gather("type", "feature", gene_id, gene_name, transcript_id) %>%
            dplyr::filter(feature %in% c(...))
        if (nrow(genetxs.features) == 0) {
            rlang::abort("No features found in object")
        } else if(!all(c(...) %in% genetxs.features$feature)){
            absent.features <- c(...)[which(! c(...) %in% genetxs.features$feature)]
            rlang::warn(sprintf("These features are not found in object: %s",
                                paste(absent.features, collapse = ", ")))
        }
        txs <- genetxs.features$tx
    }

    # select features by data
    if(type == "transcripts"){
        gtf <- methods::slot(object, "custom")
        data.to.view <- gtf[gtf$transcript_id %in% txs]

    }

    # display data
    if(in_console){
        data.to.view
    } else {
        data.to.view <- as.data.frame(data.to.view)
        if(nrow(data.to.view) > 1000){
            rlang::warn("Viewing first 1000 entries")
            tibble::view(data.to.view[1:1000], title = "factRObject")
        } else {
            tibble::view(data.to.view, title = "factRObject")
        }
    }

})
































