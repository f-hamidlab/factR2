#' @include generics.R
#'
setMethod("buildCDS", "factR", function(object, verbose = FALSE) {
    gtf <- slot(object, "transcriptome")

    if(verbose){
        gtf <- factR::buildCDS(gtf,
                               slot(object, "reference")$ranges,
                               slot(object, "reference")$genome)
    } else {
        gtf <- suppressMessages(
            factR::buildCDS(gtf,
                            slot(object, "reference")$ranges,
                            slot(object, "reference")$genome))
    }

    slot(object, "transcriptome") <- gtf

    # update cds transcripts
    cdss <- unique(gtf[gtf$type == "CDS"]$transcript_id)
    genetxs <- slot(object, "txData")
    slot(object, "txData")$cds <- ifelse(genetxs$transcript_id %in% cdss,
                                                 "yes",
                                                 genetxs$cds)
    return(object)
})














