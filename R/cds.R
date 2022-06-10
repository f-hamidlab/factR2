#' @include generics.R
#'
setMethod("buildCDS", "factR", function(object, verbose = FALSE) {

    gtf <- granges(object, set = "all") 
    gtf <- gtf[!gtf$type %in% "CDS"]

    if(verbose){
        out.gtf <- factR::buildCDS(gtf,
                               slot(object, "reference")$ranges,
                               slot(object, "reference")$genome)
    } else {
        out.gtf <- suppressMessages(
            factR::buildCDS(gtf,
                            slot(object, "reference")$ranges,
                            slot(object, "reference")$genome))
    }
    out.gtf <- out.gtf[out.gtf$type %in% "CDS"]
    gtf <- c(gtf, out.gtf)
    slot(object, "transcriptome") <- gtf

    # update cds transcripts
    cdss <- unique(gtf[gtf$type == "CDS"]$transcript_id)
    txs <- object[["transcript"]]
    object@sets$transcript@rowData$cds <- ifelse(txs$transcript_id %in% cdss,
                                                 "yes",
                                                 txs$cds)
    return(object)
})














