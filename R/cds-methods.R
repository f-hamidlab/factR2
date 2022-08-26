#' @include factRObject-class.R
#'

setGeneric("buildCDS", function(object, verbose = FALSE) standardGeneric("buildCDS"))
setMethod("buildCDS", "factR", function(object, verbose = FALSE) {

    .msgheader("Building CDS information")
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

    if(verbose){.msgsubinfo("Updating transcript feature data")}
    cdss <- unique(gtf[gtf$type == "CDS"]$transcript_id)
    object <- addMeta(object, meta = "transcript",
                      cds = ifelse(transcript_id %in% cdss, "yes", cds))

    return(object)
})














