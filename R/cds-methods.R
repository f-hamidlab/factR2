#' Build coding segments
#'
#' @description Constructs CDS information on transcripts from custom annotation
#' using reference-based approach.
#' 
#' @details
#' This function will firstly identify custom transcripts with identical 
#' exon structures as those in the reference annotation. If these reference
#' transripts contain CDS segments, the coordinates will be passed to its paired-
#' custom transcript. The putative start codon of the remaining transcripts
#' will be the first ATG sequence that is in-frame with the coding sequence
#' of reference transcripts and an in-frame stop codon will be determined.
#' 
#'
#' @param object factRObject
#' @param verbose Whether to print out messages (Default: FALSE)
#'
#' @return Updated factRObject with updated transcripts metadata.
#' @export
#' @seealso \code{\link{runfactR}}
#' @include factRObject-class.R
#' @rdname buildCDS
#' @examples
#' data(factRsample)
#' factRsample <- buildCDS(factRsample)
setGeneric("buildCDS", function(object, verbose = FALSE) standardGeneric("buildCDS"))


setMethod("buildCDS", "factR", function(object, verbose = FALSE) {

    if(verbose){.msgheader("Building CDS information")}
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














