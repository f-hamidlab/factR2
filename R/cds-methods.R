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
#' @return Updated factRObject.
#' The GTF GenomicRanges object will now contain "CDS" type features.
#'
#' The `cds` variable from the transcript metadata will return "yes" for coding
#' transcripts and "no" for non-coding transcripts
#'
#'
#' @export
#' @seealso \code{\link{runfactR}}
#' @include factRObject-class.R
#' @rdname buildCDS
#' @examples
#' ## Load sample data
#' data(factRsample)
#'
#' ## Build coding sequences
#' factRsample <- buildCDS(factRsample)
setGeneric("buildCDS", function(object, verbose = FALSE) standardGeneric("buildCDS"))


setMethod("buildCDS", "factR", function(object, verbose = FALSE) {

    if(verbose){.msgheader("Building CDS information")}
    gtf <- granges(object, set = c("gene", "transcript", "exon"))

    if(verbose){
        out.gtf <- suppressWarnings(factR::buildCDS(gtf,
                               slot(object, "reference")$ranges,
                               slot(object, "reference")$genome))
    } else {
        out.gtf <- suppressMessages(suppressWarnings(
            factR::buildCDS(gtf,
                            slot(object, "reference")$ranges,
                            slot(object, "reference")$genome)))
    }
    out.gtf <- out.gtf[out.gtf$type %in% "CDS"]

    # fix somemetadata in out.gtf
    meta <- gtf[gtf$type == "transcript"]@elementMetadata
    rownames(meta) <- meta$transcript_id
    out.gtf$match_level <- meta[out.gtf$transcript_id,]$match_level
    out.gtf$old_gene_id <- meta[out.gtf$transcript_id,]$old_gene_id

    # append start_codon and stop_codon information
    cds <- S4Vectors::split(out.gtf, ~transcript_id)
    cds.lengths <- sum(IRanges::width(cds))
    cds.start <- unlist(factR::trimTranscripts(cds, end = cds.lengths - 3))
    names(cds.start) <- NULL
    cds.start$type <- "start_codon"
    cds.end <- unlist(factR::trimTranscripts(cds, start = cds.lengths - 3))
    names(cds.end) <- NULL
    cds.end$type <- "stop_codon"

    gtf <- .sortGTF(c(granges(object, set = c("gene", "transcript", "exon", "AS")),
                      out.gtf, cds.start, cds.end))
    slot(object, "transcriptome") <- gtf

    # update cds transcripts

    if(verbose){.msgsubinfo("Updating transcript feature data")}
    cdss <- unique(gtf[gtf$type == "CDS"]$transcript_id)
    object <- addMeta(object, meta = "transcript",
                      cds = ifelse(transcript_id %in% cdss, "yes", cds))

    return(object)
})














