#' @include generics.R
#'
setMethod("RunfactR", "factR", function(object, verbose = FALSE) {
    object <- BuildCDS(object, verbose)
    object <- PredictNMD(object, verbose)
    object <- FindAltSplicing(object)
    object
})


setMethod("BuildCDS", "factR", function(object, verbose = FALSE) {
    gtf <- slot(object, "custom")

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

    slot(object, "custom") <- gtf

    # update cds transcripts
    cdss <- unique(gtf[gtf$type == "CDS"]$transcript_id)
    genetxs <- slot(object, "txdata")
    slot(object, "txdata")$cds <- ifelse(genetxs$transcript_id %in% cdss,
                                                 "yes",
                                                 genetxs$cds)
    return(object)
})

setMethod("PredictNMD", "factR", function(object, NMD_threshold = 50, verbose = FALSE) {
    gtf <- slot(object, "custom")
    genetxs <- slot(object, "txdata")

    if(! "CDS" %in% gtf$type){
        rlang::abort("No CDSs found. Please run BuildCDS() first")
    }

    if(verbose){
        nmd.out <- factR::predictNMD(gtf, NMD_threshold = NMD_threshold,
                                     progress_bar = TRUE)
    } else {
        nmd.out <- suppressMessages(
            factR::predictNMD(gtf, NMD_threshold = NMD_threshold,
                              progress_bar = FALSE))
    }

    slot(object, "nmd") <- nmd.out
    nmd.out.true <- nmd.out[nmd.out$is_NMD,]$transcript
    slot(object, "txdata")$nmd <- ifelse(genetxs$transcript_id %in% nmd.out.true,
                                                 "yes",
                                                 genetxs$nmd)
    return(object)
})











