#' @include generics.R
#'
setMethod("runfactR", "factR", function(object, verbose = FALSE) {
    object <- buildCDS(object, verbose)
    object <- predictNMD(object, verbose)
    object <- findAltSplicing(object)
    object <- getAAsequence(object, verbose)
    object
})


setMethod("buildCDS", "factR", function(object, verbose = FALSE) {
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

setMethod("predictNMD", "factR", function(object, NMD_threshold = 50, verbose = FALSE) {
    gtf <- slot(object, "custom")
    genetxs <- slot(object, "txdata")

    if(! "CDS" %in% gtf$type){
        rlang::abort("No CDSs found. Please run buildCDS() first")
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


setMethod("getAAsequence", "factR", function(object, verbose = FALSE) {
    gtf <- slot(object, "custom")
    if(! "CDS" %in% gtf$type){
        rlang::abort("No CDSs found. Please run buildCDS() first")
    }
    genetxs <- txData(object)
    txs <- genetxs[genetxs$cds == "yes",]$transcript_id

    gtf <- slot(object, "custom")
    gtf <- gtf[gtf$transcript_id %in% txs]
    cds <- S4Vectors::split(gtf[gtf$type == "CDS"], ~transcript_id)
    slot(object, "domains")$sequence <- .getSequence(cds,
                                                     slot(object, "reference")$genome,
                                                     verbose)

    return(object)
})









