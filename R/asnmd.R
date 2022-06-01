#' @include generics.R
#'
setMethod("predictNMD", "factR", function(object, NMD_threshold = 50, verbose = FALSE) {
    gtf <- slot(object, "transcriptome")
    gtf <- gtf[gtf$type != "AS"]
    genetxs <- txData(object)

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


    genetxs <- dplyr::left_join(genetxs, nmd.out, by = c("transcript_id"="transcript"))
    genetxs$nmd <- ifelse(genetxs$is_NMD & !is.na(genetxs$is_NMD), "yes", genetxs$nmd)
    slot(object, "txData") <- genetxs
    return(object)
})

setMethod("testASNMDevents", "factR", function(object) {

    genes <- txData(object)
    gtf <- slot(object, "transcriptome")
    if(all(genes$cds == "no") & all(genes$nmd == "no")){
        rlang::abort("No CDSs found. Please run runfactR() first")
    } else if(any(genes$cds == "yes") & all(genes$nmd == "no")){
        rlang::abort("No CDSs found. Please run predictNMD() first")
    }

    # check if splicing data is present
    ASevents <- gtf[gtf$type == "AS"]
    if(length(ASevents) == 1){
        rlang::abort("No AS events found. Please run findAltSplicing() first")
    }


    # check input objects
    #phastGScore <- .GScorecheck(ConsScores)

    # get reference CDS transcript for each gene
    ref <- .getbestref(object)

    # run core function
    ASNMDevents <- .runidentifynmdexons(object, ref)

    # update AS events if returned object is not null
    if(!is.null(ASNMDevents)){
        ASNMD.hits <- IRanges::findOverlaps(ASevents, ASNMDevents,
                                            type = "equal")
        ASevents$ASNMDtype <- "NA"
        ASevents[S4Vectors::queryHits(ASNMD.hits)]$ASNMDtype <- ASNMDevents[S4Vectors::subjectHits(ASNMD.hits)]$NMDtype
        ASevents$ASNMD.in.cds <- "NA"
        ASevents[S4Vectors::queryHits(ASNMD.hits)]$ASNMD.in.cds <- as.character(ASNMDevents[S4Vectors::subjectHits(ASNMD.hits)]$within.CDS)
        gtf.others <- gtf[gtf$type != "AS"]
        slot(object, "transcriptome") <- c(gtf.others, ASevents)
    }

    return(object)
})


.getbestref <- function(object) {
    #rlang::inform("Selecting best reference mRNAs")
    # get reference CDS transcript for each gene
    ## get sizes of all CDSs
    x <- slot(object, "transcriptome")
    genes <- txData(object)

    cds.sizes <- sum(BiocGenerics::width(S4Vectors::split(x[x$type == "CDS"],
                                                          ~transcript_id)))
    ## shortlist non NMD transcripts
    cds.reference <- genes %>%
        dplyr::filter(nmd == "no", novel == "no", cds == "yes") %>%
        dplyr::left_join(as.data.frame(cds.sizes) %>%
                             tibble::rownames_to_column("transcript_id"),
                         by = "transcript_id") %>%
        dplyr::group_by(gene_id) %>%
        dplyr::slice_max(n = 1, order_by = cds.sizes, with_ties = FALSE)


    return(x[x$transcript_id %in% cds.reference$transcript_id & x$type != "AS"])
}




.runidentifynmdexons <- function(object, ref) {

    #rlang::inform("Finding NMD causing exons")
    x <- slot(object, "transcriptome")
    ASevents <- x[x$type == "AS"]
    genes <- txData(object)
    NMD.pos <- genes[genes$nmd == "yes",]$transcript_id

    # get AS segments between NMD transcript and reference transcript
    ## shortlist NMD transcripts from GTF
    x.NMD <- x[x$transcript_id %in% NMD.pos]

    # shortlist NMD transcripts from genes containing reference mRNAs
    x.NMD <- x.NMD[x.NMD$gene_id %in% ref$gene_id]

    ## get AS segments and annotate its splicing nature
    ASevents$splice <- ifelse(ASevents$transcript_id %in% ref$transcript_id,
                              "skipped", "spliced")
    ASevents <- ASevents[ASevents$gene_id %in% ref$gene_id]

    # return if no AS exons were found
    if(length(ASevents) == 0) {
        rlang::warn("No alternatively spliced exons found")
        return(NULL)
    }

    # simplify df by removing redundant exons
    ## in cases where ref and NMD tx contain said exon, the skipped form will be retained
    ASevents <- ASevents %>%
        as.data.frame() %>%
        dplyr::mutate(coord = paste0(seqnames, "_", start, "_",
                                     end, "_", strand, "_", AStype)) %>%
        dplyr::arrange(splice) %>%
        dplyr::select(gene_id, transcript_id, coord, splice) %>%
        dplyr::distinct(coord, .keep_all = TRUE)

    # recreate hypothetical tx by inserting/removing AS segments
    ## create grl of ref trancripts
    ref.grl <- S4Vectors::split(ref[ref$type == "exon"], ~gene_id)
    mod.tx <- tibble::tibble(
        "seqnames" = as.character(),
        "start" = as.integer(),
        "end" = as.integer(),
        "strand" = as.character(),
        "coord" = as.character()
    )

    ## Add exons to ref transcripts
    if("spliced" %in% ASevents$splice) {
        toadd <- dplyr::filter(ASevents, splice == "spliced")
        addedtx <- ref.grl[toadd$gene_id]
        names(addedtx) <- toadd$coord
        exons.toadd <- toadd["coord"] %>%
            tidyr::separate(coord, c("seqnames", "start", "end", "strand", "AStype"),
                            sep = "_") %>%
            GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
        addedtx <- .addExonstoTx(addedtx, exons.toadd) %>%
            as.data.frame() %>%
            dplyr::select(seqnames, start, end, strand, coord = group_name)
        mod.tx <- dplyr::bind_rows(mod.tx, addedtx)
    }


    ## Remove exons from ref transcripts
    if("skipped" %in% ASevents$splice) {
        toremove <- dplyr::filter(ASevents, splice == "skipped")
        removedtx <- ref.grl[toremove$gene_id]
        names(removedtx) <- toremove$coord
        exons.toremove <- toremove["coord"] %>%
            tidyr::separate(coord, c("seqnames", "start", "end", "strand", "AStype"),
                            sep = "_") %>%
            GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
        removedtx <- .removeExonsfromTx(removedtx, exons.toremove) %>%
            as.data.frame() %>%
            dplyr::select(seqnames, start, end, strand, coord = group_name)
        mod.tx <- dplyr::bind_rows(mod.tx, removedtx)
    }


    # cleanup mod.tx by adding gene_id info and create GRanges
    mod.tx <- mod.tx %>%
        dplyr::left_join(ASevents %>% dplyr::select(coord, gene_id), by = "coord") %>%
        dplyr::mutate(transcript_id = coord, type = "exon") %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

    # test NMD again
    mod.tx <- suppressMessages(factR::buildCDS(mod.tx, ref, object@reference$genome))
    mod.NMD <- suppressMessages(factR::predictNMD(mod.tx, progress_bar = FALSE))

    # return if none of the reconstructed transcripts are NMD sensitive
    if(sum(mod.NMD$is_NMD) == 0) {
        rlang::warn("None of the alternative exons are NMD-causing")
        return(NULL)
    }

    # Report NMD exons
    ASevents <- ASevents %>%
        dplyr::left_join(mod.NMD %>% dplyr::select(transcript, is_NMD),
                         by = c("coord"="transcript")) %>%
        dplyr::mutate(NMDtype = ifelse(splice == "skipped", "Repressing", "Stimulating")) %>%
        dplyr::select(coord, gene_id, NMDtype, is_NMD) %>%
        tidyr::separate(coord, c("seqnames", "start", "end", "strand", "AStype"),
                        sep = "_") %>%
        dplyr::filter(is_NMD) %>%
        dplyr::select(seqnames:strand, gene_id, AStype, NMDtype) %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    ASevents$is_NMD <- NULL

    # Annotate location of NMD events relative to reference CDS
    ## Useful to know if exons are 3'UTR introns etc
    ref <- .getbestref(object)
    ref.cds.grl <- S4Vectors::split(ref[ref$type == "CDS"], ~gene_id)
    ASevents$within.CDS <- IRanges::overlapsAny(ASevents, range(ref.cds.grl))


    ## run conservation scoring
    # if(typeof(phast) == "S4"){
    #     rlang::inform("Quantifying exon conservation scores")
    #     ASevents <- tryCatch(
    #         {
    #             ASevents <- GenomicScores::gscores(phast, ASevents)
    #             BiocGenerics::colnames(S4Vectors::mcols(ASevents))[5] <- phast@data_pkgname
    #             return(ASevents)
    #         },
    #         error = function(cond){
    #             rlang::warn("Unable to quantify exon conservation scores. Check if annotation package matches the genome used")
    #             return(ASevents)
    #         }
    #     )
    #
    #
    # }

    return(ASevents)
}









