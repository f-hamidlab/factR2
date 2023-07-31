#' Predict NMD-sensitive transcripts
#'
#' @description Annotates NMD-sensitive transcripts based on Exon-Junction model.
#' By default, stop codons which are >50bp upstream of the last exon junction will
#' be annotated as NMD-sensitive
#'
#' @param object factRObject
#' @param NMD_threshold Minimum distance between PTC and downstream exon-exon junction to trigger NMD (Default: 50)
#' @param verbose Whether to print out messages (Default: FALSE)

#' @export
#' @seealso \code{\link{runfactR}}
#' @include factRObject-class.R
#' @rdname predictNMD
#' @examples
#' data(factRsample)
#' factRsample <- buildCDS(factRsample)
#' factRsample <- predictNMD(factRsample)
setGeneric("predictNMD", function(object, NMD_threshold = 50, verbose = FALSE) standardGeneric("predictNMD"))



setMethod("predictNMD", "factR", function(object, NMD_threshold = 50, verbose = FALSE) {

    if(verbose){.msgheader("Running transcript-level NMD prediction")}

    gtf <- slot(object, "transcriptome")
    gtf <- gtf[!gtf$type %in% c("gene", "AS")]

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

    if(verbose){.msgsubinfo("Updating transcript feature data")}
    nmd.out <- nmd.out %>%
        tibble::column_to_rownames("transcript")
    object <- addMeta(object, meta = "transcript",
                      data = nmd.out,
                      nmd = ifelse(is_NMD & !is.na(is_NMD),
                                   "yes", nmd))

    return(object)
})



#' Identify AS-NMD events
#'
#' @description Pinpoints alternative splicing events underlying NMD sensitivity.
#' Metadata of AS events will be updated to include columns "ASNMDtype" and "ASNMD.in.cds".
#' Values of "ASNMDtype" can be Repressing (skipping of exon causes transcripts to be NMD-sensitive),
#' Stimulating (splicing of exon causes transcripts to be NMD-sensitive) or
#' NA (does not cause transcripts to NMD-sensitive).
#'
#' @param object factRObject
#' @param verbose Whether to print out messages (Default: FALSE)
#'
#' @return factRObject with updated ASE metadata
#' @export
#' @seealso \code{\link{runfactR}}
#'
#' @rdname testASNMDevents
#' @examples
#' data(factRsample)
#' factRsample <- buildCDS(factRsample)
#' factRsample <- testASNMDevents(factRsample)
setGeneric("testASNMDevents", function(object, verbose = FALSE) standardGeneric("testASNMDevents"))
setMethod("testASNMDevents", "factR", function(object, verbose = FALSE) {

    if(verbose){.msgheader("Running AS-NMD testing")}

    genes <- object[["transcript"]]
    gtf <- object@transcriptome
    if(all(genes$cds == "no") & all(genes$nmd == "no")){
        rlang::abort("No CDSs found. Please run runfactR() first")
    } else if(any(genes$cds == "yes") & all(genes$nmd == "no")){
        rlang::abort("NMD not predicted. Please run predictNMD() first")
    }

    # check if splicing data is present
    ASevents <- gtf[gtf$type == "AS"]
    if(length(ASevents) == 1){
        rlang::abort("No AS events found. Please run findAltSplicing() first")
    }


    # check input objects
    #phastGScore <- .GScorecheck(ConsScores)

    # get reference CDS transcript for each gene

    if(verbose){ .msgsubinfo("Getting best reference for AS-NMD testing")}
    ref <- .getbestref(gtf, genes)

    # run core function
    if(verbose){ .msgsubinfo("Testing AS-NMD exons")}
    ASNMDevents <- .runidentifynmdexons(ASevents, genes, ref,
                                        object@reference$genome, gtf,
                                        verbose)

    # update AS events if returned object is not null
    if(!is.null(ASNMDevents)){
        ASevents$ASNMDtype <- ifelse(ASevents$AS_id %in% rownames(ASNMDevents),
                                     ASNMDevents[ASevents$AS_id,]$NMDtype, "NA")
        ASevents$ASNMD.in.cds <- ifelse(ASevents$AS_id %in% rownames(ASNMDevents),
                                        ASNMDevents[ASevents$AS_id,]$within.CDS, "NA")
        gtf.others <- gtf[gtf$type != "AS"]
        object@transcriptome <- c(gtf.others, ASevents)

        # update featureData
        if(verbose){ .msgsubinfo("Updating AS feature data")}
        ASNMDevents <- ASNMDevents %>%
            as.data.frame() %>%
            dplyr::select(ASNMDtype=NMDtype, ASNMD.in.cds=within.CDS)

        object <- addMeta(object, "AS", data = ASNMDevents)

    }

    return(object)
})


.getbestref <- function(x, genes) {
    #rlang::inform("Selecting best reference mRNAs")
    # get reference CDS transcript for each gene
    ## get sizes of cdss
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




.runidentifynmdexons <- function(ASevents, genes, ref, genome, gtf, verbose) {

    #rlang::inform("Finding NMD causing exons")
    # NMD.pos <- genes[genes$nmd == "yes",]$transcript_id
    #
    # # get AS segments between NMD transcript and reference transcript
    # ## shortlist NMD transcripts from GTF
    # x.NMD <- x[x$transcript_id %in% NMD.pos]
    #
    # # shortlist NMD transcripts from genes containing reference mRNAs
    # x.NMD <- x.NMD[x.NMD$gene_id %in% ref$gene_id]

    ## remove transcripts that are not NMD causing
    NMD.pos <- genes[genes$nmd == "yes",]$transcript_id
    NMDevents <- ASevents[ASevents$transcript_id %in% NMD.pos]
    NMDevents <- NMDevents[NMDevents$gene_id %in% ref$gene_id]

    ## get AS segments and annotate its splicing nature
    ref_overlaps <- IRanges::findOverlapPairs(NMDevents, ref[ref$type == "exon"]) %>%
        as.data.frame() %>%
        dplyr::filter(first.gene_id == second.gene_id) %>%
        dplyr::filter(first.X.start >= second.X.start) %>%
        dplyr::filter(first.X.end <= second.X.end) %>%
        dplyr::pull(first.X.AS_id)
    ASspliced <- NMDevents[!NMDevents$AS_id %in% ref_overlaps] %>%
        as.data.frame() %>%
        dplyr::mutate(coord = paste0(seqnames, "_", start, "_",
                                     end, "_", strand, "_", AStype)) %>%
        dplyr::select(AS_id, gene_id, transcript_id, coord) %>%
        dplyr::mutate(splice = "spliced")

    ## get AS segments missing from ref
    AStogene <- NMDevents %>%
        as.data.frame() %>%
        dplyr::distinct(transcript_id, gene_id)

    AS.grl <- S4Vectors::split(gtf[gtf$type == "exon"], ~transcript_id)
    AS.grl <- AS.grl[AStogene$transcript_id]

    ref.grl <- S4Vectors::split(ref[ref$type == "exon"], ~gene_id)
    ref.grl <- ref.grl[AStogene$gene_id]
    names(ref.grl) <- AStogene$transcript_id


    ASskipped <- GenomicRanges::setdiff(ref.grl, AS.grl) %>%
        as.data.frame() %>%
        dplyr::left_join(ASevents %>%
                             as.data.frame() %>%
                             dplyr::select(seqnames, start, end, strand, gene_id, AStype, AS_id) %>%
                             dplyr::distinct(),
                         by = c("seqnames", "start", "end", "strand")) %>%
        dplyr::filter(!is.na(AS_id)) %>%
        dplyr::filter(gene_id %in% ref$gene_id) %>%
        dplyr::mutate(coord = paste0(seqnames, "_", start, "_",
                                     end, "_", strand, "_", AStype)) %>%
        dplyr::select(AS_id, gene_id, transcript_id = group_name, coord) %>%
        dplyr::mutate(splice = "skipped")

    ASevents.fulltx <- dplyr::bind_rows(ASskipped, ASspliced)

    # return if no AS exons were found
    if(nrow(ASevents.fulltx) == 0) {
        if(verbose){ .msgsubwarn("No alternatively spliced exons found")}
        return(NULL)
    }

    # simplify df by removing redundant exons
    ## in cases where ref and NMD tx contain said exon, the skipped form will be retained
    ASevents <- ASevents.fulltx %>%
        dplyr::distinct(AS_id, .keep_all = TRUE)

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
    mod.tx <- suppressMessages(factR::buildCDS(mod.tx, ref, genome))
    mod.NMD <- suppressMessages(factR::predictNMD(mod.tx, progress_bar = FALSE))

    # return if none of the reconstructed transcripts are NMD sensitive
    if(sum(mod.NMD$is_NMD) == 0) {
        if(verbose){ .msgsubwarn("None of the alternative exons are NMD-causing")}
        return(NULL)
    }



    # Report NMD exons
    ASevents <- ASevents.fulltx %>%
        dplyr::left_join(mod.NMD %>% dplyr::select(transcript, is_NMD),
                         by = c("coord"="transcript")) %>%
        dplyr::mutate(NMDtype = ifelse(splice == "skipped", "Repressing", "Stimulating")) %>%
        dplyr::select(coord, AS_id, gene_id, transcript_id, NMDtype, is_NMD) %>%
        tidyr::separate(coord, c("seqnames", "start", "end", "strand", "AStype"),
                        sep = "_", convert = TRUE) %>%
        dplyr::filter(is_NMD) %>%
        dplyr::group_by(transcript_id) %>%
        dplyr::arrange(ifelse(strand == "-", dplyr::desc(start), start)) %>%
        dplyr::distinct(transcript_id, .keep_all = T) %>%
        dplyr::ungroup() %>%
        dplyr::select(seqnames:strand, gene_id, AStype, AS_id, NMDtype) %>%
        dplyr::distinct(AS_id, .keep_all = T) %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

    # Annotate location of NMD events relative to reference CDS
    ## Useful to know if exons are 3'UTR introns etc
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
    out.df <- S4Vectors::mcols(ASevents)[c("AS_id", "NMDtype", "within.CDS")]
    rownames(out.df) <- out.df$AS_id

    return(out.df)
}









