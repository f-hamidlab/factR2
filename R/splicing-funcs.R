.findAS <-  function(object) {
    gtf <- slot(object, "transcriptome")
    gtf <- c(gtf[gtf$type != "AS"], .runAS(gtf[gtf$type == "exon"]))
    slot(object, "transcriptome") <- gtf
    object@sets$AS@rowData <- as.data.frame(object@transcriptome) %>%
      dplyr::filter(type %in% "AS") %>%
      dplyr::mutate(coord = paste0(seqnames, ":", start, "-", end)) %>%
      dplyr::select(AS_id, gene_id, gene_name, coord, AStype, strand, width) %>%
      dplyr::distinct() %>%
      dplyr::mutate(AStype = factor(AStype, levels = c("CE", "AD","AA","AF","AL","RI")))
    return(object)
}


.runAS <- function(x) {
    # define global variables
    transcript_id <- pos <- seqnames <- strand <- gene_id <- NULL
    first.X.start <- second.X.start <- first.X.end <- second.X.end <- NULL
    first.pos <- AStype <- first.X.strand <- gene_name <- termini <- NULL

    # order exons by chromosome coord and label position
    ##### mutate group name herr
    x <- x %>%
        as.data.frame() %>%
        dplyr::group_by(transcript_id) %>%
        dplyr::arrange(start) %>%
        dplyr::mutate(pos = dplyr::row_number()) %>%
        dplyr::mutate(pos = ifelse(pos == 1, "First", pos)) %>%
        dplyr::mutate(pos = ifelse(pos == dplyr::n(), "Last", pos)) %>%
        dplyr::select(seqnames, start, end, strand, gene_id, gene_name, transcript_id, pos) %>%
        dplyr::group_by(gene_id) %>%
        dplyr::arrange(start, end) %>%
        dplyr::mutate(termini = dplyr::row_number()) %>%
        dplyr::mutate(termini = ifelse(termini == 1 | termini == dplyr::n(), T, F)) %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)

    # trim-off TS and TE
    x <- x %>%
        as.data.frame() %>%
        dplyr::group_by(seqnames, end, gene_id) %>%
        dplyr::mutate(start = ifelse(pos == "First", start[dplyr::n()], start)) %>%
        dplyr::group_by(seqnames, start, gene_id) %>%
        dplyr::mutate(end = ifelse(pos == "Last", end[1], end)) %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)


    # get gene sizes
    t2g <- x %>%
        as.data.frame() %>%
        dplyr::select(gene_id, transcript_id) %>%
        dplyr::distinct()

    # get reduced intron boundaries
    exonsbytx <- S4Vectors::split(x, ~transcript_id)
    intronsbytx <- GenomicRanges::psetdiff(BiocGenerics::unlist(range(exonsbytx)), exonsbytx)

    # regroup introns by gene_ids
    intronsbygroup <- t2g %>%
        dplyr::left_join(as.data.frame(intronsbytx), by = c("transcript_id" = "group_name")) %>%
        dplyr::filter(!is.na(seqnames)) %>%
        GenomicRanges::makeGRangesListFromDataFrame(split.field = "gene_id") %>%
        GenomicRanges::reduce() %>%
        as.data.frame() %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)


    # get exons that overlap with reduced introns
    altexons <- IRanges::findOverlapPairs(x, intronsbygroup)

    # annotate internal AS exons
    altannotate <- altexons %>%
        as.data.frame() %>%
        dplyr::mutate(AStype = "CE") %>%
        dplyr::mutate(AStype = ifelse(first.X.start < second.X.start & first.X.end < second.X.end, "Ad", AStype)) %>%
        dplyr::mutate(AStype = ifelse(first.X.start > second.X.start & first.X.end > second.X.end, "Aa", AStype)) %>%
        dplyr::mutate(AStype = ifelse(first.X.start < second.X.start & first.X.end <= second.X.end & first.pos == "Last", NA, AStype)) %>%
        dplyr::mutate(AStype = ifelse(first.X.start >= second.X.start & first.X.end > second.X.end & first.pos == "First", NA, AStype)) %>%
        dplyr::mutate(AStype = ifelse(first.X.start < second.X.start & first.X.end > second.X.end, "RI", AStype)) %>%
        dplyr::mutate(AStype = ifelse(first.X.start >= second.X.start & first.X.end < second.X.end & first.pos == "First", "Af", AStype)) %>%
        dplyr::mutate(AStype = ifelse(first.X.start > second.X.start & first.X.end <= second.X.end & first.pos == "Last", "Al", AStype)) %>%
        dplyr::mutate(AStype = ifelse(first.X.strand == "-", chartr("adfl", "dalf", AStype), AStype)) %>%
        dplyr::mutate(same.group = ifelse(first.gene_id == second.group_name,TRUE,FALSE))


    # get segments that fall within intron and annotate that segment
    altexons <- GenomicRanges::pintersect(altexons)
    if (length(altexons) == 0) {
        altexons$pos <- altexons$hit <- altexons$termini <- altexons$grouping<- NULL
        return(altexons)
    } else {
        altexons$type <- "AS"
        altexons$AStype <- toupper(altannotate$AStype)
        altexons <- altexons[!is.na(altexons$AStype) & altannotate$same.group]
        altexons$pos <- altexons$hit <- altexons$termini <-altexons$grouping <- NULL

        # add ID
        altexons.id <- altexons %>%
            as.data.frame() %>%
            dplyr::distinct(seqnames, start, end, strand, gene_id, gene_name, AStype) %>%
            dplyr::mutate(AS_id = sprintf("AS%05d", dplyr::row_number()))

        altexons <- altexons %>%
            as.data.frame() %>%
            dplyr::left_join(altexons.id, by = c("seqnames", "start", "end",
                                                 "strand", "gene_id", "gene_name",
                                                 "AStype")) %>%
            GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)


        return(BiocGenerics::sort(altexons))
    }
}

























