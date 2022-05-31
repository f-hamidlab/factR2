.addExonstoTx <- function(x, y, drop.unmodified = FALSE,
                         allow.external.exons = FALSE) {
    # check for length of x and y
    if (length(x) != length(y)){
        rlang::abort("x and y are not of the same length")
    }

    # drop all metadata first
    S4Vectors::mcols(x, level = "within") <- NULL

    # retain pairs with same chr and strand
    chr.x <- as.character(S4Vectors::runValue(GenomeInfoDb::seqnames(x)))
    strand.x <- as.character(S4Vectors::runValue(BiocGenerics::strand(x)))
    chr.y <- as.character(GenomeInfoDb::seqnames(y))
    strand.y <- as.character(BiocGenerics::strand(y))

    samechr.strand <-  chr.x == chr.y & strand.x == strand.y
    if(sum(samechr.strand) == 0) {
        rlang::abort("All x-y pairs are on different chromosome and/or strand")
    } else if(sum(samechr.strand) < length(x) & drop.unmodified) {
        rlang::warn(sprintf("%s x-y pairs are on different chromosome and/or strand and
  these transcripts were not modified",
                            length(x)-sum(samechr.strand)))
    }

    # check if exons in y are within transcripts in x
    internal <- as.data.frame(IRanges::pintersect(range(x), y))$hit
    if(allow.external.exons) {
        external.addition <- sum(samechr.strand & !internal)
        if(external.addition > 0) {
            rlang::warn(sprintf("%s external exons were added", external.addition))
        }
    } else {
        samechr.strand <- samechr.strand & internal
    }

    if(drop.unmodified){
        return(GenomicRanges::reduce(IRanges::punion(x,y))[samechr.strand])
    } else {
        modified <- GenomicRanges::reduce(IRanges::punion(x,y))[samechr.strand] %>%
            factR::mutateeach(modified = TRUE)
        unmodified <- x[!samechr.strand] %>% factR::mutateeach(modified = FALSE)

        return(c(modified, unmodified)[names(x)])
    }

}



.removeExonsfromTx <- function(x, y, drop.unmodified = FALSE,
                              ignore.strand = FALSE) {

    # checks for length and overlaps. returns nonoverlaps granges
    nonoverlaps.x <- .checklengthandoverlap(x, y , drop.unmodified,
                                            ignore.strand)

    # drop all metadata first
    S4Vectors::mcols(x, level = "within") <- NULL

    # merge x and y and prepare for intron removals
    x.y <- .mergexyandcorrectIR(x, y)
    # x.y$modified <- ifelse(x.y$group_name %in% names(nonoverlaps.x),
    #                        FALSE, TRUE)

    # Re-generate x and y GRanges
    x <- x.y %>%
        dplyr::select(group, group_name, ends_with("x")) %>%
        dplyr::rename_with(~ stringr::str_sub(.x,end=-3), ends_with("x")) %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)

    y <- x.y %>%
        dplyr::select(group, ends_with("y")) %>%
        dplyr::rename_with(~ stringr::str_sub(.x,end=-3), ends_with("y")) %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)


    # Core part: Removal of exon segments from transcripts in x
    diff <- suppressWarnings(
        GenomicRanges::psetdiff(x,y, ignore.strand = ignore.strand))
    diff$group_name <- x$group_name

    # prepare output grangeslist
    diff <- diff %>% as.data.frame() %>%
        dplyr::mutate(group_name = x$group_name) %>%
        dplyr::filter(width > 0) %>%
        GenomicRanges::makeGRangesListFromDataFrame(split.field = "group_name",
                                                    keep.extra.columns=TRUE) %>%
        factR::sorteach(exonorder) %>%
        GenomicRanges::reduce() %>%
        factR::mutateeach(modified = ifelse(group_name %in% names(nonoverlaps.x),
                                     FALSE, TRUE))

    # drop non-overlapping pairs if requested, or annotate modified tx
    if(drop.unmodified) {
        diff <- factR::filtereach(diff, modified) %>%
            factR::mutateeach(modified = NULL)
    }

    return(diff)
}

.checklengthandoverlap <- function(x, y, drop, ignore.strand) {
    # check for length of x and y
    if (length(x) != length(y)){
        rlang::abort("x and y are not of the same length")
    }

    # determine poverlaps between x and y and report
    temp.x <- IRanges::pintersect(x, y, ignore.strand = ignore.strand) %>%
        factR::filtereach(sum(hit) == 0)

    if(length(temp.x) == length(x)) {
        rlang::abort("All x-y pairs do not overlap")
    } else if(length(temp.x) > 0){
        if(drop) {
            rlang::warn(sprintf("%s x-y pairs do not overlap and these transcripts were not modified",
                                length(temp.x)))
        }
    }
    return(temp.x)
}

.mergexyandcorrectIR <- function(x,y) {

    # combine x and y by its group, also annotate which pair could be an IR
    temp.y <- as.data.frame(y) %>%
        dplyr::mutate(group=dplyr::row_number())
    x.y <- as.data.frame(x) %>%
        dplyr::left_join(temp.y, by = "group") %>%
        dplyr::mutate(is.internal = ifelse(start.x < start.y & end.x > end.y,
                                           TRUE, FALSE))

    # extract all IR events and change the end coord
    ## this will serve as a duplicate
    ir2 <- x.y[x.y$is.internal,] %>%
        dplyr::mutate(end.y = end.x)

    # change start coord of IR events and combine with ir2
    ## also, filter out entries of different chr/strand
    x.y <- x.y %>%
        dplyr::mutate(start.y = ifelse(is.internal, start.x, start.y)) %>%
        dplyr::bind_rows(ir2)
    return(x.y)
}
