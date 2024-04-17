.findAS2 <-  function(object) {
  gtf <- slot(object, "transcriptome")
  splice.out <- .runAS2(gtf[gtf$type == "exon"])
  object@sets$AS@rowData <- .prepColData(splice.out)
  slot(object, "transcriptome") <- c(gtf, .prepASgtf(splice.out))

  return(object)

}

.runAS2 <- function(gtf){

  ## Prefilter for genes with at least 2 multiexonic transcripts
  filtered_genes <- gtf %>%
    as.data.frame() %>%
    dplyr::group_by(gene_id, transcript_id) %>%
    dplyr::tally() %>%
    dplyr::filter(n > 2) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::tally() %>%
    dplyr::filter(n > 1)
  exons <- gtf[gtf$gene_id %in% filtered_genes$gene_id]

  ## Trim ends of transcripts to avoid TS and TE
  exons <- .trim_ends_by_gene(exons)

  ## Classify exons by position (First, Internal, Last)
  exons <- .label_exon_class(exons)

  ## Create a GenomicRanges object of all non-redundant introns
  exonsbytx <- S4Vectors::split(exons, ~transcript_id)
  intronsbytx <- GenomicRanges::psetdiff(BiocGenerics::unlist(range(exonsbytx)), exonsbytx)
  introns.nr <- unique(unlist(intronsbytx))
  names(introns.nr) <- NULL

  ## Create a disjointed version of all exons in each gene family
  disjoint.exons <- .disjoin_by_gene(exons)

  ## Pair up all disjointed exons with spliced and skipped intron junctions
  exon.juncs <- .get_juncs(disjoint.exons, introns.nr)

  ## Get junctions for Retained introns specifically
  retained.introns <- .find_retained_intron(disjoint.exons, introns.nr)


  ## Classify non-RI events and merge
  exon.juncs <- .classify_events(exon.juncs)
  full.exon.juncs <- rbind(exon.juncs, retained.introns)
  full.exon.juncs$exon_pos <- NULL


  exon.meta <- full.exon.juncs %>%
    dplyr::arrange(exon_coord) %>%
    dplyr::select(exon_coord, gene_id, gene_name, strand, transcript_ids, AStype) %>%
    dplyr::distinct() %>%
    dplyr::mutate(AS_id = sprintf("AS%05d", dplyr::row_number()), .before = exon_coord)



  # exon.junction.pairs <- full.exon.juncs %>%
  #   dplyr::mutate(exon_id = paste0(exon_coord,"_",gene_id,"_",gene_name)) %>%
  #   dplyr::select(exon_id, junc_coord, junc_type)

  # output <- list(exon.meta, exon.junction.pairs)
  # names(output) <- c("meta", "pairs")


  return(exon.meta)

}


.trim_ends_by_gene <- function(x){
  # trim-off TS and TE
  x %>%
    as.data.frame() %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::arrange(start) %>%
    dplyr::mutate(pos = dplyr::row_number()) %>%
    dplyr::mutate(pos = dplyr::case_when(pos == 1 ~ "first",
                                         pos == dplyr::n() ~ "last",
                                         .default = "a.internal")) %>%
    dplyr::group_by(seqnames, end, gene_id) %>%
    dplyr::arrange(pos, start) %>%
    dplyr::mutate(start = ifelse(pos == "first", start[1], start)) %>%
    dplyr::group_by(seqnames, start, gene_id) %>%
    dplyr::arrange(dplyr::desc(pos), dplyr::desc(end)) %>%
    dplyr::mutate(end = ifelse(pos == "last", end[dplyr::n()], end)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-pos) %>%
    dplyr::arrange(seqnames, gene_id, transcript_id, start) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)

}

#group by transcript id and label first and last exons
.label_exon_class <- function(x){
  y <- x %>%
    as.data.frame() %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::mutate(exon_pos = dplyr::case_when(
      start==min(start) ~ "first",
      end==max(end)  ~ "last",
      # start==min(start) & strand == "-" ~ "last",
      # end==max(end) & strand == "-" ~ "first",
      .default = "internal"
    ))
  x$exon_pos <- y$exon_pos
  return(x)
}

.disjoin_by_gene <- function(x){
  # actual disjoin function
  y <- GenomicRanges::disjoin(S4Vectors::split(x, ~gene_id))

  # create gene_id metacolumn
  y <- y %>%
    as.data.frame() %>%
    dplyr::mutate(gene_id = group_name) %>%
    dplyr::select(seqnames:gene_id) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  # link up disjointed exons to original exons
  y$order <- 1:length(y)
  out <- IRanges::findOverlapPairs(y, x, type="within")
  out <- out[S4Vectors::first(out)$gene_id == S4Vectors::second(out)$gene_id ]

  # clean output by labeling exon positions, nesting transcript ids and
  # get gene name
  out <- out %>%
    as.data.frame() %>%
    dplyr::select(order = first.order,
                  gene_name = second.gene_name,
                  transcript_id = second.transcript_id,
                  exon_pos = second.exon_pos) %>%
    dplyr::mutate(exon_pos = factor(exon_pos, c("internal", "first","last"))) %>%
    dplyr::group_by(order) %>%
    dplyr::arrange(exon_pos) %>%
    dplyr::summarise(gene_name = gene_name[1],
                     transcript_id = paste0(transcript_id, collapse = ";"),
                     exon_pos = exon_pos[1])
  y$gene_name <- out$gene_name
  y$transcript_ids <- out$transcript_id
  y$exon_pos <- as.character(out$exon_pos)
  y$order <- NULL

  return(y)
}

# wrapper function to run subfunctions to get spliced junc and skipped junc
.get_juncs <- function(x, y){
  exon.spljunc <- .get_spljunc(x,y)
  out <- .add_skipjunc(exon.spljunc, y)

  return(out)
}


.get_spljunc <- function(x, y){
  # make temp id of disjointed exons
  GenomicRanges::mcols(x)$index <- 1:length(x)

  # get adjacent introns for each disjointed exon
  overlap <- IRanges::findOverlapPairs(x, y, maxgap = 0L)
  adjacent <- subset(overlap,
                     IRanges::start(first) == IRanges::end(second) + 1L |
                       IRanges::end(first) == IRanges::start(second) - 1L)

  # label type of junction
  pair_df <- as.data.frame(adjacent) %>%
    dplyr::mutate(position = ifelse(first.X.start == second.end + 1, "Upstream", "Downstream"))
  GenomicRanges::mcols(adjacent)$position <- pair_df$position

  return(adjacent)
}

.add_skipjunc <- function(x, y){

  # get introns that are covering entire exon
  overlap <- IRanges::findOverlapPairs(S4Vectors::first(x), y, type = "within")

  # extract pairs object
  x.overlap <- S4Vectors::first(x)
  y.overlap <- S4Vectors::second(overlap)

  # get exons and its skipped junction pair
  exon.w.skipjunc <- as.data.frame(overlap) %>%
    dplyr::mutate(exon_coord = .get_coord(S4Vectors::first(overlap)),
                  junc_coord = .get_coord(S4Vectors::second(overlap)),
                  junc_type = "Skipped") %>%
    dplyr::distinct(first.index,junc_coord, .keep_all = TRUE) %>%
    dplyr::select(exon_coord, junc_coord,
                  gene_id = first.gene_id,
                  gene_name = first.gene_name,
                  transcript_ids = first.transcript_ids,
                  strand = first.X.strand,
                  junc_type,
                  exon_pos=first.exon_pos,
                  first.index)

  # get exons and its spliced junction pair
  exon.w.spljunc <- as.data.frame(x) %>%
    dplyr::mutate(exon_coord = .get_coord(S4Vectors::first(x)),
                  junc_coord = .get_coord(S4Vectors::second(x))) %>%
    dplyr::filter(first.index %in% exon.w.skipjunc$first.index) %>%
    dplyr::distinct(first.index,junc_coord, .keep_all = TRUE) %>%
    dplyr::select(exon_coord,junc_coord,
                  gene_id = first.gene_id,
                  gene_name = first.gene_name,
                  transcript_ids = first.transcript_ids,
                  strand=first.X.strand,
                  junc_type=position,
                  exon_pos=first.exon_pos)

  exon.w.skipjunc$first.index <- NULL

  return(rbind(exon.w.spljunc, exon.w.skipjunc))

}


.find_retained_intron <- function(x, y){
  retained_intron <- IRanges::findOverlapPairs(x, y, type = "equal")
  retained_intron <- retained_intron[GenomicRanges::width(S4Vectors::first(retained_intron)) > 50]
  # get "skipped" junction of RI and prep output
  ri_skipped <- as.data.frame(retained_intron) %>%
    dplyr::mutate(exon_coord = .get_coord(S4Vectors::first(retained_intron)),
                  junc_coord = exon_coord,
                  junc_type = "Skipped") %>%
    dplyr::select(exon_coord,junc_coord,
                  gene_id = first.gene_id,
                  gene_name = first.gene_name,
                  transcript_ids = first.transcript_ids,
                  strand=first.X.strand,
                  junc_type,
                  exon_pos=first.exon_pos)

  # get first and last coord of RI
  ri_resized_up <- GenomicRanges::resize(S4Vectors::first(retained_intron), 1)
  ri_resized_dn <- GenomicRanges::resize(S4Vectors::first(retained_intron), 1,
                                         fix = "end")

  # prepare output of "spliced" coordinates of RI
  ri_spliced <- as.data.frame(retained_intron) %>%
    dplyr::mutate(exon_coord = .get_coord(S4Vectors::first(retained_intron)),
                  Upstream = .get_coord(ri_resized_up),
                  Downstream = .get_coord(ri_resized_dn)) %>%
    dplyr::select(exon_coord,Upstream, Downstream,
                  gene_id = first.gene_id,
                  gene_name = first.gene_name,
                  transcript_ids = first.transcript_ids,
                  strand=first.X.strand,
                  exon_pos=first.exon_pos) %>%
    tidyr::pivot_longer(Upstream:Downstream, names_to = "junc_type",
                        values_to = "junc_coord")


  ri_df <- rbind(ri_skipped,ri_spliced)
  if(nrow(ri_df) > 0){
    ri_df$AStype <- "RI"
  }


  return(ri_df)
}

.classify_events <- function(x,y){

  # test if edges of exon and introns are exact
  x$common.edge <- .testedges(x$exon_coord, x$junc_coord)

  # get distinct exon_coord and junc_type and pivot junc_type
  x.pivoted <- x %>%
    dplyr::group_by(exon_coord, gene_id) %>%
    dplyr::mutate(common.edge = any(common.edge)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-junc_coord) %>%
    dplyr::distinct() %>%
    dplyr::mutate(score = TRUE) %>%
    tidyr::pivot_wider(names_from = junc_type, values_from = score,
                       values_fill = FALSE)

  # classify events based on exon_pos and presence of Upstream and Downstream junc
  # the last bit handles AD and AA events on the negative strand
  x.classified <- x.pivoted %>%
    dplyr::mutate(AStype = dplyr::case_when(
      Skipped & Downstream & Upstream ~"CE",
      Skipped & Downstream & exon_pos=="first" & !common.edge ~ "Af",
      Skipped & Upstream & exon_pos=="last" & !common.edge~ "Al",
      Skipped & Downstream & common.edge ~ "Ad",
      Skipped & Upstream & common.edge ~ "Aa"
    )) %>%
    dplyr::mutate(AStype = ifelse(strand == "-", chartr("adfl", "dalf", AStype), AStype)) %>%
    dplyr::mutate(AStype = toupper(AStype)) %>%
    dplyr::select(exon_coord, gene_id, AStype)

  x %>%
    dplyr::select(-common.edge) %>%
    dplyr::left_join(x.classified, by = c("exon_coord","gene_id"))

}

.testedges <- function(x, y){
  x <- as(x, "GRanges")
  y <- as(y, "GRanges")
  GenomicRanges::start(x) == GenomicRanges::start(y) |
    GenomicRanges::end(x) == GenomicRanges::end(y)
}

.prepColData <- function(x){
  x %>%

    dplyr::mutate(AS_id2 = AS_id) %>%
    dplyr::mutate(width = GenomicRanges::width(as(x$exon_coord, "GRanges"))) %>%
    dplyr::select(AS_id, AS_id2, gene_id, gene_name, coord = exon_coord,
                  AStype, strand, width) %>%
    dplyr::mutate(AStype = factor(AStype, levels = c("CE", "AD","AA","AF","AL","RI"))) %>%
    tibble::column_to_rownames("AS_id2")
}

.prepASgtf <- function(x){
  x %>%
    tidyr::separate(exon_coord, c("seqnames", "start", "end")) %>%
    dplyr::mutate(transcript_id = stringr::str_split(transcript_ids, ";")) %>%
    dplyr::mutate(type="AS") %>%
    dplyr::select(seqnames, start, end, gene_id, type, gene_name, transcript_id,
                  AStype, AS_id, strand) %>%
    dplyr::mutate(AStype = factor(AStype, levels = c("CE", "AD","AA","AF","AL","RI"))) %>%
    tidyr::unnest(c("transcript_id")) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
}
