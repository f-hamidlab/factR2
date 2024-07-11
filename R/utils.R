is_gtf <- function(x) {
    if (is(x, "GRanges")) {
        x <- as.data.frame(x)
        if (all(c("type", "gene_id", "transcript_id") %in% names(x))) {
            x <- x %>%
                dplyr::select(type, gene_id, transcript_id) %>%
                dplyr::distinct()
            if (nrow(x) > 1) {
                return(TRUE)
            }
        }
    }
    return(FALSE)
}

is_dir <- function(x){
    tools::file_ext(x) == ""
}
is_file <- function(x){
    tools::file_ext(x) != ""
}

#TODO: filter genomes by keywords

#' List supported genomes
#'
#' @return Dataframe of the IDs of genomes and its description
#'
#'
#' @export
#'
#' @examples
#' listSupportedGenomes()
listSupportedGenomes <- function(){
    data("genomes")
    genomes[c("ID","species","database","release.date")]
}

factR2version <- "0.99.0"


.getFeat <- function(object, ..., out = "transcript_id"){
    if(missing(...)){
        return(unique(na.omit(S4Vectors::mcols(object@transcriptome)[[out]])))
    } else {
        if(out == "AS_id"){
            object@transcriptome %>%
                as.data.frame() %>%
                dplyr::filter(type == "AS") %>%
                dplyr::mutate(AS = AS_id) %>%
                dplyr::select(gene_id, gene_name, transcript_id, AS, AS_id) %>%
                tidyr::gather("type", "feature", gene_id, gene_name, transcript_id, AS) %>%
                dplyr::filter(feature %in% c(...)) %>%
                dplyr::pull(AS_id) %>%
                na.omit() %>%
                unique()
        } else {
            object[["transcript"]] %>%
                dplyr::mutate(tx = transcript_id, gene = gene_id) %>%
                dplyr::select(gene, gene_name, tx, !!out) %>%
                tidyr::gather("type", "feature", gene, gene_name, tx) %>%
                dplyr::filter(feature %in% c(...)) %>%
                dplyr::pull(!!out) %>%
                na.omit() %>%
                unique()
        }

    }
}

.sortGTF <- function(x){
  type.order <- c("gene", "transcript", "exon",
                  "start_codon", "CDS", "stop_codon", "AS")

  gene.order <- x %>%
    as.data.frame() %>%
    dplyr::filter(type %in% "gene") %>%
    dplyr::arrange(seqnames, start) %>%
    dplyr::mutate(gene_order = dplyr::row_number()) %>%
    dplyr::select(gene_id, gene_order)

  tx.order <- x %>%
    as.data.frame() %>%
    dplyr::filter(type %in% "transcript") %>%
    dplyr::group_by(gene_id) %>%
    dplyr::arrange(ifelse(strand == "-", dplyr::desc(end), start),
                   ifelse(strand == "-", start, dplyr::desc(end) )) %>%
    dplyr::mutate(tx_order = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::select(transcript_id, tx_order)

  x$gene_order <- gene.order$gene_order[match(x$gene_id, gene.order$gene_id)]
  x$tx_order <- tx.order$tx_order[match(x$transcript_id, tx.order$transcript_id)]
  x$tx_order <- ifelse(is.na(x$tx_order), 0, x$tx_order)
  x$type <- factor(as.character(x$type), type.order)

  x.order <- x %>%
    as.data.frame() %>%
    dplyr::select(gene_order, tx_order, type, strand, start) %>%
    dplyr::mutate(i = dplyr::row_number()) %>%
    dplyr::arrange(gene_order, tx_order, type,
                   ifelse(strand == "-", dplyr::desc(start), start))
  x$gene_order <- NULL
  x$tx_order <- NULL

  x[x.order$i]
}


.get_coord <- function(x){
  paste0(GenomicRanges::seqnames(x), ":",
         GenomicRanges::start(x), "-",
         GenomicRanges::end(x))
}

.get_coord_strand <- function(x, sep = ";"){
  paste0(GenomicRanges::seqnames(x), ":",
         GenomicRanges::start(x), "-",
         GenomicRanges::end(x),
         sep, GenomicRanges::strand(x))
}

.get_coord_geneid_strand <- function(x, sep = ";"){
  paste0(GenomicRanges::seqnames(x), ":",
         GenomicRanges::start(x), "-",
         GenomicRanges::end(x),
         sep, x$gene_id,
         sep, GenomicRanges::strand(x))
}
.asinTransform <-  function(x) { 2 * asin(sqrt(x))/pi }

.msgheader <- function(x){
    message(green$bold(stringr::str_glue("\U1F846  {x}")))
}
.msgsubinfo <- function(x){
    message(blue("    \U2139 "), italic$blue(x))
}
.msgsubwarn <- function(x){
    message(red("    \U2757 "), italic$cyan(x))
}
.msgwarn <- function(x){
    message(red("\U2757 "), italic$cyan(x))
}
.msginfo <- function(x){
    message(blue("\U2139 "), italic$blue(x))
}


