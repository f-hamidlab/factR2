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



listSupportedGenomes <- function(){
    data("genomes")
    genomes[c("ID","species","database","release.date")]
}

factR2version <- "0.99.0"


.getFeat <- function(object, ..., out = "transcript_id"){
    if(missing(...)){
        return(unique(na.omit(S4Vectors::mcols(object@transcriptome)[[out]])))
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

.asinTransform <-  function(x) { 2 * asin(sqrt(x))/pi }

.msgheader <- function(x){
    message(green$bold(stringr::str_glue("\U1F846  {x}")))
}
.msgsubinfo <- function(x){
    message(blue("    \U2139 "), italic$white(x))
}
.msgsubwarn <- function(x){
    message(red("    \U2757 "), italic$cyan(x))
}



# TODO: reattaching genome?
