is_gtf <- function(...) {
    type <- gene_id <- transcript_id <- NULL
    return(unlist(lapply(list(...), function(x) {
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
    })))
}

#' List factR supported genomes
#'
#' @return data.frame with a list of factR-supported genomes
#' @export
#'
#' @examples
#' ListGenomes()
listSupportedGenomes <- function(){
    data("genomes")
    genomes[c("ID","species","database","release.date")]
}

factR2version <- "0.99.0"


.getFeat <- function(object, ..., out = "transcript_id"){
    if(missing(...)){
        return(unique(S4Vectors::mcols(object@transcriptome)[[out]]))
    } else {
        genetxs.features <- object@transcriptome %>%
            as.data.frame() %>%
            dplyr::filter(!type %in% "gene") %>%
            dplyr::mutate(tx = transcript_id, gene = gene_id) %>%
            tidyr::gather("type", "feature", gene, gene_name, tx) %>%
            dplyr::filter(feature %in% c(...))
        return(unique(genetxs.features[[out]]))
        
        
    }
}
