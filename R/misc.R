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


.getTxs <- function(object, ...){
    if(missing(...)){
        txs <- unique(obj@transcriptome$transcript_id)
    } else {
        genetxs.features <- obj@transcriptome %>%
            as.data.frame() %>%
            dplyr::mutate(tx = transcript_id) %>%
            tidyr::gather("type", "feature", gene_id, gene_name, transcript_id) %>%
            dplyr::filter(feature %in% c(...))
        if (nrow(genetxs.features) == 0) {
            rlang::abort("No features found in object")
        } else if(!all(c(...) %in% genetxs.features$feature)){
            absent.features <- c(...)[which(!c(...) %in% genetxs.features$feature)]
            rlang::warn(sprintf("These features are not found in object: %s",
                                paste(absent.features, collapse = ", ")))
        }
        txs <- genetxs.features$tx
    }
    txs
}
