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
ListSupportedGenomes <- function(){
    data("genomes")
    genomes[c("ID","species","database","release.date")]
}

factR2version <- "0.99.0"
