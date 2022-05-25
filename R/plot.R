#' @rdname Plot
#' @export Plot
Plot <- function(object, ..., rescale_introns = FALSE, ncol = 1) {

    gene_name <- gene_id <- transcript_id <- NULL
    transcript_id <- meta <- val <- n <- type <-  NULL

    # catch missing args
    mandargs <- c("object")
    passed <- names(as.list(match.call())[-1])
    if (any(!mandargs %in% passed)) {
        rlang::abort(paste(
            "missing values for",
            paste(setdiff(mandargs, passed), collapse = ", ")
        ))
    }

    # retrieve input object names
    argnames <- as.character(match.call())[-1]


    # prepare features
    x <- methods::slot(object, "custom")$ranges
    featmeta <- tryCatch(
        {
            GenomicRanges::mcols(x) %>%
                as.data.frame() %>%
                dplyr::select(gene_name, gene_id, transcript_id) %>%
                dplyr::mutate(n = dplyr::row_number()) %>%
                tidyr::gather(meta, val, -n)
        },
        error = function(e) {
            GenomicRanges::mcols(x) %>%
                as.data.frame() %>%
                dplyr::select(-type) %>%
                dplyr::mutate(n = dplyr::row_number()) %>%
                tidyr::gather(meta, val, -n)
        }
    )

    if (!missing(...)) {
        x <- tryCatch(
            {

                x[featmeta[featmeta$val %in% c(...),"n"]]
            },
            error = function(e) {
                y <- tryCatch({
                    x %>% as.data.frame() %>%
                        dplyr::filter(...) %>%
                        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
                },
                error = function(e){
                    rlang::abort(sprintf(
                        "Variables given in ... are not found in `%s`",
                        argnames[1]
                    ))
                })
            }
        )
        if (length(x) == 0) {
            rlang::abort("No transcripts to plot")
        }
    }

    # Need to have a check for plotting multiple genes.....
    ngenes <- unique(x$gene_name)
    plot <- BiocGenerics::do.call(patchwork::wrap_plots,
                                  lapply(ngenes, function(y){
                                      # Fetch gene exons and cdss
                                      exons <- S4Vectors::split(x[x$type == "exon" & x$gene_name == y], ~transcript_id)
                                      cdss <- S4Vectors::split(x[x$type == "CDS" & x$gene_name == y], ~transcript_id)
                                      as <- S4Vectors::split(x[x$type == "AS" & x$gene_name == y], ~transcript_id)
                                      if (length(cdss) == 0) {
                                          cdss <- NULL
                                      }


                                      # Control check for number of plotted transcripts
                                      if (length(exons) > 25) {
                                          exons <- exons[seq_len(25)]
                                          rlang::warn(sprintf("Plotting only first 25 transcripts for %s gene", y))
                                      }

                                      # main plot function
                                      suppressWarnings(wiggleplotr::plotTranscripts(
                                          exons = exons,
                                          cdss = cdss[names(cdss) %in% names(exons)],
                                          rescale_introns = rescale_introns
                                      )) + ggplot2::ggtitle(y)
                                  }))


    plot + patchwork::plot_layout(ncol = ncol)
}
