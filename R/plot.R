setMethod("Plot", "factR", function(object, ..., type = "transcripts",
                                    rescale_introns = FALSE, ncol = 1) {

    # check features
    genetxs <- methods::slot(object, name = "custom")$genetxs
    if(missing(...)){
        txs <- genetxs$transcript_id
    } else {
        genetxs.features <- genetxs %>%
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

    # select features by data
    if(type == "transcripts"){
        x <- methods::slot(object, "custom")$ranges
        x <- x[x$transcript_id %in% txs]

        genes <- unique(x$gene_name)

        # control for too many genes to plot
        if(length(genes) > 9){
            rlang::warn("Too many genes to plot, plotting first 9")
            genes <- genes[1:9]
            ncol <-  3
        }
        plot <- BiocGenerics::do.call(
            patchwork::wrap_plots,
              lapply(genes, function(y){
                  # Fetch gene exons and cdss
                  exons <- S4Vectors::split(x[x$type == "exon" & x$gene_name == y], ~transcript_id)
                  cdss <- S4Vectors::split(x[x$type == "CDS" & x$gene_name == y], ~transcript_id)
                  #as <- S4Vectors::split(x[x$type == "AS" & x$gene_name == y], ~transcript_id)
                  if (length(cdss) == 0) {
                      cdss <- NULL
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

})



