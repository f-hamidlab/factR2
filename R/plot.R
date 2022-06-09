setMethod("plotTranscripts", "factR", function(object, ...,
                                    rescale_introns = FALSE, ncol = 1) {


    # select features by data
    x <- ranges(object,..., set = "transcript") 

    # correct genes with no gene name
    x$gene_name <- ifelse(is.na(x$gene_name), x$gene_id, x$gene_name)
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
            )) + ggplot2::ggtitle(y) +
                ggplot2::theme(strip.background.y = ggplot2::element_blank())
        }))

    plot + patchwork::plot_layout(ncol = ncol)

})



setMethod("plotDomains", "factR", function(object, ..., ncol = 1){
    # check if CDS have been built
    gtf <- object@transcriptome
    if(! "CDS" %in% gtf$type){
        rlang::abort("No CDSs found. Please run buildCDS() first")
    }


    # get transcripts to test
    genetxs <- featureData(object, ..., set = "transcript")
    genetxs$gene_name <- ifelse(is.na(genetxs$gene_name), genetxs$gene_id, genetxs$gene_name)
    txs <- genetxs$transcript_id


    # check if all transcripts have been tested
    untested.txs <- txs[!txs %in% slot(object, "domains")$tested]
    if(length(untested.txs) > 0){
        object <- predictDomain(object, txs)

        argnames <- as.character(match.call())[-1]
        assign(argnames[1], object, envir = .GlobalEnv)
    }

    datatoplot <- object@domains$data
    datatoplot <- datatoplot[datatoplot$entryName %in% txs,]
    datatoplot <- datatoplot %>%
        dplyr::left_join(genetxs[c("transcript_id","gene_name")] ,
                         by = c("entryName"="transcript_id"))
    datatoplot.order <- datatoplot %>%
        dplyr::distinct(entryName, gene_name) %>%
        dplyr::group_by(gene_name) %>%
        dplyr::mutate(order = dplyr::row_number())
    datatoplot <- datatoplot %>% dplyr::left_join(datatoplot.order,
                                                  by = c("entryName", "gene_name"))

    ## return warning for non-coding transcripts
    if(any(!txs %in% datatoplot$entryName)){
        not.plotted <- txs[!txs %in% datatoplot$entryName]
        rlang::warn(sprintf("These transcripts are non-coding: %s",
                            paste(not.plotted, collapse = ", ")))
    }


    genes <- unique(datatoplot$gene_name)

    # control for too many genes to plot
    if(length(genes) > 9){
        rlang::warn("Too many genes to plot, plotting first 9")
        genes <- genes[1:9]
        ncol <-  3
    }
    plot <- BiocGenerics::do.call(
        patchwork::wrap_plots,
        lapply(genes, function(y){
            dat <- datatoplot[datatoplot$gene_name %in% y,]
            .draw(dat) +
            ggplot2::ggtitle(y)
        }))
    plot + patchwork::plot_layout(ncol = ncol)


})





.draw <- function(data){
    # data <- data %>%
    #     dplyr::mutate(order = factor(order, levels = unique(order)))
    data.names <- data %>%
        dplyr::distinct(order, entryName) %>%
        dplyr::arrange(order)
    domains.data <- data %>%
        dplyr::filter(type == "DOMAIN")
    data %>%
        dplyr::mutate(entryName = factor(entryName, levels = unique(entryName))) %>%
        ggplot2::ggplot(ggplot2::aes(y = order)) +
        ggplot2::geom_segment(ggplot2::aes(x=0, xend = end, yend=order), colour = "grey") +
        ggplot2::geom_rect(data = domains.data, ggplot2::aes(ymin = order-0.25,
                                                             ymax = order+0.25,
                                                             xmin =begin,
                                                             xmax = end,
                                                             fill = description)) +
        # ggplot2::geom_text(data = domains.data, ggplot2::aes(x =(begin + (end-begin)/2),
        #                                               y = order,
        #                                               label = description),
        #                    size = 3) +
        ggplot2::scale_y_continuous(
            breaks = data.names$order,
            labels = data.names$entryName) +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(c(0.01,0.05))) +
        ggplot2::labs(x="Amino acid number", y = "", fill = "Domains") +
        ggplot2::theme_classic() +
        ggplot2::theme(
            axis.ticks = ggplot2::element_blank(),
            axis.line = ggplot2::element_blank(),
            plot.background = ggplot2::element_rect(colour = "grey", fill=NA, size=0.5)
        ) +
        ggplot2::scale_fill_brewer(palette = "Set2")
}


























