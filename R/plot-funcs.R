setMethod("plotTranscripts", "factR", function(object, ...,
                                    rescale_introns = FALSE, ncol = 1) {


    # select features by data
    feat <- .getFeat(object, ...)
    x <- object@transcriptome
    x <- x[x$transcript_id %in% feat & !x$type %in% c("AS", "gene")]

    # correct genes with no gene name
    x$gene_name <- ifelse(is.na(x$gene_name), x$gene_id, x$gene_name)
    genes <- unique(x$gene_name)

    # control for too many genes to plot
    if(length(genes) > 1){
        rlang::warn("Too many genes to plot, plotting first 9")
        genes <- genes[1]
        x <- x[x$gene_name %in% genes]
    }


    .plotTx(x)
    # plot <- BiocGenerics::do.call(
    #     patchwork::wrap_plots,
    #     lapply(genes, function(y){
    #         # Fetch gene exons and cdss
    #         exons <- S4Vectors::split(x[x$type == "exon" & x$gene_name == y], ~transcript_id)
    #         cdss <- S4Vectors::split(x[x$type == "CDS" & x$gene_name == y], ~transcript_id)
    #         #as <- S4Vectors::split(x[x$type == "AS" & x$gene_name == y], ~transcript_id)
    #         if (length(cdss) == 0) {
    #             cdss <- NULL
    #         }
    #
    #         # main plot function
    #         suppressWarnings(wiggleplotr::plotTranscripts(
    #             exons = exons,
    #             cdss = cdss[names(cdss) %in% names(exons)],
    #             rescale_introns = rescale_introns
    #         )) + ggplot2::ggtitle(y) +
    #             ggplot2::theme(strip.background.y = ggplot2::element_blank())
    #     }))
    #
    # plot + patchwork::plot_layout(ncol = ncol)

})



setMethod("plotDomains", "factR", function(object, ..., ncol = 1){
    # check if CDS have been built
    gtf <- object@transcriptome
    if(! "CDS" %in% gtf$type){
        rlang::abort("No CDSs found. Please run buildCDS() first")
    }


    # get transcripts to test
    genetxs <- features(object, ..., set = "transcript")
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
    datatoplot <- datatoplot[datatoplot$transcript_id %in% txs,]
    datatoplot <- datatoplot %>%
        dplyr::left_join(genetxs[c("transcript_id","gene_name")] ,
                         by = c("transcript_id"="transcript_id"))
    datatoplot.order <- datatoplot %>%
        dplyr::distinct(transcript_id, gene_name) %>%
        dplyr::group_by(gene_name) %>%
        dplyr::mutate(order = dplyr::row_number())
    datatoplot <- datatoplot %>% dplyr::left_join(datatoplot.order,
                                                  by = c("transcript_id", "gene_name"))

    ## return warning for non-coding transcripts
    if(any(!txs %in% datatoplot$transcript_id)){
        not.plotted <- txs[!txs %in% datatoplot$transcript_id]
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
        dplyr::distinct(order, transcript_id) %>%
        dplyr::arrange(order)
    domains.data <- data %>%
        dplyr::filter(type == "DOMAIN")
    data %>%
        dplyr::mutate(transcript_id = factor(transcript_id, levels = unique(transcript_id))) %>%
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
            labels = data.names$transcript_id) +
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









.plotTx <- function(gtf){

    gtf <- as.data.frame(gtf)
    order <- gtf %>%
        dplyr::distinct(transcript_id, gene_name) %>%
        dplyr::group_by(gene_name) %>%
        dplyr::mutate(order = dplyr::row_number())
    gtf <- gtf %>% dplyr::left_join(order,
                                    by = c("transcript_id", "gene_name"))

    data.names <- gtf %>%
        dplyr::distinct(order, transcript_id, strand, start, end) %>%
        dplyr::arrange(order)

    # plot canvas
    txs <- dplyr::filter(gtf, type == "transcript")
    range <- c(min(txs$start), (max(txs$end)))
    buffer <- 0.1*(range[2]-range[1])
    range[1] <- range[1]-buffer
    range[2] <- range[2]+buffer


    plot <- ggplot2::ggplot(gtf, ggplot2::aes(y=order)) +
        ggplot2::xlim(range) +
        labs(y = "", x = sprintf("Genome position (%s)", unique(gtf$seqnames))) +
        theme_bw() +
        scale_y_continuous(
            breaks = data.names$order,
            labels = data.names$transcript_id)

    # plot transcripts
    plot <- plot +
        geom_segment(data = txs, mapping = aes(x = start, xend = end, yend = order),
                     colour = "#0000b2", size = 0.2)

    # plot exons
    exons <- dplyr::filter(gtf, type == "exon")
    plot <- plot +
        geom_rect(data = exons,
                  fill = "#0000b2",
                  mapping = aes(xmin = start, xmax = end,
                                              ymin = order-0.06, ymax = order+0.06))


    # plot CDS
    if("CDS" %in% gtf$type){
        cds <- dplyr::filter(gtf, type == "CDS")
        plot <- plot +
            geom_rect(data = cds,
                      fill = "#0000b2",
                      mapping = aes(xmin = start, xmax = end,
                                    ymin = order-0.15, ymax = order+0.15))
    }


    # plot arrows
    increments <- (range[2]-range[1])/15
    arrows <- data.names %>%
        dplyr::mutate(seqnames = unique(gtf$seqnames)) %>%
        dplyr::mutate(range = list(seq(range[1], range[2], by =increments))) %>%
        tidyr::unnest(cols = c(range)) %>%
        dplyr::mutate(rangeend = ifelse(strand == "-", range+1, range-1)) %>%
        dplyr::filter(range > (start+100) & range < (end-100))
    # plot <- plot +
    #     geom_segment(data = arrows, mapping = aes(x=range, xend=rangeend,
    #                                               yend=order), color = "#0000b2",
    #                  arrow = arrow(length = unit(0.01, "npc")),
    #                  size = 0.5, lineend = "round", linejoin = "round")


    #return(plot)


    g <- plotly::ggplotly(plot) %>%
        plotly::add_annotations(data = arrows, text = "", x=~range,
                                y = arrows$order, ay = 0.0000001, showarrow = TRUE,
                                arrowcolor = "#0000b2", arrowsize = 0.8)

    return(g)


}


















