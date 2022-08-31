setMethod("plotTranscripts", "factR", function(object, ...,
                                               collapse = FALSE,
                                               rescale_introns = FALSE,
                                               ncol = 1) {



    x <- object@transcriptome
    # handle chromosome inputs
    if(stringr::str_detect(...,":|-")){
        exon <- GenomicRanges::GRanges(...)
        xtxs <- x[x$type == "exon"]
        hits <- IRanges::findOverlaps(exon, xtxs)
        genes <- unique(xtxs[S4Vectors::subjectHits(hits)]$gene_name)
        x <- x[x$gene_name %in% genes & !x$type %in% c("AS", "gene")]
        xrange <-  stringr::str_split(..., ":|-")[[1]][c(2,3)]
    } else if(stringr::str_detect(..., "^AS[0-9]{5}")) {
        exon <- ase(object)[...,]
        x <- x[x$gene_id %in% exon$gene_id & !x$type %in% c("AS", "gene")]
        xrange <-  stringr::str_split(exon$coord, ":|-")[[1]][c(2,3)]

    }else {
        # select features by data
        feat <- .getFeat(object, ...)
        x <- x[x$transcript_id %in% feat & !x$type %in% c("AS", "gene")]
        xrange <- NULL

    }


    # correct genes with no gene name
    x$gene_name <- ifelse(is.na(x$gene_name), x$gene_id, x$gene_name)
    genes <- unique(x$gene_name)

    # control for too many genes to plot
    if(length(genes) > 1){
        rlang::warn("Too many genes to plot, plotting first 9")
        genes <- genes[1]
        x <- x[x$gene_name %in% genes]
    }

    .plotTx(x, collapse, rescale_introns, xrange)
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

# TODO: color reference transcript

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









.plotTx <- function(gtf, collapse = FALSE, rescale = FALSE, xrange = NULL){


    if(collapse){
        # make collapsed exon arch
        new.gtf <- IRanges::reduce(gtf[gtf$type == "exon"])
        new.gtf$transcript_id <- unique(gtf$gene_name)
        new.gtf$gene_name <- unique(gtf$gene_name)
        new.gtf$type <- "exon"

        # add transcript biotype
        tx.gtf <- range(new.gtf)
        tx.gtf$type <- "transcript"
        tx.gtf$transcript_id <- unique(gtf$gene_name)
        tx.gtf$gene_name <- unique(gtf$gene_name)
        new.gtf <- c(new.gtf, tx.gtf)

        # add CDS if present
        if("CDS" %in% gtf$type){
            cds.gtf <- IRanges::reduce(gtf[gtf$type == "CDS"])
            cds.gtf$transcript_id <- unique(gtf$gene_name)
            cds.gtf$gene_name <- unique(gtf$gene_name)
            cds.gtf$type <- "CDS"
            new.gtf <- c(new.gtf, cds.gtf)

        }


        gtf <- new.gtf
    }


    if(rescale){
        new.gtf <- IRanges::reduce(gtf[gtf$type == "exon"])
        new.gtf <- new.gtf %>%
            as.data.frame() %>%
            dplyr::arrange(start, end) %>%
            dplyr::mutate(space = width+50) %>%
            dplyr::mutate(new.start = dplyr::lag(cumsum(space))) %>%
            tidyr::replace_na(list(new.start = 0)) %>%
            dplyr::mutate(new.end = new.start +width) %>%
            GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)

        exons <- gtf[gtf$type != "transcript"]
        new.gtf <- IRanges::mergeByOverlaps(exons, new.gtf) %>%
           as.data.frame() %>%
            dplyr::filter(exons.start >= new.gtf.start,
                          exons.end <= new.gtf.end) %>%
            dplyr::mutate(new.start = new.start + (exons.start-new.gtf.start),
                          new.end = new.end - (new.gtf.end-exons.end)) %>%
            dplyr::select(seqnames = exons.seqnames, start = new.start, end = new.end,
                          strand = exons.strand,
                          transcript_id = exons.transcript_id, gene_name = exons.gene_name,
                          type = exons.type) %>%
            GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)

        # add transcript info
        tx.gtf <- range(S4Vectors::split(new.gtf, ~transcript_id), ignore.strand = F)
        tx.gtf <- as.data.frame(tx.gtf) %>%
            dplyr::select(seqnames, start, end, strand, transcript_id = group_name) %>%
            dplyr::mutate(type = "transcript") %>%
            dplyr::left_join(dplyr::distinct(as.data.frame(gtf)[c("transcript_id", "gene_name")]),
                             by = "transcript_id") %>%
            GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
        gtf <- c(tx.gtf, new.gtf)
    }

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
    buffer <- 0.0005*(range[2]-range[1])
    range[1] <- range[1]-buffer
    range[2] <- range[2]+buffer



    xaxis.prefix <- ifelse(rescale, "Relative position (%s)", "Genome postion (%s)")
    plot <- ggplot2::ggplot(gtf, ggplot2::aes(y=order)) +
        ggplot2::xlim(range) +
        labs(y = "", x = sprintf(xaxis.prefix, unique(gtf$seqnames))) +
        theme_bw() +
        scale_y_continuous(
            limits = c(0.6, max(data.names$order)+0.4),
            breaks = data.names$order,
            labels = data.names$transcript_id)

    # plot transcripts
    plot <- plot +
        geom_segment(data = txs, mapping = aes(x = start, xend = end, yend = order),
                     colour = "#0000b2", size = 0.2)

    # plot exons
    ntranscripts <- max(data.names$order)
    exon.height <- 0.025 + (ntranscripts^2 * 0.0005)
    exons <- dplyr::filter(gtf, type == "exon")
    plot <- plot +
        geom_rect(data = exons,
                  fill = "#0000b2",
                  mapping = aes(xmin = start, xmax = end,
                                ymin = order-exon.height,
                                ymax = order+exon.height))


    # plot CDS
    cds.height <- exon.height *2
    if("CDS" %in% gtf$type){
        cds <- dplyr::filter(gtf, type == "CDS")
        plot <- plot +
            geom_rect(data = cds,
                      fill = "#0000b2",
                      mapping = aes(xmin = start, xmax = end,
                                    ymin = order-cds.height, ymax = order+cds.height))
    }


    # plot arrows
    increments <- (range[2]-range[1])/15
    arrows <- data.names %>%
        dplyr::mutate(seqnames = unique(gtf$seqnames)) %>%
        dplyr::mutate(range = list(seq(range[1], range[2], by =increments))) %>%
        tidyr::unnest(cols = c(range)) %>%
        dplyr::mutate(dir = ifelse(strand == "-", "5", "-5")) %>%
        dplyr::filter(range > (start+100) & range < (end-100))


    # plot <- plot +
    #     geom_segment(data = arrows, mapping = aes(x=range, xend=rangeend,
    #                                               yend=order), color = "#0000b2",
    #                  arrow = arrow(length = unit(0.01, "npc")),
    #                  size = 0.5, lineend = "round", linejoin = "round")

    # TODO: annotate segments (if possible)

    #return(plot)
    if(nrow(arrows)>1){
        g <- plotly::ggplotly(plot) %>%
            plotly::add_annotations(data = arrows, text = "", x=~range,
                                    axref = "pixel", ax=~dir,
                                    y = arrows$order, ay = 0.0000001, showarrow = TRUE,
                                    arrowcolor = "#0000b2", arrowsize = 0.8)
    } else {
        g <- plotly::ggplotly(plot)
    }

    if(!is.null(xrange)){
        #g <- g %>% plotly::layout(xaxis = list(range = list(xrange[1], xrange[2])))
        #g <- g %>% plotly::layout(xaxis = list(rangeslider = list(start = 79695000, end = 79696000)))
        g <- g %>%
            plotly::rangeslider(start = xrange[1], xrange[2])
    }




    return(g)


}


















