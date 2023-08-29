#' Plot exon and domain architectures
#'
#' @param object factR object class
#' @param x Feature to plot. Can be the following:
#' \itemize{
#'  \item{gene_id: }{ID of gene to plot}
#'  \item{gene_name: }{Name of gene to plot}
#'  \item{transcript_id: }{ID of transcript to plot}
#'  \item{AS_id: }{ID of alternative splicing event (ASxxxxx)}
#'  \item{sample: }{Column name from sample metadata}
#' }
#' @param y Same as x
#' @param plot_trend Whether to plot trendline (Default: FALSE)
#'
#' @return
#' ggplot2 plot
#' @export
#'
#' @name factR-plottings
#' @rdname factR-plottings
#'
setMethod("plot2way", "factR", function(object, x, y,
                                        plot_trend=FALSE) {
    .plot2way(object, x, y , plot_trend)

})





.plot2way <- function(object, x, y, plot_trend = FALSE){

    # create db of all features and variables
    ## feature x feature; datapoint are samples
    ## feature x samples metadata

    feature.db <- data.frame(feat = object@sets$gene@rowData$gene_id,
                             id = object@sets$gene@rowData$gene_id,
                             set = "gene",
                             label = "Normalized expression",
                             label2 = "gene expression") %>%
        dplyr::bind_rows(data.frame(feat = object@sets$gene@rowData$gene_name,
                                    id = object@sets$gene@rowData$gene_id,
                                    set = "gene",
                                    label = "Normalized expression",
                                    label2 = "gene expression")) %>%
        dplyr::bind_rows(data.frame(feat = object@sets$transcript@rowData$transcript_id,
                                    id = object@sets$transcript@rowData$transcript_id,
                                    set = "transcript",
                                    label = "Normalized expression",
                                    label2 = "transcript expression"))%>%
        dplyr::bind_rows(data.frame(feat = object@sets$AS@rowData$AS_id,
                                    id = object@sets$AS@rowData$AS_id,
                                    set = "AS",
                                    label = "Percent-Spliced In (PSI)",
                                    label2 = "")) %>%
        dplyr::bind_rows(data.frame(feat = colnames(object@colData),
                                    id = colnames(object@colData),
                                    set = "sample",
                                    label = colnames(object@colData),
                                    label2 = "groupings"))

    x.feat <- feature.db[feature.db$feat %in% x,]
    y.feat <- feature.db[feature.db$feat %in% y,]

    dat <- do.call(cbind,lapply(list(x.feat, y.feat), function(dat){
         if(dat$set == "sample"){
             as.data.frame(object@colData[[dat$feat]])
         } else {
             as.data.frame(object@sets[[dat$set]]@data[dat$id,])
         }
    }))
    colnames(dat) <- c("x.dat","y.dat")



    p <- dat %>%
        ggplot2::ggplot(aes(x=x.dat,y=y.dat)) +
        geom_point(size=3) +
        theme_bw() +
        labs(x=x.feat$label[1],y=y.feat$label[1]) +
        ggtitle(stringr::str_glue(
            "{y} {y.feat$label2[1]} vs {x} {x.feat$label2[1]}"))

    if(plot_trend){
        p <- suppressMessages(p +
            ggplot2::geom_smooth(formula = y~x,
                method = "lm", se = FALSE, color = "cornflowerblue"))
    }

    p

}
