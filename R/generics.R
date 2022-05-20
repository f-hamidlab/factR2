#' Plot custom transcripts
#'
#' @description
#' A wrapper around wiggleplotr's plotTranscripts function.
#' See the documentation for (\code{\link[wiggleplotr]{plotTranscripts}})
#' for more information.
#'
#' @param object
#' factRObject
#'
#' @param ...
#' Character value of features to plot. Multiple features can be plotted by
#' entering comma-delimited values. Features will be extracted from
#' metadata gene_name, gene_id and transcript_id of the GTF. Can also be a
#' conditional statement to filter values from variables in the GTF (e.g. gene_name == "Ptbp1")
#'
#' @param rescale_introns
#' Specifies if the introns should be scaled to fixed length or not. (default: FALSE)
#'
#' @param ncol
#' Number of columns to patch the output plots (default: 1)
#'
#' @return ggplot2 object. If multiple genes are detected, plots will be
#' combined using patchwork
#' @export plot
#' @rdname plot
#'
#' @author Fursham Hamid
#'
#' @examples
#' ## ---------------------------------------------------------------------
#' ## EXAMPLE USING SAMPLE DATASET
#' ## ---------------------------------------------------------------------
#' # Load datasets
#' data(query_gtf, ref_gtf)
#'
#' viewTranscripts(query_gtf)
#' viewTranscripts(query_gtf, "transcript1")
#' viewTranscripts(ref_gtf)
#'
#' ## ---------------------------------------------------------------------
#' ## EXAMPLE USING TRANSCRIPT ANNOTATION
#' ## ---------------------------------------------------------------------
#' \donttest{
#' library(AnnotationHub)
#'
#' ## Retrieve GRCm38 trancript annotation
#' ah <- AnnotationHub()
#' GRCm38_gtf <- ah[["AH60127"]]
#'
#' ## Plot transcripts from Ptbp1 gene
#' viewTranscripts(GRCm38_gtf, "Ptbp1")
#'
#' # Plot transcripts from Ptbp1 and Ptbp2 genes
#' viewTranscripts(GRCm38_gtf, "Ptbp1", "Ptbp2")
#' }
#'
plot <- function(object, ..., rescale_introns = FALSE, ncol = 1) {
    UseMethod(generic = 'plot', object = object)
}
