#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Create factRObject
#'
#' @description
#' factRObject stores all imported inputs, intermediate objects and results of
#' factR workflow. To construct a factRObject, you need a custom transcriptome
#' in GTF format, a reference annotation and a genome sequence. See examples
#' [NOT COMPLETE]
#'
#'
#' @slot ranges a list of GenomicRanges GTF objects
#' @slot genome genome sequence
#' @slot domains a dataFrame of predicted protein domains
#' @slot nmd a list of NMD-related objects
#' @slot misc  a list of miscellaneous information
#' @slot version version of factR this object was built under
#'
#' @return A factRObject
#' @name factRObject-class
#' @rdname factRObject-class
#' @exportClass factR
#'
#' @importFrom dplyr %>%
#'
setClass("factR",
         slots = c(
             custom = "list",
             reference = "list",
             ASplicings = "list",
             domains = "data.frame",
             nmd = "data.frame",
             misc = "list",
             version = "character"
         )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Preview custom transcriptome from factR object
#'
#' @description Test test
#' @param object factRObject
#' @param ... Optional: a list of features to view. Input
#' can be a mixture of names from gene_id, gene_name or transcript_id metadata.
#' If missing, function will return all features
#' @param type data to show. Can be "transcripts" or "proteins"
#' @param in_console whether to print custom transcriptome in console
#'
#' @export
#' @rdname preview-methods
#' @author Fursham Hamid
setGeneric("View", function(object, ...,
                            type = "transcripts",
                            in_console = FALSE) standardGeneric("View"))

#' @param object factRObject
#' @param n an integer of the length to display
#'
#' @return GenomicRanges object
#' @export
#' @rdname preview-methods
setGeneric("head", function(object, n = 6L) standardGeneric("head"))

#' @param object factRObject
#' @param n an integer of the length to display
#'
#' @export
#' @rdname preview-methods
setGeneric("tail", function(object,  n = 6L) standardGeneric("tail"))





#' Visualize RNA transcripts or translated proteins
#'
#' @description
#' Plots out transcript architecture or translated protein domains from
#' selected genes
#'
#' @param object  factRObject
#' @param ... Optional: a list of features to plot. Input
#' can be a mixture of names from gene_id, gene_name or transcript_id metadata.
#' If missing, function will plot transcripts/proteins from the first 9 genes
#' of the custom transcriptome
#' @param type data to show. Can be "transcripts" or "proteins"
#' @param rescale_introns when plotting transcripts, whether to rescale introns
#' @param ncol number of columns to combine multiple feature plots to
#'
#' @return ggplot2 object. If multiple genes are detected, plots will be
#' combined using patchwork
#' @export
#' @rdname Plot
#'
#' @author Fursham Hamid
#'
#' @examples
#' data(factRsample)
#' View(factRsample, "Selenop")
setGeneric("Plot", function(object, ..., type = "transcripts",
                            rescale_introns = FALSE, ncol = 1) standardGeneric("Plot"))



#' Test
#'
#' @param object
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
setGeneric("RunfactR", function(object, verbose = FALSE) standardGeneric("RunfactR"))























