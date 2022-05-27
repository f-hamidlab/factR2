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
             custom = "GenomicRanges",
             txdata = "data.frame",
             ASplicings = "GenomicRanges",
             domains = "data.frame",
             nmd = "data.frame",
             misc = "list",
             reference = "list",
             version = "character"
         )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Preview factR object
#'
#' @description A set of functions to view the contents of a factR object.
#' \itemize{
#'  \item{view(): }{view contents of custom transcriptome in spreadsheet-style format}
#'  \item{head(): }{prints out the first n lines of custom transcriptome}
#'  \item{tail(): }{prints out the last n lines of custom transcriptome}
#'  \item{txData(): }{prints out transcript metadata}
#' }
#'
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
setGeneric("view", function(object, ...,
                            in_console = FALSE) standardGeneric("view"))

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

#' @param object factRObject
#'
#' @export
#' @rdname preview-methods
setGeneric("txData", function(object) standardGeneric("txData"))





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
#' Plot(factRsample, "Selenop")
setGeneric("Plot", function(object, ..., type = "transcripts",
                            rescale_introns = FALSE, ncol = 1) standardGeneric("Plot"))



#' Run factR workflow
#'
#' @description This wrapper function will perform the core factR workflow which
#' include (1) building of CDSs, (2) predicting of NMD-sensitive transcripts
#' and (3) labelling of alternative splicing segments.
#'
#' @param object factRObject
#' @param verbose whether to display workflow progress
#'
#' @return Updated factRObject
#' @export
#' @seealso \code{\link{BuildCDS}}, \code{\link{PredictNMD}}, \code{\link{FindAltSplicing}}
#'
#' @rdname RunfactR
#' @examples
#' data(factRsample)
#' factRsample <- RunfactR(factRsample)
setGeneric("RunfactR", function(object, verbose = FALSE) standardGeneric("RunfactR"))


#' Build coding sequence on custom transcriptome
#'
#' @description This function constructs coding sequences on custom transcriptomes.
#' This is done based on the core functionality of factR's buildCDS function.
#'
#' @param object factRObject
#' @param verbose whether to display workflow progress
#'
#' @return Updated factRObject
#' @export
#' @seealso \code{\link[factR]{buildCDS}}, \code{\link{RunfactR}}
#'
#' @rdname BuildCDS
#' @examples
#' data(factRsample)
#' factRsample <- BuildCDS(factRsample)
setGeneric("BuildCDS", function(object, verbose = FALSE) standardGeneric("BuildCDS"))


#' Predict NMD sensitivity on coding transcripts
#'
#' @description This function predicts the sensitivity of protein-coding
#' transcripts to nonsense-mediated decay.
#'
#' @param object factRObject
#' @param NMD_threshold Minimum distance of stop_codon to last exon junction (EJ)
#' which triggers NMD. Default = 50bp
#' @param verbose whether to display workflow progress
#'
#' @return Updated factRObject
#' @export
#' @seealso \code{\link[factR]{predictNMD}}, \code{\link{RunfactR}}
#'
#' @rdname PredictNMD
#' @examples
#' data(factRsample)
#' factRsample <- BuildCDS(factRsample)
#' factRsample <- PredictNMD(factRsample)
setGeneric("PredictNMD", function(object, NMD_threshold = 50, verbose = FALSE) standardGeneric("PredictNMD"))

#' Find alternative splicing events
#'
#' @description This function identifies all alternative splicing events that
#' occur in the custom transcriptome.
#'
#' @param object factRObject
#'
#' @return Updated factRObject
#' @export
#' @seealso \code{\link{RunfactR}}
#'
#' @rdname FindAltSplicing
#' @examples
#' data(factRsample)
#' factRsample <- FindAltSplicing(factRsample)
setGeneric("FindAltSplicing", function(object) standardGeneric("FindAltSplicing"))























