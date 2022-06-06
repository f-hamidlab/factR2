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
#' @importFrom methods slot
#'
setClass("factR",
         slots = c(
             sets = "list",
             transcriptome = "GenomicRanges",
             colData = "data.frame",
             domains = "list",
             design = "formula",
             active.set = "character",
             reference = "list",
             misc = "list",
             version = "character"
         )
)

setClass("factRset",
         slots = c(
             counts = "matrix",
             data = "matrix",
             rowData = "data.frame",
             comparisons = "data.frame"
         )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

################################################################################
## General object preview
################################################################################




################################################################################
## Sets-related functions
################################################################################

setGeneric("activeSet", function(object) standardGeneric("activeSet"))
setGeneric("activeSet<-", function(object, value) standardGeneric("activeSet<-"))
setGeneric("listSets", function(object) standardGeneric("listSets"))

################################################################################
## Previewing rangesData
################################################################################

#' @export
#' @rdname preview-methods
setGeneric("rangesData", function(object, ..., set = NULL)
    standardGeneric("rangesData"))

################################################################################
## Previewing and modifying featureData
################################################################################

#' Preview factR object
#'
#' @description A set of functions to view the contents of a factR object.
#' \itemize{
#'  \item{view(): }{view contents of custom transcriptome in spreadsheet-style format}
#'  \item{txRanges(): }{prints out transcript ranges}
#'  \item{featureData(): }{prints out transcript metadata}
#' }
#'
#' @param object factRObject
#'
#' @export
#' @rdname preview-methods
setGeneric("featureData", function(object, ..., set = NULL) standardGeneric("featureData"))

setGeneric("featureData$", function(object, ...) standardGeneric("featureData$"))
setGeneric("featureData<-", function(object, value) standardGeneric("featureData<-"))
setGeneric("featureData$<-", function(object, ..., value) standardGeneric("featureData$<-"))
setGeneric("addFeatureData", function(object, ..., set = NULL) standardGeneric("addFeatureData"))

################################################################################
## Previewing and modifying sampleData
################################################################################

setGeneric("sampleData", function(object) standardGeneric("sampleData"))
setGeneric("sampleData$", function(object, ...) standardGeneric("sampleData$"))
setGeneric("sampleData<-", function(object, value) standardGeneric("sampleData<-"))
setGeneric("sampleData$<-", function(object, ..., value) standardGeneric("sampleData$<-"))
setGeneric("addSampleData", function(object, ..., set = NULL) standardGeneric("addSampleData"))




#' Visualize custom transcriptome and proteins
#'
#' @description
#' Plots out transcripts or protein domains from custom transcriptome
#'
#' @param object  factRObject
#' @param ... Optional: a list of features to plot. Input
#' can be a mixture of names from gene_id, gene_name or transcript_id metadata.
#' If missing, function will plot transcripts from the first 9 genes
#' of the custom transcriptome
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
#' # plot transcripts
#' data(factRsample)
#' plotTranscripts(factRsample, "Selenop")
setGeneric("plotTranscripts", function(object, ...,
                            rescale_introns = FALSE,
                            ncol = 1) standardGeneric("plotTranscripts"))



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
#' @seealso \code{\link{buildCDS}}, \code{\link{predictNMD}}, \code{\link{findAltSplicing}}
#'
#' @rdname RunfactR
#' @examples
#' data(factRsample)
#' factRsample <- runfactR(factRsample)
setGeneric("runfactR", function(object, verbose = FALSE) standardGeneric("runfactR"))
setMethod("runfactR", "factR", function(object, verbose = FALSE) {
    object <- buildCDS(object, verbose)
    object <- predictNMD(object, verbose)
    object <- getAAsequence(object, verbose)
    object
})


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
#' @seealso \code{\link[factR]{buildCDS}}, \code{\link{runfactR}}
#'
#' @rdname buildCDS
#' @examples
#' data(factRsample)
#' factRsample <- buildCDS(factRsample)
setGeneric("buildCDS", function(object, verbose = FALSE) standardGeneric("buildCDS"))


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
#' @seealso \code{\link[factR]{predictNMD}}, \code{\link{runfactR}}
#'
#' @rdname predictNMD
#' @examples
#' data(factRsample)
#' factRsample <- buildCDS(factRsample)
#' factRsample <- predictNMD(factRsample)
setGeneric("predictNMD", function(object, NMD_threshold = 50, verbose = FALSE) standardGeneric("predictNMD"))

#' Find alternative splicing events
#'
#' @description This function identifies all alternative splicing events that
#' occur in the custom transcriptome.
#'
#' @param object factRObject
#'
#' @return Updated factRObject
#' @export
#' @seealso \code{\link{runfactR}}
#'
#' @rdname findAltSplicing
#' @examples
#' data(factRsample)
#' factRsample <- findAltSplicing(factRsample)
setGeneric("findAltSplicing", function(object) standardGeneric("findAltSplicing"))

#' Predict functional domains
#'
#' @description This function predicts the functional domains found on
#' proteins
#'
#' @param object factRObject
#' @param ... a list of features to perform prediction on. Input
#' can be a mixture of names from gene_id, gene_name or transcript_id metadata.
#' If missing, function will predict domains on all proteins in the GTF
#' @param database character string specifying the database to search. Currently
#' support "superfamily" (Default) and "pfam"
#' @param ncores number of cores to run function (Default: 4)
#'
#' @return Updated factRObject
#' @export
#'
#' @rdname predictDomain
#' @examples
#' data(factRsample)
#' factRsample <- buildCDS(factRsample)
#' factRsample <- predictNMD(factRsample)
setGeneric("predictDomain", function(object, ...,
                                     database = "superfamily",
                                     ncores = 4) standardGeneric("predictDomain"))



#' @param object  factRObject
#' @param ... Optional: a list of features to plot. Input
#' can be a mixture of names from gene_id, gene_name or transcript_id metadata.
#' If missing, function will plot transcripts from the first 9 genes
#' of the custom transcriptome
#' @param ncol number of columns to combine multiple feature plots to
#'
#' @export
#' @rdname Plot
#'
#' @author Fursham Hamid
#'
#' @examples
#'
#' # plot protein domains
#' factRsample <- runfactR(factRsample)
#' plotDomains(factRsample, "Prkaa1")
setGeneric("plotDomains", function(object, ..., ncol = 1) standardGeneric("plotDomains"))


#' Get protein coding sequences
#'
#' @description This function will translate coding transcript sequences into
#' amino acid sequences and stores them as a dataframe in the object's
#' domain slot
#'
#' @param object factRObject
#' @param verbose whether to display workflow progress
#'
#' @return Updated factRObject
#' @export
#' @seealso \code{\link{runfactR}}
#'
#' @rdname getAAsequence
#' @examples
#' data(factRsample)
#' factRsample <- buildCDS(factRsample)
#' factRsample <- getAAsequence(factRsample)
setGeneric("getAAsequence", function(object, verbose = FALSE) standardGeneric("getAAsequence"))




#' Identify AS-NMD events
#'
#' @description This function will xxx
#'
#' @param object factRObject
#'
#' @return Updated factRObject
#' @export
#' @seealso \code{\link{runfactR}}
#'
#' @rdname testASNMDevents
#' @examples
#' data(factRsample)
#' factRsample <- buildCDS(factRsample)
#' factRsample <- testASNMDevents(factRsample)
setGeneric("testASNMDevents", function(object) standardGeneric("testASNMDevents"))







