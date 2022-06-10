#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' factRObject
#'
#' @description
#' A container for storing and processing custom transcriptome.
#'
#' @details
#' The factRObject is a representation of custom-assembled transcriptome data.
#' Coordinate and metadata information from an input transcriptome (GTF) is
#' parsed and grouped into 3 data sets which corresponds to the levels of gene,
#' transcript and alternative splicing (AS). With the aid of a compatible
#' reference genome and expression data, functional and
#' comparative analyses can be performed to better understand the complexity
#' of a transcriptome.
#'
#' @section Constructor:
#' \describe{
#'    \item{factRObject can be easily constructed as such:}{
#'       \code{createfactRObject(gtf, reference)}
#'    }
#' }
#'
#'
#'
#' @section Accessors:
#' Test
#'
#'
#
#'
#' @docType class
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







#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
## Display ====
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

### General ####
setMethod("show", "factR", function(object) show.factR(object))
setMethod("summary", "factR", function(object) object )

# head and tail previews featureData of current set
setMethod("head", "factR", function(x, n = 6L){
    utils::head(x[[]], n = n)
})
setMethod("tail", "factR", function(x, n = 6L){ utils::tail(x[[]], n = n) })

setMethod("dim", "factR", function(x){ dim(x@sets[[x@active.set]]@data) })
setMethod("nrow", "factR", function(x){ base::nrow(x@sets[[x@active.set]]@data) })
setMethod("ncol", "factR", function(x){ base::ncol(x@sets[[x@active.set]]@data) })





### Sets ####
setGeneric("activeSet", function(object) standardGeneric("activeSet"))
setMethod("activeSet", "factR", function(object){
    methods::slot(object, "active.set")
})
setGeneric("listSets", function(object) standardGeneric("listSets"))
setMethod("listSets", "factR", function(object){ names(object@sets) })





### GRanges ####
setGeneric("granges", function(object, ..., set = NULL) standardGeneric("granges"))
setMethod("granges", "factR", function(object, ..., set = NULL) {
    granges.factR(object, ..., set = set)
})


### Features ####
setMethod("[[", "factR", function(x, i){
    if(missing(i)){
        x@sets[[x@active.set]]@rowData
    } else if(i %in% names(x@sets)){
        x@sets[[i]]@rowData
    } else {
        rlang::abort("Incorrect set name or index")
    }
})

# feature preview with option to subset data
setGeneric("features", function(object, ..., set = NULL) standardGeneric("features"))
setMethod("features", "factR", function(object, ..., set = NULL) {
    if(is.null(set)){
        set <- slot(object, "active.set")
    } else if(!set %in% listSets(object)){
        rlang::warn("Incorrect set name or index, using active set")
        set <- slot(object, "active.set")
    }

    dat <- slot(object@sets[[set]], "rowData")
    out.type <- ifelse(set == "transcript", "transcript_id", "gene_id")
    feat <- .getFeat(object, ..., out = out.type)
    return(dat[dat[[out.type]] %in% feat,])
})


### Samples ####
setGeneric("samples", function(object) standardGeneric("samples"))
setMethod("samples", "factR", function(object) { object@colData })

### Counts ####

### Design ####
setGeneric("design", function(object) standardGeneric("design"))
setMethod("design", "factR", function(object) object@design)


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
## Setters ====
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

### Set active set ####
setGeneric("activeSet<-", function(object, value) standardGeneric("activeSet<-"))
setMethod("activeSet<-", "factR", function(object, value){
    methods::slot(object, "active.set") <- value
    object
})

### Set design ####
setGeneric("design<-", function(object, value) standardGeneric("design<-"))


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
## Subsetters ====
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
## Plotters ====
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

### Transcripts ####
setGeneric("plotTranscripts", function(object, ...,
                                       rescale_introns = FALSE,
                                       ncol = 1) standardGeneric("plotTranscripts"))

### Domains ####
setGeneric("plotDomains", function(object, ..., ncol = 1) standardGeneric("plotDomains"))


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
## Validty ====
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX












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


setGeneric("addTxCounts", function(object, countData, sampleData = NULL, design = NULL, verbose = FALSE) standardGeneric("addTxCounts"))


































#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# END
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
