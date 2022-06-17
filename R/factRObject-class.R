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

### Add features/sample data ####
setGeneric("mutate", function(object, ..., data = "samples") standardGeneric("mutate"))
setMethod("mutate", "factR", function(object, ..., data = "samples"){
    mutate.factR(object, ..., data = data)
})


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





























setGeneric("addTxCounts", function(object, countData, sampleData = NULL, design = NULL, verbose = FALSE) standardGeneric("addTxCounts"))


































#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# END
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%