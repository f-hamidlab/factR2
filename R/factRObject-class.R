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
#' @import ggplot2
#' @importFrom crayon green italic blue bold cyan white red
#'
setClass("factR",
         slots = c(
             sets = "list",
             transcriptome = "GenomicRanges",
             colData = "data.frame",
             domains = "list",
             project = "character",
             active.set = "character",
             active.ident = "character",
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
             misc = "list"
         )
)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
## Display ====
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

### General ####
#' Preview factR object
#'
#' @param object factRObject
#' @export
#'
#'
setMethod("show", "factR", function(object) show.factR(object))

#' Preview factR object
#' #' @param object factRObject
#' @export
#'
#'
setMethod("summary", "factR", function(object) object )

# head and tail previews featureData of current set

#' @export
setMethod("head", "factR", function(x, n = 6L){
    utils::head(x[[]], n = n)
})
#' @export
setMethod("tail", "factR", function(x, n = 6L){ utils::tail(x[[]], n = n) })
#' @export
setMethod("dim", "factR", function(x){ dim(x@sets[[x@active.set]]@data) })

#' @export
setMethod("nrow", "factR", function(x){ base::nrow(x@sets[[x@active.set]]@data) })

#' @export
setMethod("ncol", "factR", function(x){ base::ncol(x@sets[[x@active.set]]@data) })





### Sets ####
#' @export
setGeneric("activeSet", function(object) standardGeneric("activeSet"))
setMethod("activeSet", "factR", function(object){
    methods::slot(object, "active.set")
})

#' @export
setGeneric("activeSet<-", function(object, value) standardGeneric("activeSet<-"))
setMethod("activeSet<-", "factR", function(object, value){
    object@active.set <- value
})

#' @export
setGeneric("listSets", function(object) standardGeneric("listSets"))
setMethod("listSets", "factR", function(object){ names(object@sets) })





### GRanges ####
#' @export
setGeneric("granges", function(object, ..., set = NULL) standardGeneric("granges"))
setMethod("granges", "factR", function(object, ..., set = NULL) {
    granges.factR(object, ..., set = set)
})


### Features ####
#' @export
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
#' @export
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
    dat <- dat[dat[[out.type]] %in% feat,]
    #rownames(dat) <- NULL
    return(dat)
})

# wrappers to quickly get genes, transcripts and AS
#' @export
setGeneric("genes", function(object, ...) standardGeneric("genes"))
setMethod("genes", "factR", function(object, ...) {

    return(tibble::as_tibble(features(object,..., set="gene")))
})

#TODO: set shortname for genes to gns


#' @export
setGeneric("transcripts", function(object, ...) standardGeneric("transcripts"))
setMethod("transcripts", "factR", function(object, ...) {

    return(tibble::as_tibble(features(object,..., set="transcript")))
})

#' @export
setGeneric("txs", function(object, ...) standardGeneric("txs"))
setMethod("txs", "factR", function(object, ...) {

    return(tibble::as_tibble(features(object,..., set="transcript")))
})

#' @export
setGeneric("ase", function(object, ...) standardGeneric("ase"))
setMethod("ase", "factR", function(object, ...) {

    return(tibble::as_tibble(features(object,..., set="AS")))
})



#' @export
setGeneric("rowNames", function(object) standardGeneric("rowNames"))
setMethod("rowNames", "factR", function(object) {
    rownames(object@sets[[object@active.set]]@rowData) })


### Samples ####
#' @export
setGeneric("samples", function(object) standardGeneric("samples"))
setMethod("samples", "factR", function(object) { object@colData })

#' @export
setGeneric("colNames", function(object) standardGeneric("colNames"))
setMethod("colNames", "factR", function(object) { rownames(object@colData) })

#' @export
setGeneric("colNames<-", function(object, value) standardGeneric("colNames<-"))
setMethod("colNames<-", "factR", function(object, value) {
    object@colData$old.names <- rownames(object@colData)
    rownames(object@colData) <- value
    return(.updatefactR(object))
})


#' @export
setGeneric("ident", function(object) standardGeneric("ident"))
setMethod("ident", "factR", function(object) object@active.ident)

### Counts ####
#' @export
setGeneric("counts", function(object, ..., set = NULL, slot = "data") standardGeneric("counts"))
setMethod("counts", "factR", function(object, ..., set = NULL) {
    if(is.null(set)){
        set <- slot(object, "active.set")
    } else if(!set %in% listSets(object)){
        rlang::warn("Incorrect set name or index, using active set")
        set <- slot(object, "active.set")
    }

    if(set == "AS"){
        dat <- slot(object@sets[[set]], "data")
        feat <- .getFeat(object, ..., out = "AS_id")
        return(dat[feat,])

    } else {
        dat <- slot(object@sets[[set]], slot)
        out.type <- ifelse(set == "transcript", "transcript_id", "gene_id")
        feat <- .getFeat(object, ..., out = out.type)
        return(dat[feat,])
    }
})



### Domains ####
#' @export
setGeneric("domains", function(object, ...) standardGeneric("domains"))
setMethod("domains", "factR", function(object, ...){
    feat <- .getFeat(object, ...)
    domains <- object@domains$data
    domains[domains$transcript_id %in% feat & domains$type == "DOMAIN",]
})


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
## Setters ====
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

### Set active set ####
#' @export
setGeneric("activeSet<-", function(object, value) standardGeneric("activeSet<-"))
setMethod("activeSet<-", "factR", function(object, value){
    methods::slot(object, "active.set") <- value
    object
})

### Set design ####
#' @export
setGeneric("ident<-", function(object, value) standardGeneric("ident<-"))
setMethod("ident<-", "factR", function(object, value){
    methods::slot(object, "active.ident") <- value
    object
})


### Add features/sample data ####
#' @export
setGeneric("addMeta", function(object, meta="samples", data=NULL, ...) standardGeneric("addMeta"))
setMethod("addMeta", "factR", function(object, meta="samples", data=NULL, ...){
    mutate.factR(object, meta, data, ...)
})


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
## Subsetters ====
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#' @export
setGeneric("select", function(object, ...) standardGeneric("select"))
setMethod("select", "factR", function(object, ..., data = "samples"){
    select.factR(object, ...)
})
setMethod("[", "factR", function(x, i, j, set = "sample"){
    if(set=="sample"){
        x@colData[[i]]
    }
    #select.factR(x, j)
})


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
## Plotters ====
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

### Transcripts ####
#' @export
setGeneric("plotTranscripts", function(object, ...,
                                       rescale_introns = FALSE,
                                       ncol = 1) standardGeneric("plotTranscripts"))

### Domains ####
#' @export
setGeneric("plotDomains", function(object, ..., ncol = 1) standardGeneric("plotDomains"))


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
## Validty ====
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#' @export
setGeneric("checkfactR", function(object) standardGeneric("checkfactR"))
setMethod("checkfactR", "factR", function(object){
    checkfactR.factR(object)
})



#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
## Exporters ====
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## TODO: do up export functions (exportGTF, exportTable)
#' @export
setGeneric("exportGTF", function(object, out=getwd()) standardGeneric("exportGTF"))
setMethod("exportGTF", "factR", function(object, out=getwd()){
    export.factR(object, out, "gtf")
})

#' @export
setGeneric("exportTable", function(object, out=getwd(), data = "AS") standardGeneric("exportTable"))
setMethod("exportTable", "factR", function(object, out=getwd(), data = "AS"){
    export.factR(object, out, data)
})

#' @export
setGeneric("exportAll", function(object, path=getwd()) standardGeneric("exportAll"))
setMethod("exportAll", "factR", function(object, path=getwd()){
    export.factR(object, path, "gtf")
    export.factR(object, path, "gene")
    export.factR(object, path, "transcript")
    export.factR(object, path, "AS")
})




















































#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# END
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
