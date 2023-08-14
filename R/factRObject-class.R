# TODO: add accessor functions for expression data
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
#' parsed and grouped into 3 data Sets which corresponds to the levels of gene (genes),
#' transcript (txs) and alternative splicing (ase). With the aid of a compatible
#' reference genome and expression data, functional and
#' comparative analyses can be performed to better understand the complexity
#' of a transcriptome.
#'
#' @section Constructor:
#' \describe{
#'      \item{factRObject can be easily constructed as such:}{
#'          \code{\code{\link{createfactRObject}}(gtf, reference)}
#'    }
#' }
#'
#'
#'
#' @section Accessors:
#' Interact with a factRObject (x) the following ways:
#' \describe{
#'    \item{}{
#'       \code{summary(x): Prints preview of factRObject}
#'    }
#'    \item{}{
#'       \code{head(x): Prints the first 6 metadata entries of the active Set}
#'    }
#'    \item{}{
#'       \code{tail(x): Prints the last 6 metadata entries of the active Set}
#'    }
#'    \item{}{
#'       \code{dim(x): Prints the dimension of the active set (features x samples)}
#'    }
#'    \item{}{
#'       \code{nrow(x): Returns the number of features of the active set}
#'    }
#'    \item{}{
#'       \code{ncol(x): Returns the number of samples}
#'    }
#' }
#'
#' A factRObject consists of 3 `Sets` of data: gene, transcript and AS. At any
#' one time, one of these `Sets` will be made the active Set and some factR2
#' functions will be called only on this active Set. Below are some ways to
#' check and change the active Set:
#' \describe{
#'    \item{}{
#'       \code{\code{\link{listSets}}(x): Lists Sets in object}
#'    }
#'    \item{}{
#'       \code{\code{\link{activeSet}}(x): Returns the active Set of object}
#'    }
#'    \item{}{
#'       \code{\code{\link{activeSet}}(x)<-: Change the active Set}
#'    }
#' }
#'
#' Each `Set` of an factRObject hold rich information of the features that it
#' contains. Below are some convenient wrappers to display the metadata in
#' each `Set`.
#' \describe{
#'    \item{}{
#'       \code{\code{\link{features}}(x): Displays metadata of active Set}
#'    }
#'    \item{}{
#'       \code{\code{\link{genes}}(x) or \code{\link{gns}}(x): Displays gene metadata}
#'    }
#'    \item{}{
#'       \code{\code{\link{transcripts}}(x) or \code{\link{txs}}(x): Displays transcript metadata}
#'    }
#'    \item{}{
#'       \code{\code{\link{ase}}(x): Displays alternative splicing events metadata}
#'    }
#' }
#'
#' factR2 also contain functions to plot exon-level and domain-level
#' structures as an interactive plot:
#' \describe{
#'    \item{}{
#'       \code{\code{\link{plotTranscripts}}(x, "gene of interest"): Plots transcripts}
#'    }
#'    \item{}{
#'       \code{\code{\link{plotDomains}}(x, "gene of interest")): Plots domains}
#'    }
#' }
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
#' @export
setMethod("show", "factR", function(object) show.factR(object))

#' @export
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




### GRanges ####

setGeneric("granges", function(object, ..., set = NULL) standardGeneric("granges"))
#' Display factR2 object data
#'
#' @description
#' A factRObject-class contains different types of data at the gene, transcript,
#'  alternative splicing (AS) and protein domain levels. The functions below are
#'  designed to display specific contents of a factRObject.
#'
#'
#' @param object factRObject
#' @param ... One or more features to display. Can be the following:
#' \itemize{
#'  \item{gene_id: }{ID of gene to plot}
#'  \item{gene_name: }{Name of gene to plot}
#'  \item{transcript_id: }{ID of transcript to plot}
#' }
#' @param set Set metadata to display. Can be "gene", "transcript" or "AS".
#'
#' @export
#' @return
#' \itemize{
#'  \item{`granges`: }{GenomicRanges object of selected features}
#'  \item{`activeSet` and `listSets`: }{Character value/vector}
#'  \item{All other functions: }{Tibble dataframe containing metadata of selected features}
#' }
#'
#'
#' @name factR-meta
#' @rdname factR-meta
#' @examples
#' ### Load sample factRObject
#' data("factRsample")
#'
#' ## Prints out activeSet
#' activeSet(factRsample)
#'
#' ## Change activeSet
#' activeSet(factRsample) <- "transcript"
#'
#' ## Returns coordinates and metadata of features as a GenomicRanges object
#' granges(factRsample)   # from activeSet
#' granges(factRsample, "Dab2")   # specific features
#' granges(factRsample, "Dab2", set = "gene")   # specific features from different Set
#'
#' ## Returns metadata of features
#' features(factRsample)   # from activeSet
#' features(factRsample, "Dab2")   # specific features
#' features(factRsample, "Dab2", set = "gene")   # specific features from different Set
#'
#' ### This is the same as:
#' genes(factRsample, "Dab2")
#'
#'
#' ## To return protein-coding domains, the protein-coding domains need to be predicted first:
#' factRsample <- buildCDS(factRsample)
#' factRsample <- getAAsequence(factRsample)
#' factRsample <- predictDomains(factRsample, "Dab2")
#'
#' ## Then, the domains of the selected gene can be printed as such:
#' domains(factRsample, "Dab2")
#'
#' ## All outputs can be assigned to a variable and manipulated further using other functions:
#' ase(factRsample) %>% dplyr::filter(AStype == "CE")
#'
#'

setMethod("granges", "factR", function(object, ..., set = NULL) {
    granges.factR(object, ..., set = set)
})



### Sets ####

setGeneric("activeSet", function(object) standardGeneric("activeSet"))
#' @param object factRObject
#' @export
#' @rdname factR-meta
setMethod("activeSet", "factR", function(object){
    methods::slot(object, "active.set")
})

setGeneric("activeSet<-", function(object, value) standardGeneric("activeSet<-"))
#' @param object factRObject
#' @param value Character value of one of the following: "gene", "transcript" or "AS"
#' @export
#' @rdname factR-meta
setMethod("activeSet<-", "factR", function(object, value){
    object@active.set <- value
})

#' @export
setGeneric("listSets", function(object) standardGeneric("listSets"))
#' @param object factRObject
#' @export
#' @rdname factR-meta
setMethod("listSets", "factR", function(object){ names(object@sets) })








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
setGeneric("features", function(object, ..., set = NULL,
                                show_more = FALSE) standardGeneric("features"))
#' @param object factRObject
#' @param ... One or more features to display. Can be the following:
#' \itemize{
#'  \item{gene_id: }{ID of gene to plot}
#'  \item{gene_name: }{Name of gene to plot}
#'  \item{transcript_id: }{ID of transcript to plot}
#' }
#' @param set Set metadata to display. Can be "gene", "transcript" or "AS".
#' @export
#' @rdname factR-meta
setMethod("features", "factR", function(object, ..., set = NULL,
                                        show_more = FALSE) {
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
    if(!show_more){
        .msginfo("Set `show_more to TRUE to show more info`")
        return(tibble::as_tibble(dat))
    } else {
        return(dat)
    }


})

# wrappers to quickly get genes, transcripts and AS
setGeneric("genes", function(object, ...,
                             show_more = FALSE) standardGeneric("genes"))
#' @param object factRObject
#' @param ... One or more features to display. Can be the following:
#' \itemize{
#'  \item{gene_id: }{ID of gene to plot}
#'  \item{gene_name: }{Name of gene to plot}
#'  \item{transcript_id: }{ID of transcript to plot}
#' }
#' @export
#' @rdname factR-meta
setMethod("genes", "factR", function(object, ...,
                                     show_more = FALSE) {

    return((features(object,..., set="gene")))
})

setGeneric("gns", function(object, ...,
                           show_more = FALSE) standardGeneric("gns"))
#' @param object factRObject
#' @param ... One or more features to display. Can be the following:
#' \itemize{
#'  \item{gene_id: }{ID of gene to plot}
#'  \item{gene_name: }{Name of gene to plot}
#'  \item{transcript_id: }{ID of transcript to plot}
#' }
#' @export
#' @rdname factR-meta
setMethod("gns", "factR", function(object, ...,
                                   show_more = FALSE) {

    return((features(object,..., set="gene",
                     show_more = show_more)))
})


setGeneric("transcripts", function(object, ...,
                                   show_more = FALSE) standardGeneric("transcripts"))
#' @param object factRObject
#' @param ... One or more features to display. Can be the following:
#' \itemize{
#'  \item{gene_id: }{ID of gene to plot}
#'  \item{gene_name: }{Name of gene to plot}
#'  \item{transcript_id: }{ID of transcript to plot}
#' }
#' @export
#' @rdname factR-meta
setMethod("transcripts", "factR", function(object, ...,
                                           show_more = FALSE) {

    return((features(object,..., set="transcript",
                     show_more = show_more)))
})

setGeneric("txs", function(object, ...,
                           show_more = FALSE) standardGeneric("txs"))
#' @param object factRObject
#' @param ... One or more features to display. Can be the following:
#' \itemize{
#'  \item{gene_id: }{ID of gene to plot}
#'  \item{gene_name: }{Name of gene to plot}
#'  \item{transcript_id: }{ID of transcript to plot}
#' }
#' @export
#' @rdname factR-meta
setMethod("txs", "factR", function(object, ...,
                                   show_more = FALSE) {

    return((features(object,..., set="transcript",
                     show_more = show_more)))
})


setGeneric("ase", function(object, ...,
                           show_more = FALSE) standardGeneric("ase"))
#' @param object factRObject
#' @param ... One or more features to display. Can be the following:
#' \itemize{
#'  \item{gene_id: }{ID of gene to plot}
#'  \item{gene_name: }{Name of gene to plot}
#'  \item{transcript_id: }{ID of transcript to plot}
#' }
#' @export
#' @rdname factR-meta
setMethod("ase", "factR", function(object, ...,
                                   show_more = FALSE) {

    return((features(object,..., set="AS",
                     show_more = show_more)))
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
setGeneric("domains", function(object, ...) standardGeneric("domains"))
#' @param object factRObject
#' @param ... One or more features to display. Can be the following:
#' \itemize{
#'  \item{gene_id: }{ID of gene to plot}
#'  \item{gene_name: }{Name of gene to plot}
#'  \item{transcript_id: }{ID of transcript to plot}
#' }
#' @export
#' @rdname factR-meta
setMethod("domains", "factR", function(object, ...){
    feat <- .getFeat(object, ...)
    domains <- object@domains$data
    domains <- domains[domains$transcript_id %in% feat & domains$type == "DOMAIN",]
    tibble::as_tibble(domains)
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
setGeneric("plotTranscripts", function(object, ...,
                                       rescale_introns = FALSE,
                                       ncol = 1) standardGeneric("plotTranscripts"))

### Domains ####
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

setGeneric("exportGTF", function(object, out=getwd()) standardGeneric("exportGTF"))

#' Exporting data from factRObject
#' @param object factRObject
#' @param out Path to output directory or path to output file. If the former,
#' the output file will be named "factR.gtf"
#'
#' @export
#' @rdname factR-export
setMethod("exportGTF", "factR", function(object, out=getwd()){
    export.factR(object, out, "gtf")
})

setGeneric("exportTable", function(object, out=getwd(), data = "AS") standardGeneric("exportTable"))
#' @param object factRObject
#' @param out Path to output directory or path to output file. If the former,
#' the output file will be named "factR.gtf"
#' @param data Set metadata to export Can be "gene", "transcript" or "AS".
#'
#' @export
#' @rdname factR-export
setMethod("exportTable", "factR", function(object, out=getwd(), data = activeSet(object)){
    export.factR(object, out, data)
})

setGeneric("exportAll", function(object, path=getwd()) standardGeneric("exportAll"))
#' @param object factRObject
#' @param path Path to output directory or path to output file.
#'
#' @export
#' @rdname factR-export
setMethod("exportAll", "factR", function(object, path=getwd()){
    export.factR(object, path, "gtf")
    export.factR(object, path, "gene")
    export.factR(object, path, "transcript")
    export.factR(object, path, "AS")
})




















































#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# END
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
