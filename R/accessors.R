#' @include generics.R

setMethod("show", "factR", function(object){
    cat(sprintf("class: factRObject [version %s]\n", object@version))
    cat(sprintf("# transcriptome:\n"))
    ngenes <- length(unique(object@transcriptome$gene_id))
    #ntxs <- length(unique(object@transcriptome$transcript_id))
    #nnovel <- sum(object@txData$novel == "yes")
    #ncds <- sum(object@txData$cds == "yes")
    cat(sprintf("## %s genes\n", ngenes))
    #cat(sprintf("## %s transcripts [%s novel]\n", ntxs, nnovel))
    #cat(sprintf("## %s coding transcripts \n", ncds))
    cat(sprintf("# active set: %s\n", object@active.set))
    cat(sprintf("# samples:\n"))
})

setMethod("summary", "factR", function(object) object )



setMethod("head", "factR", function(x, n = 6L){
    utils::head(x@custom, n)
})
setMethod("tail", "factR", function(x, n = 6L){
    utils::tail(x@custom, n)
})


# setMethod("view", "factR", function(object, ...) {
# 
#     # check features
#     genetxs <- methods::slot(object, name = "txData")
#     txs <- .getTxs(object, ...)
# 
# 
#     # display data
#     gtf <- methods::slot(object, "custom")
#     data.to.view <- gtf[gtf$transcript_id %in% txs]
#     data.to.view <- as.data.frame(data.to.view)
#     if(nrow(data.to.view) > 1000){
#         rlang::warn("Viewing first 1000 entries")
#         tibble::view(data.to.view[1:1000], title = "factRObject")
#     } else {
#         tibble::view(data.to.view, title = "factRObject")
#     }
# 
# })
setMethod("rangesData", "factR", function(object, ..., set = NULL) {
    gtf <- methods::slot(object, "transcriptome")
    txs <- tryCatch(.getTxs(object, ...),
                error = function(e) rlang::abort("Feature not found"))
    gtf[gtf$transcript_id %in% txs]
})




## switch sets on the fly
setMethod("[[", "factR", function(x, i, j){
    if(i %in% names(obj@sets)){
        x@active.set <- i
    } else if(i %in% 1:3){
        x@active.set <- names(obj@sets)[i]
    } else {
        rlang::abort("Incorrect set name or index")
    }
    x
})



setMethod("txData", "factR", function(object) {
    methods::slot(object, name = "txData")
})
setMethod("txData$", "factR", function(object,name) {
    slot(object, "txData")[[name]]
})
setMethod("txData<-", "factR", function(object,value) {
    slot(object, "txData")<- value
    object
})
setMethod("txData$<-", "factR", function(object,...,value) {
    slot(object, "txData")[[...]] <- value
    object
})


setMethod("sampleData", "factR", function(object) {
    methods::slot(object, name = "colData")
})



.DollarNames.factR <- function(x, pattern = "")
    grep(pattern, colnames(as.data.frame(x@custom)), value=TRUE)

setMethod("$", signature("factR"), function (x, name) {
    # check features
    x <- as.data.frame(x@custom)
    x[[name]]
    
})

setMethod("$<-", "factR", function(x, name, value){
    mcols(x@custom)[[name]] <- value
    x
})























