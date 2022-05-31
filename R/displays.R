#' @include generics.R

setMethod("show", "factR", function(object){
    print(object@custom)
})

setMethod("summary",
          "factR",
          function(object) {
              ngenes <- length(unique(object@txdata$gene_id))
              ntxs <- length(unique(object@txdata$transcript_id))
              nnovel <- sum(object@txdata$novel == "yes")
              ncds <- sum(object@txdata$cds == "yes")
              cat(sprintf("factR object [version %s]\n", object@version))
              cat(sprintf("## Total number of genes: %s\n", ngenes))
              cat(sprintf("## Total number of transcripts: %s [%s novel]\n", ntxs, nnovel))
              cat(sprintf("## Total number of CDSs: %s \n", ncds))
              cat(sprintf("## Reference species: %s\n", object@reference$species))
              cat(sprintf("## Reference database: %s\n", object@reference$db))
          }
)



# setMethod("head", "factR", function(object, n = 6L) {
#     return(object@custom[1:n])
# })
#
# setMethod("tail", "factR", function(object, n = 6L) {
#     return(utils::tail(object@custom, n))
# })

setMethod("view", "factR", function(object, ...) {

    # check features
    genetxs <- methods::slot(object, name = "txdata")
    txs <- .getTxs(object, ...)


    # display data
    gtf <- methods::slot(object, "custom")
    data.to.view <- gtf[gtf$transcript_id %in% txs]
    data.to.view <- as.data.frame(data.to.view)
    if(nrow(data.to.view) > 1000){
        rlang::warn("Viewing first 1000 entries")
        tibble::view(data.to.view[1:1000], title = "factRObject")
    } else {
        tibble::view(data.to.view, title = "factRObject")
    }

})
setMethod("txRanges", "factR", function(object, ...) {
    gtf <- methods::slot(object, "custom")
    txs <- tryCatch(.getTxs(object, ...),
                error = function(e) rlang::abort("Feature not found"))
    gtf[gtf$transcript_id %in% txs]
})
setMethod("txData", "factR", function(object) {
    methods::slot(object, name = "txdata")
})




setMethod("[", signature("factR"), function (x, i, j){
      # check features
      gtf <- methods::slot(x, "custom")
      genetxs <- methods::slot(x, name = "txdata")
      if(missing(i)){
          gtf <- gtf
      } else if(typeof(i) %in% c("integer", "double")){
          gtf <- gtf[i]
      } else if(typeof(i) %in% c("logical")){
          gtf <- gtf[i]
      } else {
          txs <- tryCatch(.getTxs(x, i),
                          error = function(e) rlang::abort("Feature not found"))
          gtf <- gtf[gtf$transcript_id %in% txs]
      }
      gtf[,j]
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



setMethod("head", "factR", function(x, n = 6L){
    utils::head(x@custom, n)
})
setMethod("tail", "factR", function(x, n = 6L){
    utils::tail(x@custom, n)
})




















