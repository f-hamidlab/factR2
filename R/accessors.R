#' @include generics.R

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Preview #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# show and summary prints a summary of the object
setMethod("show", "factR", function(object){
    cat(sprintf("class: factRObject [version %s]\n", object@version))
    cat(sprintf("# transcriptome: "))
    ngenes <- length(object[["gene"]]$gene_id)
    ntxs <- length(object[["transcript"]]$transcript_id)
    nnovel <- sum(object[["transcript"]]$novel == "yes")
    ncds <- sum(object[["transcript"]]$cds == "yes")
    cat(sprintf("%s genes; ", ngenes))
    cat(sprintf("%s transcripts [%s novel]; ", ntxs, nnovel))
    cat(sprintf("%s coding transcripts \n", ncds))
    cat(sprintf("# active set: %s\n", object@active.set))

    nsamples <- nrow(object@colData)
    samples <- rownames(object@colData)
    if(nsamples > 4){
        samples <- paste(c(samples[1], samples[2], " ... ", samples[nsamples-1],
                          samples[nsamples]), collapse = "  ")
    } else {
        samples <- paste(samples, collapse = ", ")
    }
    cat(sprintf("# samples (%s): %s\n", nsamples, samples))
})
setMethod("summary", "factR", function(object) object )

# head and tail previews featureData of current set
setMethod("head", "factR", function(x, n = 6L){
    utils::head(x[[]], n = n)
})
setMethod("tail", "factR", function(x, n = 6L){
    utils::tail(x[[]], n = n)
})

setMethod("dim", "factR", function(x){
    dim()
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Sets Data #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



## dedicated function to switch sets
setMethod("activeSet", "factR", function(object){
    methods::slot(object, "active.set")
})
setMethod("activeSet<-", "factR", function(object, value){
    methods::slot(object, "active.set") <- value
    object
})
## list sets
setMethod("listSets", "factR", function(object){
    names(object@sets)
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Ranges Data #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMethod("granges", "factR", function(object, ..., set = NULL) {
    if(is.null(set)){
        set <- slot(object, "active.set")
    } else if(!set %in% c(listSets(object), "all")){
        rlang::warn("Incorrect set name or index, using active set")
        set <- slot(object, "active.set")
    }
    gtf <- methods::slot(object, "transcriptome")
    out.type <- ifelse(set %in% c("transcript"), "transcript_id", "gene_id")
    feat <- .getFeat(object, ..., out = out.type)

    if(set == "all"){
        return(gtf[S4Vectors::mcols(gtf)[[out.type]] %in% feat])
    } else {
        return(gtf[S4Vectors::mcols(gtf)[[out.type]] %in% feat & gtf$type %in% set])
    }

})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Feature Data #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## print features quickly
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

# # access and modify columns in featureData
# setMethod("featureData$", "factR", function(object,name) {
#     object@sets[[object@active.set]]@rowData[[name]]
# })
# setMethod("featureData<-", "factR", function(object,value) {
#     object@sets[[object@active.set]]@rowData <- value
#     object
# })
# setMethod("featureData$<-", "factR", function(object, name, value) {
#     object@sets[[object@active.set]]@rowData[[name]] <- value
#     object
# })

# official function to add feature metadata
setMethod("addFeatureData", "factR", function(object, data, colname = NULL, set = NULL) {
    if(is.null(set)){
        set <- slot(object, "active.set")
    } else if(!set %in% listSets(object)){
        rlang::warn("Incorrect set name or index, using active set")
        set <- slot(object, "active.set")
    }

    if(any(class(data) %in% c("data.frame", "tbl", "tbl_df"))){
        # check for dimensions
        if(nrow(data) != nrow(object@sets[[set]]@rowData)){
            rlang::abort(sprintf("Replacement length (%s) not equal to data length (%s)",
                                 nrow(data), nrow(object@sets[[set]]@rowData)))
        }

        # check for rownames
        if(all(rownames(data) %in% rownames(object@sets[[set]]@rowData))){
            featurelvls <- rownames(object@sets[[set]]@rowData)
            object@sets[[set]]@rowData <- cbind(object@sets[[set]]@rowData,
                                                data[featurelvls,])
        } else {
            object@sets[[set]]@rowData <- cbind(object@sets[[set]]@rowData,
                                                data)
        }
    } else {
        if(is.null(colname)){
            colname <- "featureData"
        }
        if(colname %in% colnames(object@sets[[set]]@rowData)){
            colname <- paste0(colname, ".1")
        }
        object@sets[[set]]@rowData[[colname]] <- data
    }
    object
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Sample Data #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# sample preview with option to subset data
setMethod("samples", "factR", function(object) {
    object@colData
})


# official function to add sample metadata
# setMethod("addSampleData", "factR", function(object, data, colname = NULL, set = NULL) {
#
#     if(any(class(data) %in% c("data.frame", "tbl", "tbl_df"))){
#         # check for dimensions
#         if(nrow(data) != nrow(object@colData)){
#             rlang::abort(sprintf("Replacement length (%s) not equal to data length (%s)",
#                                  nrow(data), nrow(object@colData)))
#         }
#
#         # check for rownames
#         if(all(rownames(data) %in% rownames(object@colData))){
#             featurelvls <- rownames(object@colData)
#             object@colData <- cbind(object@colData,
#                                                 data[featurelvls,])
#         } else {
#             object@colData <- cbind(object@colData,
#                                                 data)
#         }
#     } else {
#         if(is.null(colname)){
#             colname <- "sampleData"
#         }
#         if(colname %in% colnames(object@colData)){
#             colname <- paste0(colname, ".1")
#         }
#         object@colData[[colname]] <- data
#     }
#     object
# })

























