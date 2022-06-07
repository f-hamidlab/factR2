#' @include generics.R

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Preview #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# show and summary prints a summary of the object
setMethod("show", "factR", function(object){
    cat(sprintf("class: factRObject [version %s]\n", object@version))
    cat(sprintf("# transcriptome: "))
    ngenes <- length(unique(object@transcriptome$gene_id))
    ntxs <- length(unique(object@transcriptome$transcript_id))
    nnovel <- sum(featureData(object)$novel == "yes")
    ncds <- sum(featureData(object)$cds == "yes")
    cat(sprintf("%s genes; ", ngenes))
    cat(sprintf("%s transcripts [%s novel]; ", ntxs, nnovel))
    cat(sprintf("%s coding transcripts \n", ncds))
    cat(sprintf("# active set: %s\n", object@active.set))

    nsamples <- ncol(object@colData)
    samples <- colnames(object@colData)
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
    utils::head(featureData(x), n = n)
})
setMethod("tail", "factR", function(x, n = 6L){
    utils::tail(featureData(x), n = n)
})

setMethod("dim", "factR", function(x){
    utils::tail(featureData(x), n = n)
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Sets Data #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## switch sets on the fly
setMethod("[[", "factR", function(x, i, j){
    if(i %in% names(x@sets)){
        x@active.set <- i
    } else if(i %in% 1:3){
        x@active.set <- names(x@sets)[i]
    } else {
        rlang::abort("Incorrect set name or index")
    }
    x
})

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

setMethod("rangesData", "factR", function(object, ..., set = NULL) {
    if(is.null(set)){
        set <- slot(object, "active.set")
    } else if(!set %in% listSets(object)){
        rlang::warn("Incorrect set name or index, using active set")
        set <- slot(object, "active.set")
    }
    gtf <- methods::slot(object, "transcriptome")
    txs <- tryCatch(.getFeat(object, ...),
                    error = function(e) rlang::abort("Feature not found"))

    if(set == "all"){
        return(gtf[gtf$transcript_id %in% txs])
    } else {
        return(gtf[gtf$transcript_id %in% txs & gtf$type %in% set])
    }

})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Feature Data #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# feature preview with option to subset data
setMethod("featureData", "factR", function(object, ..., set = NULL) {
    if(is.null(set)){
        set <- slot(object, "active.set")
    } else if(!set %in% listSets(object)){
        rlang::warn("Incorrect set name or index, using active set")
        set <- slot(object, "active.set")
    }

    dat <- slot(object@sets[[set]], "rowData")
    if(set %in% c("gene", "AS")){
        feat <- tryCatch(.getFeat(object, ..., out = "gene_id"),
                        error = function(e) rlang::abort("Feature not found"))
        return(dat[dat$gene_id %in% feat,])

    } else{
        feat <- tryCatch(.getFeat(object, ...),
                        error = function(e) rlang::abort("Feature not found"))
        return(dat[dat$transcript_id %in% feat,])
    }

})

# access and modify columns in featureData
setMethod("featureData$", "factR", function(object,name) {
    print('check')
    object@sets[[object@active.set]]@rowData[[name]]
})
setMethod("featureData<-", "factR", function(object,value) {
    object@sets[[object@active.set]]@rowData <- value
    object
})
setMethod("featureData$<-", "factR", function(object, name, value) {
    object@sets[[object@active.set]]@rowData[[name]] <- value
    object
})

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
setMethod("sampleData", "factR", function(object) {
    object@colData
})

# access and modify columns in sampleData
setMethod("sampleData$", "factR", function(object,name) {
    object@colData[[name]]
})
setMethod("sampleData<-", "factR", function(object,value) {
    object@colData <- value
    object
})
setMethod("sampleData$<-", "factR", function(object, name, value) {
    object@colData[[name]] <- value
    object
})

# official function to add sample metadata
setMethod("addSampleData", "factR", function(object, data, colname = NULL, set = NULL) {

    if(any(class(data) %in% c("data.frame", "tbl", "tbl_df"))){
        # check for dimensions
        if(nrow(data) != nrow(object@colData)){
            rlang::abort(sprintf("Replacement length (%s) not equal to data length (%s)",
                                 nrow(data), nrow(object@colData)))
        }

        # check for rownames
        if(all(rownames(data) %in% rownames(object@colData))){
            featurelvls <- rownames(object@colData)
            object@colData <- cbind(object@colData,
                                                data[featurelvls,])
        } else {
            object@colData <- cbind(object@colData,
                                                data)
        }
    } else {
        if(is.null(colname)){
            colname <- "sampleData"
        }
        if(colname %in% colnames(object@colData)){
            colname <- paste0(colname, ".1")
        }
        object@colData[[colname]] <- data
    }
    object
})

























