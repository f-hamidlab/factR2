select.factR <- function(object, ...){
    samples <- object@colData
    samples$old.names <- rownames(samples)
    samples <- as.data.frame(t(samples))

    # perform col selection
    samples <- as.data.frame(t(dplyr::select(samples, ...)))

    # correct counts data
    for(i in listSets(object)){

        if(!i %in% "AS"){
            object@sets[[i]]@counts <- object@sets[[i]]@counts[,samples$old.names]
            colnames(object@sets[[i]]@counts) <- rownames(samples)
        }
        object@sets[[i]]@data <- object@sets[[i]]@data[,samples$old.names]
        colnames(object@sets[[i]]@data) <- rownames(samples)
    }

    return(object)
}
