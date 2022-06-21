.updatefactR <- function(object){

    # perform sample check if colData was modifieed
    if("old.names" %in% colnames(object@colData)){
        old.n <- object@colData$old.names
        new.n <- rownames(object@colData)
        # correct counts data
        for(i in listSets(object)){

            if(!i %in% "AS"){
                object@sets[[i]]@counts <- object@sets[[i]]@counts[,old.n]
                colnames(object@sets[[i]]@counts) <- new.n
            }
            object@sets[[i]]@data <- object@sets[[i]]@data[,old.n]
            colnames(object@sets[[i]]@data) <- new.n
        }
        object@colData$old.names <- NULL
    }
    
    
    return(object)

    #all(sapply(list(colnames(counts())), FUN = identical, A))
}
