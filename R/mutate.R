mutate.factR <- function(object, ..., data = "samples"){

    # check for type of data to modify
    if(!data %in% c("gene", "transcript", "AS", "samples")){
        rlang::abort("Incorrect input for `data`")
    }

    # get data
    if(data == "samples"){
        df <- object@colData
    } else {
        df <- object[[data]]
    }

    # store rownames and perform
    names <- rownames(df)
    df <- dplyr::mutate(df, ...)
    rownames(df) <- names  # restore rownames


    # update data and return object
    if(data == "samples"){
        object@colData <- df
    } else {
        object@sets[[data]]@rowData <- df
    }
    return(object)
}
