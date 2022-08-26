#TODO: Check or refine this. Maybe change name
mutate.factR <- function(object, meta = "samples", data=NULL, ...){

    # check for type of meta to modify
    if(!meta %in% c("gene", "transcript", "AS", "samples")){
        rlang::abort("Incorrect input for `meta`")
    }

    # get meta
    if(meta == "samples"){
        df <- object@colData
    } else {
        df <- object[[meta]]
    }

    # Add dataframe to metadata if given
    if(!is.null(data)){

        # check for rownames
        if(is.null(rownames(data))){
            if(nrow(data) != nrow(df)){
                rlang::abort("Length mismatch")
            }
        } else if(!all(rownames(data) %in% rownames(df))){
            # create new data.frame with rownames from df
            tokeep <- rownames(data)[rownames(data) %in% rownames(df)]
            newdata <- data.frame(row.names = rownames(df))
            newdata[tokeep,colnames(data)] <- data[tokeep,colnames(data)]
            data <- newdata
        } else{
            data <- data[rownames(df),]
        }

        df <- dplyr::mutate(df, data)
    } else{
        # do dplyr mutate style transformation

        df <- dplyr::mutate(df, ...)
    }


    # update meta and return object
    if(meta == "samples"){
        object@colData <- df
    } else {
        object@sets[[meta]]@rowData <- df
    }
    return(object)
}
