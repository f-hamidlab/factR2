.compLvls <- function(object, A, B = NULL, ident = NULL,
                      col_prefix = NULL,
                      set = NULL, return_df = FALSE){

    # get levels data based on set
    if(is.null(set)){
        set <- object@active.set
    } else{
        # check if given set is part of object
        if(!set %in% names(object@sets)){
            set <- object@active.set
            .msgwarn(stringr::str_glue("Set `{set}` not found, using active.set"))
        }
    }

    # check ident input
    if(is.null(ident)){
        ident <- object@active.ident
    }

    # check comparisons
    if(!all(A %in% samples(object)[[ident]])){
        rlang::abort("Sample group not found in current ident")
    } else {
        compgroupings <- ifelse(samples(object)[[ident]] %in% A, 1, 0)
        if(!is.null(B)){
            compgroupings <- ifelse(samples(object)[[ident]] %in% B, 2, compgroupings)
        } else {
            compgroupings <- ifelse(compgroupings==0, 2, compgroupings)
        }
    }
    # check col_prefix name
    if(is.null(col_prefix)){
        if(is.null(B)){
            B <- "rest"
        }
        name_A <- paste0(A, collapse = ":")
        name_B <- paste0(B, collapse = ":")

        col_prefix <- stringr::str_glue(
            "{name_A}_vs_{name_B}")
    }

    # run Double GLM model fitting
    dge <- suppressWarnings(DoubleExpSeq::DBGLM1(object@sets$AS@misc$inc, object@sets$AS@misc$total,
                                groups = compgroupings,
                                contrast = c(1,2)))

    out <- dge$All %>%
        as.data.frame() %>%
        dplyr::mutate(deltaPSI = MLE_1-MLE_2) %>%
        dplyr::select(deltaPSI, pval = pVal, padj = Adj.pVal) %>%
        dplyr::mutate(across(.fns = signif, digits=3)) %>%
        dplyr::rename_with(~stringr::str_glue("{col_prefix}.{.x}"))
    if(return_df){
        return(out)
    } else {
        ## append feature dataframe
        object@sets$AS@rowData[colnames(out)] <-  out[rownames(object@sets$AS@rowData),]
        return(object)
    }

}
