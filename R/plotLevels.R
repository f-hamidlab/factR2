.plotLvls <- function(object, ..., set = NULL, group.by = NULL, slot = "data",
                        ncol=1){

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
    # check slot input
    if(!slot %in% c("counts","data")){
        slot <- "data"
        .msgwarn(stringr::str_glue("Data `{slot}` not found, using data slot"))
    }
    if(set=="AS" & slot=="counts"){
        slot <- "data"
    }
    # check group.by input
    if(is.null(group.by)){
        group.by <- object@active.ident
    } else if(!group.by %in% colnames(object@colData)){
        rlang::abort("grouping not found")
    }

    dat <- slot(object@sets[[set]], slot)  # get actual data

    # filter data by features
    settoout <- c("gene"="gene_id", "transcript"="transcript_id", "AS"="AS_id")
    feat <- .getFeat(object, ..., out = settoout[[set]])
    if(length(feat)==1){
        dat <- as.data.frame(dat[feat,])
        colnames(dat) <- feat
    } else {
        dat <- t(dat[feat,])
    }



    # join samples to levels data
    samples <- object@colData
    annodat <- cbind(samples, dat)

    # # pivot data longer
    # annodat.pivot <- annodat %>%
    #     tidyr::pivot_longer(cols=feat, names_to = "features",
    #                         values_to = "levels")

    # run main plot
    ylab <- c("gene"="Expression", "transcript"="Expression", "AS"="PSI")
    plot <- BiocGenerics::do.call(
        patchwork::wrap_plots,
        lapply(feat, function(y){
            annodat %>%
                ggplot(aes_string(x=group.by, y=y, fill=group.by)) +
                geom_violin() +
                geom_jitter(alpha = 0.3, fill = "grey20")+
                ggplot2::ggtitle(y) +
                ylab(ylab[set]) +
                theme_minimal() +
                theme(legend.position = "none")
        }))
    plot + patchwork::plot_layout(ncol = ncol)



}
