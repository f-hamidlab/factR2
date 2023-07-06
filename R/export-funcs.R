export.factR <- function(object, con, data){
    # check for data input
    if(!data %in% c("gtf", "gene", "transcript", "AS")){
        rlang::abort("`data` to output does not exists")
    }

    # check for con
    if(is_dir(con)){
        dir.create(con, showWarnings = FALSE)  # create dir in case
        if(data == "gtf"){
            con <- file.path(con, "factR2.gtf")
        } else{
            con <- file.path(con, paste0("factR2_",data,".tsv"))
        }
    }


    if(data == "gtf"){
        rtracklayer::export(object@transcriptome, con)
    } else{
        write.table(object[[data]], file = con, sep = "\t", quote = FALSE)
    }
}
