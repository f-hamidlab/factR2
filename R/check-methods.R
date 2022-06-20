checkfactR.factR <- function(object){

    # perform sample check if there is sample info
    if(nrow(object@colData) >0){
        ## check if all samples in coldata matches with counts data
        samp.check.out <- do.call(c, lapply(listSets(object), function(x)
            !identical(colnames(counts(object, set = x)),
                      rownames(samples(object)))))

        if(any(samp.check.out)){
            sets.to.chg <- listSets(object)[samp.check.out]
            for(x in sets.to.chg){
                # t
            }

        }
    }

    #all(sapply(list(colnames(counts())), FUN = identical, A))
}
