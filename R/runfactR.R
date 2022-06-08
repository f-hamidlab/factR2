setMethod("runfactR", "factR", function(object, verbose = FALSE) {
    object <- prepTranscriptome(object, verbose)
    object <- buildCDS(object, verbose)
    object <- predictNMD(object, verbose)
    object <- getAAsequence(object, verbose)
    object <- testASNMDevents(object)
    object
})

















