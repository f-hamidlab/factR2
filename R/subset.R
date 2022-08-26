#TODO: Check or refine this. Maybe change name
select.factR <- function(object, ...){
    samples <- object@colData
    samples$old.names <- rownames(samples)
    samples <- as.data.frame(t(samples))

    samples <- tryCatch(
        samples[,...],
        error = function(e)
            dplyr::select(samples, ...))

    # perform col selection
    object@colData <- as.data.frame(t(samples))

    return(.updatefactR(object))
}
