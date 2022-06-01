

.DollarNames.factR <- function(x, pattern = "")
    grep(pattern, colnames(as.data.frame(x@custom)), value=TRUE)
setMethod("$", signature("factR"), function (x, name) {
    # check features
    x <- as.data.frame(x@custom)
    x[[name]]

})

setMethod("$<-", "factR", function(x, name, value){
    mcols(x@custom)[[name]] <- value
    x
})



