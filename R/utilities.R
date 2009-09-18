.STRAND_LEVELS <- levels(strand())

.show_classname <-
    function(x) cat("class: ", class(x), "\n", sep="")

.uunlist <- function(x) unlist(x, use.names=FALSE)
