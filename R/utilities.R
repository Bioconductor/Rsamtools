.STRAND_LEVELS <- levels(strand())

.show_classname <-
    function(x) cat("class: ", class(x), "\n", sep="")

.normalizePath <-
    function(path)  normalizePath(path.expand(path))

.uunlist <-
    function(x) unlist(x, use.names=FALSE)
