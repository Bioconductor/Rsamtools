.yieldReduce_iterate <-
    function(X, MAP, REDUCE, DONE, ..., init)
{
    result <- if (missing(init)) {
        value <- MAP(X, ...)
        if (DONE(value))
            return(list())
        value
    } else
        init

    repeat {
        value <- MAP(X, ...)
        if (DONE(value))
            break
        result <- REDUCE(result, value)
    }
    result
}

.yieldReduce_all <-
    function(X, MAP, REDUCE, DONE, ...)
{
    N_GROW <- 100L
    n <- 0
    result <- vector("list", n)
    i <- 0L
    repeat {
        value <- MAP(X, ...)
        if (DONE(value))
            break
        i <- i + 1L
        if (i > n) {
            n <- n + N_GROW
            length(result) <- n
        }
        result[[i]] <- value
    }
    length(result) <- i
    REDUCE(result)
}

yieldReduce <-
    function(X, MAP, REDUCE, DONE, ..., init, ITERATE=TRUE)
{
    if (missing(REDUCE)) 
        REDUCE <- if (ITERATE) c else identity
    if (missing(DONE))
        DONE <- function(VALUE) length(VALUE) == 0L
    if (!ITERATE && !missing(init))
        warning("'init' ignored when ITERATE==FALSE")

    if (!isOpen(X)) {
        open(X)
        on.exit(close(X))
    }
    if (ITERATE)
        .yieldReduce_iterate(X, MAP, REDUCE, DONE, ..., init=init)
    else
        .yieldReduce_all(X, MAP, REDUCE, DONE, ...)
}    
