.build_path <- 
    function(path)
{
    if (.Platform$OS.type == "windows") {
        path <- normalizePath(path)
        if (grepl(' ', path, fixed=TRUE))
            path <- utils::shortPathName(path)
        path <- gsub("\\\\", "/", path)
    }
    path
}

.pkgconfig <-
    function()
{
    path <- if (.Platform$OS.type == "windows") {
        system.file(package="Rsamtools", "usrlib", .Platform[["r_arch"]],
                    "Rsamtools.mk", mustWork=TRUE)
    } else {
        system.file(package="Rsamtools", "usrlib", .Platform[["r_arch"]],
                    mustWork=TRUE)
    }
    .build_path(path)
}
