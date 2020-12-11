library("tidyverse")

load_bioconductor_package <- function(path_to_bioc_pkg, pkg_name) {

    lib <- str_remove(path_to_bioc_pkg, pkg_name)

    # ensure that dependencies of the pkg are also found at same location
    .libPaths( c( lib , .libPaths() ) )

    library(pkg_name, character.only = TRUE)

    print(str_c("loaded package ", pkg_name))

    # ensure that library() calls outside this function don't go looking in the
    # location needed here
    .libPaths( .libPaths()[-1] )
}