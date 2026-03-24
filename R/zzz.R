#' @importFrom rJava .jpackage .jcall
#' @importFrom rjd3toolkit get_java_version minimal_java_version
.onAttach <- function(libname, pkgname) {
    if (rjd3toolkit::get_java_version() < rjd3toolkit::minimal_java_version) {
        packageStartupMessage(sprintf(
            "Your java version is %s. %s or higher is needed.",
            rjd3toolkit::get_java_version(),
            rjd3toolkit::minimal_java_version
        ))
    }
}
.onLoad <- function(libname, pkgname) {
    result <- .jpackage(pkgname, lib.loc = libname)
    if (!result) {
        stop("Loading java packages failed", call. = FALSE)
    }

    # reload extractors

    if (rjd3toolkit::get_java_version() >= rjd3toolkit::minimal_java_version) {
        # reload providers
        try({
            rJava::.jcall(
                obj = "jdplus/toolkit/base/api/information/InformationExtractors",
                "V",
                "reloadExtractors"
            )
        })
    }
}
