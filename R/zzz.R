#' @include utils.R tdarima_plot.R tdarima_print.R
#' @importFrom rjd3jars check_java_version
#' @importFrom rjd3toolkit reload_dictionaries
NULL

.onLoad <- function(libname, pkgname) {
  if (!requireNamespace("rjd3xjars", quietly = TRUE)) {
    stop("Loading {rjd3xjars} failed", call. = FALSE)
  }
  if (!requireNamespace("rjd3toolkit", quietly = TRUE)) {
    stop("Loading {rjd3toolkit} failed", call. = FALSE)
  }
  # Loading Java class
  jar_dir <- file.path(libname, pkgname, "inst", "java")
  jars_inst <- list.files(
    jar_dir,
    pattern = "\\.jar$",
    full.names = TRUE,
    all.files = TRUE
  )
  result <- rJava::.jpackage(
    pkgname,
    lib.loc = libname,
    morePaths = jars_inst
  )
  if (!result)
    stop("Loading java packages failed")


  assign("advancedsa", list(), rjd3toolkit::.jd3_env)

  # reload extractors
  if (rjd3jars::check_java_version()) rjd3jars::reload_dictionaries()

}
