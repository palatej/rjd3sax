#' @importFrom rjd3toolkit sa_decomposition
#' @export
sa_decomposition.JD3_LTDARIMA_RSLTS <- function(x, ...) {
  if (is.null(x)) {
    return(NULL)
  }
  x<-x$decomposition
  if (is.null(x)) {
    return(NULL)
  }
  decomp<-x$finals
  if (is.null(decomp)) decomp<-x$components
  return(rjd3toolkit::sadecomposition(
    decomp$y,
    decomp$sa,
    decomp$trend,
    decomp$seas,
    decomp$irregular,
    FALSE
  ))
}


#' @export
sa_decomposition.JD3_LTDARIMA_RSLTS <- function(x, ...) {
  return(rjd3toolkit::sa_decomposition(x$result))
}
