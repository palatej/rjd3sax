#' @importFrom rjd3toolkit sa_decomposition
NULL

#' @export
sa_decomposition.JD3_LTDARIMA <- function(x, ...) {
  if (is.null(x)) {
    return(NULL)
  }
  x<-x$final$decomposition
  if (is.null(x)) {
    return(NULL)
  }
  return (.sa(x$finals))
}

#' @export
sa_decomposition.JD3_LTDARIMA_TD <- function(x, ...) {
  return (.sa(x))
}

#' @export
sa_decomposition.JD3_LTDARIMA_FIXED <- function(x, ...) {
  return (.sa(x))
}

.sa<-function(x){
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
