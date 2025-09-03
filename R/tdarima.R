#' @include utils.R
NULL

LTDARIMA<-'JD3_LtdArima'



#' Title
#'
#' @param data
#' @param period
#' @param th
#' @param bth
#' @param se
#'
#' @returns
#' @export
#'
#' @examples
tdairline_decomposition<-function(data, th, bth, se=FALSE){
  if (! is.ts(data)) stop("data should be a time series (ts)")
  jmatrix<-.jcall('jdplus/advancedsa/base/r/TimeVaryingArimaModels', 'Ljdplus/toolkit/base/api/math/matrices/Matrix;', 'airlineDecomposition',
            as.numeric(data), as.integer(frequency(data)), as.numeric(th), as.numeric(bth), as.logical(se))
  return (rjd3toolkit::.jd2r_matrix(jmatrix))
}

#' Title
#'
#' @param data
#' @param mean
#' @param x
#' @param regular
#' @param seasonal
#' @param fixed_phi
#' @param fixed_bphi
#' @param fixed_theta
#' @param fixed_btheta
#' @param fixed_var
#' @param eps
#'
#' @returns
#' @export
#'
#' @examples
ltdarima_estimation<-function(data, mean=FALSE, X=NULL, regular, seasonal, fixed_phi=TRUE, fixed_bphi=TRUE,
                                fixed_theta=FALSE, fixed_btheta=FALSE, fixed_var=TRUE, eps=1e-7){
  if (! is.ts(data)) stop("data should be a time series (ts)")
  jrslt<-.jcall('jdplus/advancedsa/base/r/TimeVaryingArimaModels', 'Ljdplus/advancedsa/base/core/tdarima/LtdArimaResults;', 'estimate',
                  as.numeric(data), as.integer(frequency(data)), as.logical(mean), rjd3toolkit::.r2jd_matrix(X), as.integer(regular), as.integer(seasonal)
                  , as.logical(fixed_phi), as.logical(fixed_bphi), as.logical(fixed_theta), as.logical(fixed_btheta), as.logical(fixed_var), as.numeric(eps))
  return(rjd3toolkit::.jd3_object(jrslt, result=TRUE))

}




