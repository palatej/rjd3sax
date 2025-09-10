#' @include utils.R
NULL

LTDARIMA<-'JD3_LTDARIMA_RSLTS'
LTDARIMA_LL<-'JD3_LTDARIMA_LIKELIHOOD'



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

.jestimation<-function(data, mean=FALSE, X=NULL, regular, seasonal, fixed_phi=TRUE, fixed_bphi=TRUE,
                              fixed_theta=FALSE, fixed_btheta=FALSE, fixed_var=TRUE, eps=1e-7, parametrization="mean_delta"){
  if (! is.ts(data)) stop("data should be a time series (ts)")
  jrslt<-.jcall('jdplus/advancedsa/base/r/TimeVaryingArimaModels', 'Ljdplus/advancedsa/base/core/tdarima/LtdArimaResults;', 'estimate',
                as.numeric(data), as.integer(frequency(data)), as.logical(mean), rjd3toolkit::.r2jd_matrix(X), as.integer(regular), as.integer(seasonal)
                , as.logical(fixed_phi), as.logical(fixed_bphi), as.logical(fixed_theta), as.logical(fixed_btheta), as.logical(fixed_var)
                , as.numeric(eps), as.character(parametrization))
  return (jrslt)
}

#' Title
#'
#' @param data
#' @param regular
#' @param seasonal
#' @param p0
#' @param p1
#' @param var1
#' @param se
#'
#' @returns
#' @export
#'
#' @examples
tdarima_decomposition<-function(data, regular, seasonal, p0, p1, var1=1, se=FALSE){
  if (! is.ts(data)) stop("data should be a time series (ts)")
  m<-length(p0)
  n<-length(data)
  if (length(p1) != m || n == 0) stop("invalid parameters")
  p<-NULL
  if (m == 1){
    p<-matrix(.linear(p0, p1, n), nrow = 1, ncols=n)
  }else if (m>1){
    p<-.linear(p0, p1, n)
  }
  if (var1 != 1){
    var<-matrix(.linear(1, var1, n), nrow = 1, ncols=n)
    p<-rbind(p, var)
  }
  if (is.null(p)){
    jp<-.jnull("jdplus/toolkit/base/api/math/matrices/Matrix")
  }
  else{
    jp<-.r2jd_matrix(p)
  }

  jmatrix<-.jcall('jdplus/advancedsa/base/r/TimeVaryingArimaModels', 'Ljdplus/toolkit/base/api/math/matrices/Matrix;', 'arimaDecomposition',
                  as.numeric(data), as.integer(frequency(data)), as.integer(regular), as.integer(seasonal), as.logical(var1 != 1), jp, as.logical(se))
  return (rjd3toolkit::.jd2r_matrix(jmatrix))
}

.linear<-function(p0, p1, n){
  delta<-(p1-p0)/(n-1)
  sapply(0:(n-1), function(i){p0+i*delta})
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
                                fixed_theta=FALSE, fixed_btheta=FALSE, fixed_var=TRUE, eps=1e-7, parametrization=c("mean_delta", "start_end")){

  parametrization=match.arg(parametrization)
  jrslt<-.jestimation(data, mean, X, regular, seasonal, fixed_phi, fixed_bphi, fixed_theta, fixed_btheta, fixed_var, eps, parametrization)

  freq<-frequency(data)
  start=start(data)

  regs0<-.proc_vector(jrslt, "regression.effect0")
  regs1<-.proc_vector(jrslt, "regression.effect1")
  tsregs0 <-NULL;
  tsregs1<-NULL;

  if ( ! is.null(regs0)){
    tsregs0<-ts(data=regs0, frequency = freq, start = start)
    tsregs1<-ts(data=regs1, frequency = freq, start = start)
  }

  # initial (fixed) model
  initial=list(
    model=list(
      period=freq,
      regular=regular,
      seasonal=seasonal,
      parameters=.proc_vector(jrslt, "model.pfixed"),
      covariance=.proc_matrix(jrslt, "model.pfixed_cov")),
    likelihood=structure(.proc_likelihood(jrslt, "ll0."), class=LTDARIMA_LL),
    regression=list(
      coefficients=.proc_vector(jrslt, "regression.c0"),
      covariance=.proc_matrix(jrslt, "regression.cov0"),
      linearized=ts(data=.proc_vector(jrslt, "regression.y_lin0"), frequency = freq, start = start),
      regression_effect=tsregs0),
    residuals=list(
      res=.proc_vector(jrslt, "res0.res"),
      type=.proc_data(jrslt, "res0.type"),
      df=.proc_int(jrslt, "res0.df"),
      dfc=.proc_int(jrslt, "res0.dfc"),
      mean=.proc_test(jrslt, "res0.mean"),
      skewness=.proc_test(jrslt, "res0.skewness"),
      kurtosis=.proc_test(jrslt, "res0.kurtosis"),
      normality=.proc_test(jrslt, "res0.doornikhansen"),
      ljung_box=.proc_test(jrslt, "res0.lb"),
      seasonal_ljung_box=.proc_test(jrslt, "res0.seaslb"),
      nruns=.proc_test(jrslt, "res0.nruns"),
      lruns=.proc_test(jrslt, "res0.lruns"),
      nudruns=.proc_test(jrslt, "res0.nudruns"),
      ludruns=.proc_test(jrslt, "res0.ludruns")
    )

  )

  # final model
  final=list(
    model=list(
      period=freq,
      regular=regular,
      seasonal=seasonal,
      parameters=.proc_vector(jrslt, "model.pall"),
      parima_0=.proc_vector(jrslt, "model.p0"),
      parima_1=.proc_vector(jrslt, "model.p1"),
      parima_mean=.proc_vector(jrslt, "model.pmean"),
      parima_delta=.proc_vector(jrslt, "model.pdelta"),
      covariance=.proc_matrix(jrslt, "model.pall_cov"),
      scores=.proc_vector(jrslt, "regression.ml.score1"),
      information=.proc_matrix(jrslt, "regression.ml.information1")),
    likelihood=structure(rjd3toolkit::.proc_likelihood(jrslt, "ll1."), class=LTDARIMA_LL),
    regression=list(
      coefficients=.proc_vector(jrslt, "regression.c1"),
      covariance=.proc_matrix(jrslt, "regression.cov1"),
      linearized=ts(data=.proc_vector(jrslt, "regression.y_lin1"), frequency = freq, start = start),
      regression_effect=tsregs1
    ),
    residuals=list(
      res=.proc_vector(jrslt, "res1.res"),
      type=.proc_data(jrslt, "res1.type"),
      df=.proc_int(jrslt, "res1.df"),
      dfc=.proc_int(jrslt, "res1.dfc"),
      mean=.proc_test(jrslt, "res1.mean"),
      skewness=.proc_test(jrslt, "res1.skewness"),
      kurtosis=.proc_test(jrslt, "res1.kurtosis"),
      normality=.proc_test(jrslt, "res1.doornikhansen"),
      ljung_box=.proc_test(jrslt, "res1.lb"),
      seasonal_ljung_box=.proc_test(jrslt, "res1.seaslb"),
      nruns=.proc_test(jrslt, "res1.nruns"),
      lruns=.proc_test(jrslt, "res1.lruns"),
      nudruns=.proc_test(jrslt, "res1.nudruns"),
      ludruns=.proc_test(jrslt, "res1.ludruns")
    )
  )

  return(structure(list(initial=initial, final=final), class=LTDARIMA))
}




