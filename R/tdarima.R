#' @include utils.R
NULL

LTDARIMA<-'JD3_LTDARIMA_RSLTS'
LTDARIMA_LL<-'JD3_LTDARIMA_LIKELIHOOD'



#' Estimation by means of the Kalman smoother of a time-dependent canonical decomposition
#' of airline models
#'
#' @param data The series (a ts).
#' @param th The regular moving average parameter (same length as data)
#' @param bth The seasonal moving average parameter (same length as data)
#' @param se Specifies if the standard deviations of the components are computed
#'
#' @returns A mts object is returned, with the trend, the seasonal, the irregular components and - if requested - their standard deviations
#' @export
#'
#' @examples
#' s<-rjd3toolkit::ABS$X0.2.09.10.M
#' th<-seq(.2,-.9, -1.1/(length(s)-1))
#' bth<-seq(-.9,-.2, .7/(length(s)-1))
#' sa<-tdairline_decomposition(s, th, bth, se=TRUE)
#' ts.plot(ts.union(s, sa[,1], sa[,1]+sa[,3]), col=c('gray', 'blue', 'red'))
#' ts.plot(sa[,c(4,5,6)], col=c('red', 'blue', 'magenta'))
tdairline_decomposition<-function(data, th, bth, se=FALSE){
  if (! is.ts(data)) stop("data should be a time series (ts)")
  jmatrix<-.jcall('jdplus/advancedsa/base/r/TimeVaryingArimaModels', 'Ljdplus/toolkit/base/api/math/matrices/Matrix;', 'airlineDecomposition',
            as.numeric(data), as.integer(frequency(data)), as.numeric(th), as.numeric(bth), as.logical(se))
  z<-rjd3toolkit::.jd2r_matrix(jmatrix)
  names<-c("trend", "seasonal", "irregular")
  if (se){
    names<-c(names, paste0(names, "-stdev"))
  }
  colnames(z)<-names
  all<-ts(z, frequency=frequency(data), start=start(data))
  return (all)
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

#' Estimation by means of the Kalman smoother of the canonical decomposition
#' of linear time-dependent SARIMA models
#'
#' @param data The series (a ts).
#' @param regular The orders of the regular part of the SARIMA models (p,d,q)
#' @param seasonal The orders of the seasonal part of the SARIMA models (bp,bd,bq). NULL is (0,0,0)
#' @param p0 The arima parameters at the beginning of the series
#' @param p1 The arima parameters at the end of the series
#' @param var1 The (unscaled) variance of the innovations at the end of the series (1 at the beginning of the series)
#' @param se Specifies if the standard deviations of the components are computed
#'
#' @returns A mts object is returned, with the trend, the seasonal, the irregular components and - if requested - their standard deviations
#' @export
#'
#' @details The arima parameters are organized in the following order: phi, bphi, theta, btheta
#'
#' @examples
#' s<- rjd3toolkit::Retail$RetailSalesTotal
#'
#' q<-rjd3sax::ltdarima_estimation(s, regular=c(0,1,1), seasonal=c(0,1,1),
#'  fixed_var=FALSE, eps=1e-15, parametrization = "mean_delta")
#' sa<-rjd3sax::ltdarima_decomposition(s, c(0,1,1), c(0,1,1),
#' q$final$model$parima_0, q$final$model$parima_1, var1=q$final$model$parameters[5], se=TRUE)
#' ts.plot(ts.union(s, sa[,1], sa[,1]+sa[,3]), col=c('gray', 'blue', 'magenta'))
ltdarima_decomposition<-function(data, regular, seasonal, p0, p1, var1=1, se=FALSE){
  if (! is.ts(data)) stop("data should be a time series (ts)")
  m<-length(p0)
  n<-length(data)
  if (length(p1) != m || n == 0) stop("invalid parameters")
  p<-NULL
  if (m == 1){
    p<-matrix(.linear(p0, p1, n), nrow = 1, ncol=n)
  }else if (m>1){
    p<-.linear(p0, p1, n)
  }
  if (var1 != 1){
    var<-matrix(.linear(1, var1, n), nrow = 1, ncol=n)
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
  z<-rjd3toolkit::.jd2r_matrix(jmatrix)
  ncmps<-dim(z)[2]
  if (se)    ncmps<-ncmps/2
  if (ncmps == 3)  names<-c("trend", "seasonal", "irregular")else names <-c("trend", "irregular")
  if (se){
    names<-c(names, paste0(names, "-stdev"))
  }
  colnames(z)<-names
  all<-ts(z, frequency=frequency(data), start=start(data))
  return (all)
}

.pcount<-function(regular, seasonal){
  n<-regular[1]+regular[3]
  if (! is.null(seasonal))
    n=n+seasonal[1]+seasonal[3]
  return (n)
}

#' Estimation by means of the Kalman smoother of a time-dependent canonical decomposition (generic solution)
#'
#' @param data The series (a ts).
#' @param regular The orders of the regular part of the SARIMA models (p,d,q)
#' @param seasonal The orders of the seasonal part of the SARIMA models (bp,bd,bq). NULL is (0,0,0)
#' @param parameters The time dependent parameters. The matrix contains for each column
#' all the parameters of the model, with the columns corresponding to an observation. So,
#' the number of columns is identical to the length of the series and the number of rows corresponds
#' to the arima parameters (regular + seasonal), with eventually the innovation variances
#'
#' @param se Specifies if the standard deviations of the components are computed
#'
#' @returns A mts object is returned, with the trend, the seasonal, the irregular components and - if requested - their standard deviations
#' @export
#'
#' @details The arima parameters are organized in the following order: phi, bphi, theta, btheta
#' @export
#'
#' @examples
#' s<-rjd3toolkit::ABS$X0.2.09.10.M
#' n<-length(s)
#' r<-runif(n*2,-1,-.6)
#' p<-matrix(r, nrow=2, ncol=n)
#' p<-rbind(p, runif(n))
#' sa<-rjd3sax::tdarima_decomposition(s, c(0,1,1), c(0,1,1),parameter=p, se=TRUE)
#' ts.plot(ts.union(s, sa[,1], sa[,1]+sa[,3]), col=c('gray', 'blue', 'magenta'))
#' ts.plot(sa[,c(4,5,6)], col=c('red', 'blue', 'magenta'))
tdarima_decomposition<-function(data, regular, seasonal, parameters, se=FALSE){
  if (! is.ts(data)) stop("data should be a time series (ts)")
  if (! is.matrix(parameters)) stop("parameters should be a matrix")
  d=dim(parameters)
  n<-length(data)
  if (d[2] != n || n == 0) stop("invalid parameters")
  if (is.null(seasonal)) seasonal<-c(0,0,0)
  m=.pcount(regular, seasonal)
  if (d[1] != m && d[1] != m+1) stop("invalid parameters")
  var = d[1]>m
  jp<-.r2jd_matrix(parameters)
  jmatrix<-.jcall('jdplus/advancedsa/base/r/TimeVaryingArimaModels', 'Ljdplus/toolkit/base/api/math/matrices/Matrix;', 'arimaDecomposition',
                  as.numeric(data), as.integer(frequency(data)), as.integer(regular), as.integer(seasonal), as.logical(var), jp, as.logical(se))
  z<-rjd3toolkit::.jd2r_matrix(jmatrix)
  ncmps<-dim(z)[2]
  if (se)    ncmps<-ncmps/2
  if (ncmps == 3)  names<-c("trend", "seasonal", "irregular")else names <-c("trend", "irregular")
  if (se){
    names<-c(names, paste0(names, "-stdev"))
  }
  colnames(z)<-names
  all<-ts(z, frequency=frequency(data), start=start(data))
  return (all)
}

.linear<-function(p0, p1, n){
  delta<-(p1-p0)/(n-1)
  sapply(0:(n-1), function(i){p0+i*delta})
}



#' Estimation and decomposition of a regression model with linear time-dependent SARIMA noises
#'
#' @param data The series (a ts).
#' @param mean Mean correction
#' @param X The regression variables. A matrix with the same number of rows as the
#' data in the series
#' @param regular The orders of the regular part of the SARIMA models (p,d,q)
#' @param seasonal The orders of the seasonal part of the SARIMA models (bp,bd,bq). NULL is (0,0,0)
#' @param fixed_phi Indicates that the regular auto-regressive parameters are fixed or not.
#' @param fixed_bphi Indicates that the seasonal auto-regressive parameters are fixed or not.
#' @param fixed_theta Indicates that the regular moving average parameters are fixed or not.
#' @param fixed_btheta Indicates that the seasonal moving average parameters are fixed or not.
#' @param fixed_var Indicates that the variances of the innovations are fixed or not.
#' @param eps The precision of the optimization procedure
#' @param parametrization Type of the parametrization of the linear time-dependent parameters:
#' first and last parameters or average and variation of the parameters
#' @param decomposition Indicates that a decomposition of the series must be processed
#' @param regeffects In case of regression variables (excluding the mean correction),
#' the user must specify to which component the regression effects are associated
#' (1 for the trend, 2 for the seasonal and 3 for the irregular). If not specified, the
#' final components are not computed.
#' @param clean Cleans missing values at the beginning/end of the series. Should be
#' set to false if you want forecasts/backcasts
#'
#' @returns A JD3_LTDARIMA_RSLTS object
#' @export
#'
#' @examples
#' s<-rjd3toolkit::ABS$X0.2.09.10.M
#' td<-rjd3toolkit::td(s=s)
#' q<-rjd3sax::ltdarima_estimation(s, regular=c(0,1,1), seasonal=c(0,1,1), X=td,
#'  fixed_var = FALSE, eps=1e-15, parametrization = "mean_delta")


ltdarima_estimation<-function(data, mean=FALSE, X=NULL, regular=c(0,1,1), seasonal=c(0,1,1), fixed_phi=TRUE, fixed_bphi=TRUE,
                                fixed_theta=FALSE, fixed_btheta=FALSE, fixed_var=TRUE,
                              eps=1e-7, parametrization=c("mean_delta", "start_end"),
                              decomposition=TRUE, regeffects=NULL, clean=TRUE){

  parametrization=match.arg(parametrization)
  if (clean)
    data=clean_extremities(data)
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
  flin<-ts(data=.proc_vector(jrslt, "regression.y_lin1"), frequency = freq, start = start)
  parameters<-.proc_vector(jrslt, "model.pall")
  covariance<-.proc_matrix(jrslt, "model.pall_cov")
  pdetails<-.pdetails(parametrization=="mean_delta", regular, seasonal, fixed_phi, fixed_bphi, fixed_theta, fixed_btheta, fixed_var)
  pderived<-.pderived(parametrization=="mean_delta", parameters, covariance, pdetails$didx)
  m<-length(data)

  idx<-pdetails$sidx
  if (! is.null(idx)){
    np<-length(parameters)
    parameters[idx]<- parameters[idx]/(m-1)
    covariance[,idx]<-covariance[,idx]/(m-1)
    covariance[idx,]<-covariance[idx,]/(m-1)
  }else{
    ndp<-length(pderived$dp)
    i<-2
    while (i<= ndp){
      pderived$dp[i]<-pderived$dp[i]/(m-1)
      pderived$edp[i]<-pderived$edp[i]/(m-1)
      i<-i+2
    }
  }

  # final model
  final=list(
    model=list(
      type=parametrization,
      n=m,
      period=freq,
      regular=regular,
      seasonal=seasonal,
      parameters_names=pdetails$names,
      parameters=parameters,
      parameters_stde=suppressWarnings(sqrt(diag(covariance))),
      derived_parameters_names=pdetails$dnames,
      derived_parameters=pderived$dp,
      derived_parameters_stde=pderived$edp,
      parima_0=.proc_vector(jrslt, "model.p0"),
      parima_1=.proc_vector(jrslt, "model.p1"),
      parima_mean=.proc_vector(jrslt, "model.pmean"),
      parima_delta=.proc_vector(jrslt, "model.pdelta"),
      covariance=covariance,
      scores=.proc_vector(jrslt, "regression.ml.score1"),
      information=.proc_matrix(jrslt, "regression.ml.information1")),
    likelihood=structure(rjd3toolkit::.proc_likelihood(jrslt, "ll1."), class=LTDARIMA_LL),
    regression=list(
      coefficients=.proc_vector(jrslt, "regression.c1"),
      covariance=.proc_matrix(jrslt, "regression.cov1"),
      linearized=flin,
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

  sa<-NULL
  if (decomposition == TRUE){
    s<-data
    var<-1
    if (! fixed_var){
      fp<-final$model$parameters
      var<-fp[length(fp)]
    }
    cmps<-ltdarima_decomposition(flin, regular = regular, seasonal = seasonal,
                                 p0=final$model$parima_0, p1=final$model$parima_1, var1 = var, se=TRUE)
    components<-list(
      y=s,
      sa=cmps[,1]+cmps[,3],
      trend=cmps[,1],
      seas=cmps[,2],
      irregular=cmps[,3],
      trend_stdev=cmps[,4],
      seas_stdev=cmps[,5],
      irregular_sted=cmps[,6]
    )
    finals<-NULL
    if (!is.null(X) && ! is.null(regeffects)){

    }
    sa<-list(
      components=components,
      finals=finals
    )
  }

  return(structure(list(initial=initial, final=final, decomposition=sa), class=LTDARIMA))
}

.pdetails<-function(meandelta, regular, seasonal, fixed_phi, fixed_bphi,
                   fixed_theta, fixed_btheta, fixed_var){
  p = regular[1]
  d = regular[2]
  q = regular[3]
  bp = seasonal[1]
  bd = seasonal[2]
  bq = seasonal[3]
  names <- NULL
  dnames <- NULL
  didx <- NULL
  sidx <- NULL
  icur <- 1
    # first, p or p0 or pmean
  if (p > 0) {
    if (fixed_phi){
      names <- c(names, paste0("phi(", 1:p, ")"))
    }else{
      if (meandelta){
        names <- c(names, paste0("phi-mean(", 1:p, ")"))
      }else{
        names <- c(names, paste0("phi-start(", 1:p, ")"))
      }
      didx<-c(didx, icur:(icur+p-1))
    }
    icur<-icur+p
  }
  if (bp > 0) {
    if (fixed_bphi){
      names <- c(names, paste0("bphi(", 1:bp, ")"))
    }else{
      if (meandelta){
        names <- c(names, paste0("bphi-mean(", 1:bp, ")"))
      }else{
        names <- c(names, paste0("bphi-start(", 1:bp, ")"))
      }
      didx<-c(didx, icur:(icur+bp-1))
    }
    icur<-icur+bp
  }
  if (q > 0) {
    if (fixed_theta){
      names <- c(names, paste0("theta(", 1:q, ")"))
    }else{
      if (meandelta){
        names <- c(names, paste0("theta-mean(", 1:q, ")"))
      }else{
        names <- c(names, paste0("theta-start(", 1:q, ")"))
      }
      didx<-c(didx, icur:(icur+q-1))
    }
    icur<-icur+q
  }
  if (bq > 0) {
    if (fixed_btheta){
      names <- c(names, paste0("btheta(", 1:bq, ")"))
    }else{
      if (meandelta){
        names <- c(names, paste0("btheta-mean(", 1:bq, ")"))
      }else{
        names <- c(names, paste0("btheta-start(", 1:bq, ")"))
      }
      didx<-c(didx, icur:(icur+bq-1))
    }
    icur<-icur+bq
  }
  # then, p1 or pdelta + derived
  if (p > 0 && ! fixed_phi) {
    if (meandelta){
      names <- c(names, paste0("phi-delta(", 1:p, ")"))
      sidx<-c(sidx, icur:(icur+p-1))
      dnames <- c(dnames, paste0("phi-start(", 1:p, ")[derived]"), paste0("phi-end(", 1:p, ")[derived]"))
    }else{
      names <- c(names, paste0("phi-end(", 1:p, ")"))
      dnames <- c(dnames, paste0("phi-mean(", 1:p, ")[derived]"), paste0("phi-delta(", 1:p, ")[derived]"))
    }
    didx<-c(didx, icur:(icur+p-1))
    icur<-icur+p
  }
  if (bp > 0 && ! fixed_bphi) {
    if (meandelta){
      names <- c(names, paste0("bphi-delta(", 1:bp, ")"))
      sidx<-c(sidx, icur:(icur+bp-1))
      dnames <- c(dnames, paste0("bphi-start(", 1:bp, ")[derived]"), paste0("bphi-end(", 1:bp, ")[derived]"))
    }else{
      names <- c(names, paste0("bphi-end(", 1:bp, ")"))
      dnames <- c(dnames, paste0("bphi-mean(", 1:bp, ")[derived]"), paste0("bphi-delta(", 1:bp, ")[derived]"))
    }
    didx<-c(didx, icur:(icur+bp-1))
    icur<-icur+bp
  }
  if (q > 0 && ! fixed_theta) {
    if (meandelta){
      names <- c(names, paste0("theta-delta(", 1:q, ")"))
      sidx<-c(sidx, icur:(icur+q-1))
      dnames <- c(dnames, paste0("theta-start(", 1:q, ")[derived]"), paste0("theta-end(", 1:q, ")[derived]"))
    }else{
      names <- c(names, paste0("theta-end(", 1:q, ")"))
      dnames <- c(dnames, paste0("theta-mean(", 1:q, ")[derived]"), paste0("theta-delta(", 1:q, ")[derived]"))
    }
    didx<-c(didx, icur:(icur+q-1))
    icur<-icur+q
  }
  if (bq > 0 && ! fixed_btheta) {
    if (meandelta){
      names <- c(names, paste0("btheta-delta(", 1:bq, ")"))
      sidx<-c(sidx, icur:(icur+bq-1))
      dnames <- c(dnames, paste0("btheta-start(", 1:bq, ")[derived]"), paste0("btheta-end(", 1:bq, ")[derived]"))
    }else{
      names <- c(names, paste0("btheta-end(", 1:bq, ")"))
      dnames <- c(dnames, paste0("btheta-mean(", 1:bq, ")[derived]"), paste0("btheta-delta(", 1:bq, ")[derived]"))
    }
    didx<-c(didx, icur:(icur+bq-1))
    icur<-icur+bq
  }
  if (! fixed_var){
    if (meandelta){
       names <- c(names, paste0("var-delta"))
      sidx<-c(sidx,icur)
    }
  }
  return (list(
    names=names,
    dnames=dnames,
    sidx=sidx,
    didx=didx
  ))
}

.pderived<-function(meandelta, p, cov, idx){
  n<-length(idx)
  n2<-n/2
  dp<-array(dim=n2)
  edp<-array(dim=n/2)
  if (meandelta){
    for (i in 1:n2){
      dp[2*i-1]<-p[idx[i]]-p[idx[i+n2]]/2
      dp[2*i]<-p[idx[i]]+p[idx[i+n2]]/2
      if (! is.null(cov)){
        edp[2*i-1]<-cov[idx[i], idx[i]]+cov[idx[i+n2],idx[i+n2]]/4-cov[idx[i],idx[i+n2]]
        edp[2*i]<-cov[idx[i], idx[i]]+cov[idx[i+n2],idx[i+n2]]/4+cov[idx[i],idx[i+n2]]
      }
    }

  }else{
    for (i in 1:n2){
      dp[2*i-1]<-(p[idx[i]]+p[idx[i+n2]])/2
      dp[2*i]<-p[idx[i+n2]]-p[idx[i]]
      if (! is.null(cov)){
        edp[2*i-1]<-(cov[idx[i], idx[i]]+cov[idx[i+n2],idx[i+n2]])/4+cov[idx[i],idx[i+n2]]/2
        edp[2*i]<-(cov[idx[i], idx[i]]+cov[idx[i+n2],idx[i+n2]])-2*cov[idx[i],idx[i+n2]]
      }
    }
  }
  return (list(dp=dp,edp=suppressWarnings(sqrt(edp))))
}





