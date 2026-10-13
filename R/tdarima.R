#' @include utils.R
NULL

LTDARIMA<-'JD3_LTDARIMA'
LTDARIMA_LL<-'JD3_LTDARIMA_LIKELIHOOD'

LTDARIMA_FIXED<-'JD3_LTDARIMA_FIXED'
LTDARIMA_TD<-'JD3_LTDARIMA_TD'

TDARIMA_DECOMPOSITION<-'JD3_TDARIMA_DECOMPOSITION_RSLTS'

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
                              fixed_theta=FALSE, fixed_btheta=FALSE, fixed_var=TRUE, parametrization="mean_delta", eps=1e-7, simplificationThreshold=0){
  if (! is.ts(data)) stop("data should be a time series (ts)")
  jrslt<-.jcall('jdplus/advancedsa/base/r/TimeVaryingArimaModels', 'Ljdplus/advancedsa/base/core/tdarima/LtdResults;', 'estimate',
                as.numeric(data),
                as.integer(frequency(data)),
                as.logical(mean),
                rjd3toolkit::.r2jd_matrix(X),
                as.integer(regular),
                as.integer(seasonal),
                as.logical(fixed_phi),
                as.logical(fixed_bphi),
                as.logical(fixed_theta),
                as.logical(fixed_btheta),
                as.logical(fixed_var),
                as.character(parametrization),
                as.numeric(eps),
                as.numeric(simplificationThreshold))
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
#' ts.plot(ts.union(s, sa$components[,1], sa$components[,1]+sa$components_stdev[,1]), col=c('gray', 'blue', 'magenta'))
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
    # the stderr are linear, not the var !!
    et<-.linear(1, sqrt(var1), n)
    var<-matrix(et*et, nrow = 1, ncol=n)
    p<-rbind(p, var)
  }
  if (is.null(p)){
    jp<-.jnull("jdplus/toolkit/base/api/math/matrices/Matrix")
  }
  else{
    jp<-.r2jd_matrix(p)
  }

  jrslts<-.jcall('jdplus/advancedsa/base/r/TimeVaryingArimaModels', 'Ljdplus/advancedsa/base/r/TimeVaryingArimaModels$DecompositionResults;', 'arimaDecomposition',
                  as.numeric(data), as.integer(frequency(data)), as.integer(regular), as.integer(seasonal), as.logical(var1 != 1), jp, as.logical(se))
  if (is.jnull(jrslts)) return (NULL)
  jcmps<-.jcall(jrslts, 'Ljdplus/toolkit/base/api/math/matrices/Matrix;', 'getComponents')
  cmps<-rjd3toolkit::.jd2r_matrix(jcmps)
  ncmps<-dim(cmps)[2]
  if (ncmps == 3)  names<-c("trend", "seasonal", "irregular")else names <-c("trend", "irregular")
  colnames(cmps)<-names
  cmps<-ts(cmps, frequency=frequency(data), start=start(data))
  if (se){
    jecmps<-.jcall(jrslts, 'Ljdplus/toolkit/base/api/math/matrices/Matrix;', 'getComponentsStdev')
    ecmps<-rjd3toolkit::.jd2r_matrix(jecmps)
    enames<-paste0(names, "-stdev")
    colnames(ecmps)<-enames
    ecmps<-ts(ecmps, frequency=frequency(data), start=start(data))
  }

  return (structure(list(
    components=cmps,
    components_stdev=ecmps,
    direct_likelihood=structure(.proc_likelihood(jrslts, "direct_ll."), class=LTDARIMA_LL),
    decomposition_likelihood=structure(.proc_likelihood(jrslts, "decomposition_ll."), class=LTDARIMA_LL),
    direct_res=.proc_vector(jrslts, "direct_res"),
    decomposition_res=.proc_vector(jrslts, "decomposition_res")
  ), class=TDARIMA_DECOMPOSITION))
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
#' ts.plot(ts.union(s, sa$components[,1], sa$components[,1]+sa$components[,3]), col=c('gray', 'blue', 'magenta'))
#' ts.plot(sa$components_stdev, col=c('red', 'blue', 'magenta'))
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
  jrslts<-.jcall('jdplus/advancedsa/base/r/TimeVaryingArimaModels', 'Ljdplus/advancedsa/base/r/TimeVaryingArimaModels$DecompositionResults;', 'arimaDecomposition',
                 as.numeric(data), as.integer(frequency(data)), as.integer(regular), as.integer(seasonal), as.logical(var), jp, as.logical(se))
  if (is.jnull(jrslts)) return (NULL)
  jcmps<-.jcall(jrslts, 'Ljdplus/toolkit/base/api/math/matrices/Matrix;', 'getComponents')
  cmps<-rjd3toolkit::.jd2r_matrix(jcmps)
  ncmps<-dim(cmps)[2]
  if (ncmps == 3)  names<-c("trend", "seasonal", "irregular")else names <-c("trend", "irregular")
  colnames(cmps)<-names
  cmps<-ts(cmps, frequency=frequency(data), start=start(data))
  if (se){
    jecmps<-.jcall(jrslts, 'Ljdplus/toolkit/base/api/math/matrices/Matrix;', 'getComponentsStdev')
    ecmps<-rjd3toolkit::.jd2r_matrix(jecmps)
    enames<-paste0(names, "-stdev")
    colnames(ecmps)<-enames
    ecmps<-ts(ecmps, frequency=frequency(data), start=start(data))
  }

  return (structure(list(
    components=cmps,
    components_stdev=ecmps,
    direct_likelihood=structure(.proc_likelihood(jrslts, "direct_ll."), class=LTDARIMA_LL),
    decomposition_likelihood=structure(.proc_likelihood(jrslts, "decomposition_ll."), class=LTDARIMA_LL),
    direct_res=.proc_vector(jrslts, "direct_res"),
    decomposition_res=.proc_vector(jrslts, "decomposition_res")
  ), class=TDARIMA_DECOMPOSITION))
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
#' @param parametrization Type of the parametrization of the linear time-dependent parameters:
#' first and last parameters or average and variation of the parameters
#' @param eps The precision of the optimization procedure
#' @param simplification_threshold p-value for non significant parameters (simplified)
#' @param decomposition Indicates that a decomposition of the series must be processed
#' @param regeffects In case of regression variables (excluding the mean correction),
#' the user must specify to which component the regression effects are associated
#' (1 for the trend, 2 for the seasonal and 3 for the irregular). If not specified, the
#' final components are not computed.
#' @param clean Cleans missing values at the beginning/end of the series. Should be
#' set to false if you want forecasts/backcasts
#'
#' @returns A JD3_LTDARIMA object
#' @export
#'
#' @examples
#' s<-rjd3toolkit::ABS$X0.2.09.10.M
#' td<-rjd3toolkit::td(s=s)
#' q<-rjd3sax::ltdarima_estimation(s, regular=c(0,1,1), seasonal=c(0,1,1), X=td,
#'  fixed_var = FALSE, eps=1e-15, parametrization = "mean_delta")


ltdarima_estimation<-function(data, mean=FALSE, X=NULL, regular=c(0,1,1), seasonal=c(0,1,1), fixed_phi=TRUE, fixed_bphi=TRUE,
                                fixed_theta=FALSE, fixed_btheta=FALSE, fixed_var=TRUE, parametrization=c("mean_delta", "start_end"),
                              eps=1e-7, simplification_threshold=0,
                              decomposition=TRUE, regeffects=NULL, clean=TRUE){

  parametrization=match.arg(parametrization)
  if (clean)
    data<-clean_extremities(data)
  jrslt<-.jestimation(data, mean, X, regular, seasonal, fixed_phi, fixed_bphi, fixed_theta, fixed_btheta, fixed_var, parametrization, eps,simplification_threshold)

  initial<-.estimation_initial(jrslt, data, decomposition, X, regeffects, regular, seasonal)
  final<-NULL
  niter<-.proc_int(jrslt, "niters")
  tv<-.proc_bool(jrslt, "time_varying")
  if (tv)
    niter<-niter-1
  steps<-NULL
  final<-NULL
  if (niter>0){
    steps<-list()
    for (iter in 1:niter){
      name<-paste0("iter-", as.character(iter))
      id<-paste0("iter(", as.character(iter-1), ")")
      steps[[name]]<-.estimation_step(jrslt, id, data, FALSE, X, regeffects, regular, seasonal, parametrization, fixed_var)
    }
  }
  if(tv){
    final<-.estimation_step(jrslt, "final", data, decomposition, X, regeffects, regular, seasonal, parametrization, fixed_var)
  }

  return(structure(list(initial=initial, steps=steps, final=final), class=LTDARIMA))
}

.item<-function(...){
  return (paste(..., sep="."))
}


S0<-"sarima"
S1<-"final"
REG<-"regression"
LIN<-"y_lin"
M<-"model"
E<-"residuals"
LL<-"likelihood."


.estimation_initial<-function(jrslt, data, decomposition, X, regeffects, regular, seasonal){
  freq<-frequency(data)
  start=start(data)

  regs<-.proc_vector(jrslt, .item(S0, REG, "effect"))
  tsregs <-NULL;
  if ( ! is.null(regs)){
    tsregs<-ts(data=regs, frequency = freq, start = start)
  }
  parameters<-.proc_vector(jrslt, .item(S0, M, "p"))
  flin<-ts(data=.proc_vector(jrslt, .item(S0, REG, LIN)), frequency = freq, start = start)
  # initial (fixed) model
  initial=list(
    model=list(
      period=freq,
      regular=regular,
      seasonal=seasonal,
      parameters=parameters,
      covariance=.proc_matrix(jrslt, .item(S0, M, "p_cov"))),
    likelihood=structure(.proc_likelihood(jrslt, .item(S0, LL)), class=LTDARIMA_LL),
    regression=list(
      coefficients=.proc_vector(jrslt, .item(S0, REG, "c")),
      covariance=.proc_matrix(jrslt, .item(S0, REG, "cov")),
      linearized=ts(data=flin),
      regression_effect=tsregs),
    residuals=list(
      res=.proc_vector(jrslt, .item(S0, E, "res")),
      type=.proc_data(jrslt, .item(S0, E, "type")),
      df=.proc_int(jrslt, .item(S0, E, "df")),
      dfc=.proc_int(jrslt, .item(S0, E, "dfc")),
      mean=.proc_test(jrslt, .item(S0, E, "mean")),
      skewness=.proc_test(jrslt, .item(S0, E, "skewness")),
      kurtosis=.proc_test(jrslt, .item(S0, E, "kurtosis")),
      normality=.proc_test(jrslt, .item(S0, E, "doornikhansen")),
      ljung_box=.proc_test(jrslt, .item(S0, E, "lb")),
      seasonal_ljung_box=.proc_test(jrslt, .item(S0, E, "seaslb")),
      nruns=.proc_test(jrslt, .item(S0, E, "nruns")),
      lruns=.proc_test(jrslt, .item(S0, E, "lruns")),
      nudruns=.proc_test(jrslt, .item(S0, E, "nudruns")),
      ludruns=.proc_test(jrslt, .item(S0, E, "ludruns"))
    )
  )
  if (decomposition){
    sa0<-NULL
    s<-data
    cmps<-.ucm(flin, regular, seasonal, parameters)
    if (! is.null(cmps)){
      components<-list(
        y=flin,
        sa=cmps[,1]+cmps[,3],
        trend=cmps[,1],
        seas=cmps[,2],
        irregular=cmps[,3],
        trend_stdev=cmps[,4],
        seas_stdev=cmps[,5],
        irregular_stdev=cmps[,6]
      )
      finals<-NULL
      if (!is.null(X) && ! is.null(regeffects)){

      }
      sa0<-list(
        components=components,
        finals=finals
      )
      initial[["decomposition"]]=sa0
    }
  }
  return (structure(initial, class=LTDARIMA_FIXED))
}

.estimation_step<-function(jrslt, step, data, decomposition, X, regeffects, regular, seasonal, parametrization, fixed_var){
  freq<-frequency(data)
  start=start(data)

  regs<-.proc_vector(jrslt, .item(step, REG, "effect"))
  tsregs <-NULL;
  if ( ! is.null(regs)){
    tsregs<-ts(data=regs, frequency = freq, start = start)
  }

  flin<-ts(data=.proc_vector(jrslt, .item(step, REG, "y_lin")), frequency = freq, start = start)
  parameters<-.proc_vector(jrslt, .item(step, M, "pall"))
  pcovariance<-.proc_matrix(jrslt, .item(step, M, "pall_cov"))
  dpcovariance<-.proc_matrix(jrslt, .item(step, M, "pderived_cov"))
  m<-length(data)
  # step model
  rslt<-list(
    model=list(
      type=parametrization,
      n=m,
      period=freq,
      regular=regular,
      seasonal=seasonal,
      parameters_names=.proc_vector(jrslt, .item(step, M, "pnames")),
      parameters=parameters,
      parameters_covariance=pcovariance,
      parameters_stde=suppressWarnings(sqrt(diag(pcovariance))),
      derived_parameters_names=.proc_vector(jrslt, .item(step, M, "pderived_names")),
      derived_parameters=.proc_vector(jrslt, .item(step, M, "pderived")),
      derived_parameters_covariance=dpcovariance,
      derived_parameters_stde=suppressWarnings(sqrt(diag(dpcovariance))),
      p0=.proc_vector(jrslt, .item(step, M, "p0")),
      p1=.proc_vector(jrslt, .item(step, M, "p1")),
      pmean=.proc_vector(jrslt, .item(step, M, "pmean")),
      pdelta=.proc_vector(jrslt, .item(step, M, "pdelta")),
      covariance=pcovariance,
      lr_test=.proc_test(jrslt, .item(step, M, "lrtest")),
      stationarity_test=.proc_test(jrslt, .item(step, M, "test_stationarity"))
    ),
    # scores=.proc_vector(jrslt, "regression.ml.score1"),
    # information=.proc_matrix(jrslt, "regression.ml.information1")),
    likelihood=structure(rjd3toolkit::.proc_likelihood(jrslt, .item(step, LL)), class=LTDARIMA_LL),
    regression=list(
      coefficients=.proc_vector(jrslt, .item(step, REG, "c")),
      covariance=.proc_matrix(jrslt, .item(step, REG, "cov")),
      linearized=flin,
      regression_effect=tsregs
    ),
    residuals=list(
      res=.proc_vector(jrslt, .item(step, E, "res")),
      type=.proc_data(jrslt, .item(step, E, "type")),
      df=.proc_int(jrslt, .item(step, E, "df")),
      dfc=.proc_int(jrslt, .item(step, E, "dfc")),
      mean=.proc_test(jrslt, .item(step, E, "mean")),
      skewness=.proc_test(jrslt, .item(step, E, "skewness")),
      kurtosis=.proc_test(jrslt, .item(step, E, "kurtosis")),
      normality=.proc_test(jrslt, .item(step, E, "doornikhansen")),
      ljung_box=.proc_test(jrslt, .item(step, E, "lb")),
      seasonal_ljung_box=.proc_test(jrslt, .item(step, E, "seaslb")),
      nruns=.proc_test(jrslt, .item(step, E, "nruns")),
      lruns=.proc_test(jrslt, .item(step, E, "lruns")),
      nudruns=.proc_test(jrslt, .item(step, E, "nudruns")),
      ludruns=.proc_test(jrslt, .item(step, E, "ludruns"))
    )
  )

  if (decomposition == TRUE){
    sa1<-NULL
    s<-data
    var<-1
    if (! fixed_var){
      fp<-rslt$model$parameters
      e<-fp[length(fp)]
      if ( parametrization == "mean_delta"){
        e<-1+e*(m-1)
        var<-e*e
      }
      var<-e*e
    }
    cmps<-ltdarima_decomposition(flin, regular = regular, seasonal = seasonal,
                                 p0=rslt$model$p0, p1=rslt$model$p1, var1 = var, se=TRUE)
    if (! is.null(cmps)){
      components<-list(
        y=flin,
        sa=cmps$components[,1]+cmps$components[,3],
        trend=cmps$components[,1],
        seas=cmps$components[,2],
        irregular=cmps$components[,3],
        trend_stdev=cmps$components_stdev[,1],
        seas_stdev=cmps$components_stdev[,2],
        irregular_stdev=cmps$components_stdev[,3]
      )
      finals<-NULL
      if (!is.null(X) && ! is.null(regeffects)){

      }
      sa1<-list(
        components=components,
        finals=finals,
        direct_likelihood=cmps$direct_likelihood,
        decomposition_likelihood=cmps$decomposition_likelihood,
        direct_res=cmps$direct_res,
        decomposition_res=cmps$decomposition_res
      )
    }
    rslt[["decomposition"]]=sa1
  }
  return (structure(rslt, class=LTDARIMA_TD))
}


.sarima<-function(period, regular, seasonal, w){
  phi<-NULL
  bphi<-NULL
  theta<-NULL
  btheta<-NULL

  p = regular[1]
  d = regular[2]
  q = regular[3]
  bp = seasonal[1]
  bd = seasonal[2]
  bq = seasonal[3]
  i<-1
  if (p>0){
    phi<-w[i:(i+p-1)]
    i<-i+p
  }
  if (bp>0){
    bphi<-w[i:(i+bp-1)]
    i<-i+bp
  }
  if (q>0){
    theta<-w[i:(i+q-1)]
    i<-i+q
  }
  if (bq>0){
    btheta<-w[i:(i+bq-1)]
  }

  return (rjd3toolkit::sarima_model("", period, phi, d, theta, bphi, bd, btheta))
}

.ucm<-function(data, regular, seasonal, p){
  freq<-frequency(data)
  start<-start(data)
  sarima<-.sarima(freq, regular, seasonal, p)
  ucm<-rjd3toolkit::sarima_decompose(sarima)
  if (is.null(ucm)) return (NULL)
  return (ts(rjd3toolkit::ucarima_estimate(data, ucm), frequency = freq, start = start) )
}

