#' @importFrom stats printCoefmat end time
#' @importFrom utils capture.output

#' @export
print.JD3_LTDARIMA_RSLTS <- function(x, digits = max(3L, getOption("digits") - 3L), summary_info = getOption("summary_info"),
                                       ...) {
  cat("Model: Time dependent SARIMA", "\n", sep = "")
  if (summary_info) {
    cat("\nFor a more detailed output, use the 'summary()' function.\n")
  }
  printInitialModel(x$initial, digits=digits, ...)
  return(invisible(x))
}

#' @export
summary.JD3_LTDARIMA_RSLTS  <- function(x, digits = max(3L, getOption("digits") - 3L), signif.stars = getOption("show.signif.stars"), ...) {
  cat("Model: Time dependent SARIMA", "\n", sep = "")
  summaryInitialModel(x$initial, ...)
}


printInitialModel <- function(x, digits = digits, ...) {
  sarima<-.sarima_coef_table(x, ...)
  cat("\n", "Initial SARIMA model", "\n", sep = "")
  cat(.arima_node(sarima$sarima_orders$p, sarima$sarima_orders$d, sarima$sarima_orders$q),
      .arima_node(sarima$sarima_orders$bp, sarima$sarima_orders$bd, sarima$sarima_orders$bq),"\n")
  if (!is.null(sarima$coef_table)) {
    print(sarima$coef_table, digits = digits, na.print = "NA", ...)
  }
  invisible(x)
}

summaryInitialModel <- function(x, digits = max(3L, getOption("digits") - 3L), signif.stars = getOption("show.signif.stars"), ...) {
  cat("\n", "Initial SARIMA model", "\n", sep = "")
  cat("\nSARIMA coefficients:\n")
}

#' @rdname jd3_print
#' @export
print.JD3_LTDARIMA_LIKELIHOOD <- function(x, ...) {
  ll <- x
  cat("Number of observations:", ll$nobs, "\n")
  cat("Number of effective observations:", ll$neffective, "\n")
  cat("Number of parameters:", ll$nparams, "\n\n")
  cat("Loglikelihood:", ll$ll, "\n")
  cat("Standard error of the regression (ML estimate):", sqrt(ll$ssq / ll$neffective), "\n")
  cat("AIC:", ll$aic, "\n")
  cat("AICC:", ll$aicc, "\n")
  cat("BIC:", ll$bic, "\n\n")
  invisible(x)
}
#' @export
summary.JD3_LTDARIMA_LIKELIHOOD <- function(object, ...) {
  print(object)
}

.sarima_coef_table <- function(x, ...) {
  model<-x$model
  ll<-x$likelihood
  ndf <- ll$neffective - ll$nparams
  p = model$regular[1]
  d = model$regular[2]
  q = model$regular[3]
  bp = model$seasonal[1]
  bd = model$seasonal[2]
  bq = model$seasonal[3]
  period = model$period
  sarima_orders <- list(
    p=p,d=d,q=q,bp=bp,bd=bd,bq=bq,period=period
  )
  estimate <- model$parameters
  names <- NULL
  if (p > 0) {
    names <- c(names, paste0("phi(", 1:p, ")"))
  }
  if (bp > 0) {
    names <- c(names, paste0("bphi(", 1:bp, ")"))
  }
  if (q > 0) {
    names <- c(names, paste0("theta(", 1:q, ")"))
  }
  if (bq > 0) {
    names <- c(names, paste0("btheta(", 1:bq, ")"))
  }

  if (length(estimate) > 0) {
    stde <- sqrt(diag(model$covariance))
    t <- estimate / stde
    pval <- 2 * pt(abs(t), ndf, lower.tail = FALSE)
    table <- data.frame(estimate, stde, t, pval,
                        stringsAsFactors = FALSE)
    colnames(table) <- c(
      "Estimate", "Std. Error",
      "T-stat", "Pr(>|t|)"
    )
    rownames(table) <- names
  } else {
    table <- NULL
  }
  return (list(
    sarima_orders = sarima_orders,
    coef_table = table)
  )
}

.arima_node <- function(p, d, q) {
  s <- paste(p, d, q, sep = ",")
  return(paste0("(", s, ")"))
}


