#' @include tdgeneric.R


#' @export
plot.JD3_LTDARIMA_RSLTS <- function(x, first_date = NULL, last_date = NULL,
                                      type_chart = c("sa-trend", "seas-irr"),
                                      caption = c(
                                        "sa-trend" = "Y, Sa, trend",
                                        "seas-irr" = "Sea., irr."
                                      )[type_chart],
                                      colors = c(
                                        y = "#F0B400", t = "#1E6C0B", sa = "#155692",
                                        s = "#1E6C0B", i = "#155692"
                                      ),
                                      ...) {
  if (is.null(x$final$decomposition)){
    return(NULL)
  }
  plot(sa_decomposition(x),
       first_date = first_date, last_date = last_date,
       type_chart = type_chart,
       caption = caption,
       colors = colors,
       ...
  )
}
