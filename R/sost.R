#' Size of state
#'
#' The sost function calculates the size of state following:
#' Cabezas, H., Campbell, D., Eason, T., Garmestani, A. S., Heberling, M. T., Hopton, M. E., Templeton, J., White, D., Zanowick, M., & Sparks, R. T. (2010). San Luis Basin sustainability metrics project: A methodology for evaluating regional sustainability. USEPA. USA.
#' Eason, T., & Cabezas, H. (2012). Evaluating the sustainability of a regional system using Fisher information in the San Luis Basin, Colorado. Journal of Environmental Management, 94(1), 41–49. https://doi.org/10.1016/j.jenvman.2011.08.003
#' Ahmad, N., Derrible, S., Eason, T., & Cabezas, H. (2016). Using Fisher information to track stability in multivariate systems. Royal Society Open Science, 3(11), 160582. https://doi.org/10.1098/rsos.160582
#'
#' @param df A dataframe containing time information and variables of interest
#' @param w_size Window size
#'
#' @return A vector containing size of state values with length equal to the number of variable columns in your input df


sost = function(df, w_size = 8) {

  if (any(unlist(lapply(df, is.numeric)) == FALSE)) {
    stop("All columns must be numeric, check structure. First column should be time and all following columns the variables of interest")
  } else {
    cat("Structure seems good, figuring out size of state")
  }

  if (is_tibble(df)) {
    warning("Tibble format not suited to fishR, converting to data.frame")
    df <- as.data.frame(df)
  }

  sos = c()
  df[is.na(df)] = 0
  for (j in 2:ncol(df)) {
    sos_temp = c()
    for (i in 1:nrow(df)) {
      A=na.omit(df[i:(i+w_size - 1),j])
      if (length(A) == w_size && !any(A == 0)) {
        sos_temp = c(sos_temp, sd(A))
      }
    }
    if (length(sos_temp) == 0) {
      sos = c(sos, 0)
    } else {
      sos = c(sos, min(sos_temp) * 2)
    }
  }
  sos
}
