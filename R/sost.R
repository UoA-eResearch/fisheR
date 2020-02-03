# The sost function calculates the size of state following:
# Cabezas, H., Campbell, D., Eason, T., Garmestani, A. S., Heberling, M. T., Hopton, M. E., Templeton, J., White, D., Zanowick, M., & Sparks, R. T. (2010). San Luis Basin sustainability metrics project: A methodology for evaluating regional sustainability. USEPA. USA.
# Eason, T., & Cabezas, H. (2012). Evaluating the sustainability of a regional system using Fisher information in the San Luis Basin, Colorado. Journal of Environmental Management, 94(1), 41–49. https://doi.org/10.1016/j.jenvman.2011.08.003
# Ahmad, N., Derrible, S., Eason, T., & Cabezas, H. (2016). Using Fisher information to track stability in multivariate systems. Royal Society Open Science, 3(11), 160582. https://doi.org/10.1098/rsos.160582

sost = function(df, w_size = 8) {
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
