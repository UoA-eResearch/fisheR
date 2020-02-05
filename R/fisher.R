#' Fisher function
#'
#' The fisher function calculates the Fisher's Information following:
#' Cabezas, H., Campbell, D., Eason, T., Garmestani, A. S., Heberling, M. T., Hopton, M. E., Templeton, J., White, D., Zanowick, M., & Sparks, R. T. (2010). San Luis Basin sustainability metrics project: A methodology for evaluating regional sustainability. USEPA. USA.
#' Eason, T., & Cabezas, H. (2012). Evaluating the sustainability of a regional system using Fisher information in the San Luis Basin, Colorado. Journal of Environmental Management, 94(1), 41–49. https://doi.org/10.1016/j.jenvman.2011.08.003
#' Ahmad, N., Derrible, S., Eason, T., & Cabezas, H. (2016). Using Fisher information to track stability in multivariate systems. Royal Society Open Science, 3(11), 160582. https://doi.org/10.1098/rsos.160582
#'
#' @param df A dataframe containing time information and variables of interest
#' @param sos An option vector containing size of state
#' @param w_size Window size
#' @param w_incr Window increment
#' @param smooth_step Window size for the smoothing operation
#' @param xtick_step Number of ticks on the x axis (currently unused)
#' @returns A dataframe where the last three columns are the Fisher's Information means, Fisher's Information smoothed and time-steps

fisher = function(df, sos = c(), w_size = 8, w_incre = 1, smooth_step = 3, xtick_step = 1) {
  start_time <- as.numeric(Sys.time())
  if (length(sos) == 0) {
    sos = sost(df)
  }
  df[is.na(df)] = 0
  FI_final = c()
  k_init = c()
  for (i in seq(1, nrow(df), w_incre)) {
    Data_win = na.omit(df[i:(i+w_size - 1),2:ncol(df)])
    if (nrow(Data_win) == w_size) {
      Bin = c()
      for (m in 1:w_size) {
        Bin_temp = c()
        for (n in 1:w_size) {
          if (m == n) {
            Bin_temp = c(Bin_temp, "I")
          } else {
            Bin_temp_1 = 0
            for (k in 1:ncol(Data_win)) {
              if (abs(Data_win[m,k] - Data_win[n,k]) <= sos[k]) {
                Bin_temp_1 = Bin_temp_1 + 1
              }
            }
            Bin_temp = c(Bin_temp, Bin_temp_1)
          }
        }
        Bin = rbind(Bin, Bin_temp)
      }
      #print(Bin)
      FI = c()
      for (tl in 1:100) {
        tl1 = length(sos) * tl / 100
        Bin_1 = c()
        Bin_2 = c()
        for (j in 1:w_size) {
          if (!(j %in% Bin_2)) {
            Bin_1_temp = c(j)
            for (i in 1:ncol(Bin)) {
              if (Bin[j,i] != "I" && as.numeric(Bin[j,i]) >= tl1 && !(i %in% Bin_2)) {
                Bin_1_temp = c(Bin_1_temp, i)
              }
            }
            Bin_1 = c(Bin_1, list(Bin_1_temp))
            Bin_2 = c(Bin_2, Bin_1_temp)
          }
        }
        prob = c(0, lengths(Bin_1) / length(Bin_2), 0)
        prob_q = sqrt(prob)
        FI_temp = 0
        for (i in 1:(length(prob_q) - 1)) {
          FI_temp = FI_temp + (prob_q[i] - prob_q[i + 1]) ** 2
        }
        FI_temp = 4 * FI_temp
        FI = c(FI, FI_temp)

      }

      k_init = c(k_init, which(FI != 8)[1])

      FI_final = rbind(FI_final, FI)
    }
  }
  if (length(k_init) == 0) {
    k_init = c(1)
  }
  FI_means = rowMeans(FI_final[,min(k_init):ncol(FI_final)])
  time_windows = df[1:nrow(FI_final) * w_incre + w_size - 1, 1]

  FI_smth = c()
  for (i in seq(smooth_step + 1, length(FI_means)+smooth_step, smooth_step)) {
    for (j in 1:smooth_step) {
      FI_smth = c(FI_smth, mean(FI_means[(i-smooth_step):(i - 1)], na.rm=TRUE))
    }
  }
  FI_smth = FI_smth[1:length(FI_means)]
  FI_final = cbind(FI_final, FI_means, FI_smth, time_windows)

  rownames(FI_final) = NULL
  df_FI = as.data.frame(FI_final)

  #if (write_out_csv == TRUE) {}
  write.table(df_FI, "FI.csv", sep=",", col.names = FALSE, row.names = FALSE)

  #if (write_out_rds == TRUE)
  #saveRDS(df_FI, "FI.RData")

  #if (plot_out == TURE) {}
  plot(df_FI$time_windows, df_FI$FI_means, type="l", col="blue", xlab = "Time Step", ylab = "Fisher Information")
  lines(df_FI$time_windows, df_FI$FI_smth, type="l", col="red")
  print(sprintf("Completed in %.2f seconds", as.numeric(Sys.time()) - start_time))
  df_FI

}
