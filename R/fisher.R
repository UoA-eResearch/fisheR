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
#library(matrixStats)

IVAL <- -Inf

fisher = function(df, sos = c(), w_size = 8, w_incre = 1, smooth_step = 3, RedRum = FALSE, write_out_csv = FALSE, write_out_rds = FALSE, display_plot = FALSE) {
  start_time <- as.numeric(Sys.time())

  if (class(df)[1] != "data.frame") {
    warning("fisher requires a dataframe (tibble format not suitable), attempting to converting to data.frame")
    df <- as.data.frame(df)
  }

  if (any(unlist(lapply(df, is.numeric)) == FALSE)) {
    stop("All columns must be numeric, check structure. First column should be time and all following columns the variables of interest")
  } else {
    cat("Structure seems good, Cookie Monster says \"num, num, num\", let's go fishing... \n")
  }

  if (length(sos) == 0) {
    sos = sost(df, w_size = w_size)
  }
  df[is.na(df)] = 0
  FI_final = c()
  k_init = c()
  window_seq = seq(1, nrow(df), w_incre)
  number_of_states_per_tl = matrix(ncol = 100, nrow = length(window_seq) - (w_size / w_incre - 1))
  rownames(number_of_states_per_tl) = paste0("wi", window_seq[1:nrow(number_of_states_per_tl)])
  colnames(number_of_states_per_tl) = paste0("tl", 1:100)
  for (i in window_seq) {
    window_index_str = paste0("wi", as.character(i))
    Data_win = as.matrix(na.omit(df[i:(i+w_size - 1),2:ncol(df)]), byrow = TRUE)
    if (nrow(Data_win) == w_size) {
      Bin = c()
      for (m in 1:w_size) {
        Bin_temp = c()
        for (n in 1:w_size) {
          if (m == n) {
            Bin_temp = c(Bin_temp, IVAL)
          } else {
            Bin_temp_1 <- sum( abs(Data_win[m,] - Data_win[n,]) <= sos )
            Bin_temp = c(Bin_temp, Bin_temp_1)
          }
        }
        Bin = rbind(Bin, Bin_temp)
      }
      #print(Bin)
      FI = c()
      i_vals <- 1:ncol(Bin)
      for (tl in 1:100) {
        if (RedRum){
          print("All work and no play makes Jack a dull boy")
        }
        tl1 = length(sos) * tl / 100
        Bin_1 = c()
        Bin_2 = c()
        for (j in 1:w_size) {
          if (!(j %in% Bin_2)) {
            binj <- Bin[j,]
            # note: the last test is redundant since IVAL = -Inf
            Bin_1_temp <- c(j, which(!(i_vals %in% Bin_2) & (binj >= tl1) & (binj != IVAL)))
            Bin_1 = c(Bin_1, list(Bin_1_temp))
            Bin_2 = c(Bin_2, Bin_1_temp)
          }
        }

        number_of_states_per_tl[window_index_str, tl] = length(Bin_1)

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
  number_of_states_per_tl <- number_of_states_per_tl[rowSums(is.na(number_of_states_per_tl)) != ncol(number_of_states_per_tl), ]
  FI_means = rowMeans(FI_final[,min(k_init):ncol(FI_final)])
  time_windows = df[1:nrow(FI_final) * w_incre + w_size - 1, 1]
  mean_no_states = rowMeans(number_of_states_per_tl[,min(k_init):ncol(number_of_states_per_tl)]) # function needs to be checked
  median_no_states = matrixStats::rowMedians(number_of_states_per_tl[,min(k_init):ncol(number_of_states_per_tl)]) # function needs to be checked
  sum_of_states = rowSums(number_of_states_per_tl[,min(k_init):ncol(number_of_states_per_tl)]) # function needs to be checked


  FI_smth = c()
  for (i in seq(smooth_step + 1, length(FI_means)+smooth_step, smooth_step)) {
    for (j in 1:smooth_step) {
      FI_smth = c(FI_smth, mean(FI_means[(i-smooth_step):(i - 1)], na.rm=TRUE))
    }
  }
  FI_smth = FI_smth[1:length(FI_means)]
  FI_final = cbind(FI_final, FI_means, FI_smth, sum_of_states, mean_no_states, median_no_states, time_windows)

  rownames(FI_final) = NULL
  df_FI = as.data.frame(FI_final)
  df_FI = cbind(number_of_states_per_tl, df_FI)


  if (write_out_csv == TRUE) {
  write.table(df_FI, "FI.csv", sep=",", col.names = FALSE, row.names = FALSE)
  }

  if (write_out_rds == TRUE) {
  saveRDS(df_FI, "FI.RData")
  }

  if (display_plot == TRUE) {
  plot(df_FI$time_windows, df_FI$FI_means, type="l", col="blue", xlab = "Time Step", ylab = "Fisher Information")
  lines(df_FI$time_windows, df_FI$FI_smth, type="l", col="red")
  }
  print(sprintf("Completed in %.2f seconds", as.numeric(Sys.time()) - start_time))
  print(sprintf("Checksum: %.10f", sum(df_FI)))
  #list(df = df_FI, number_of_states_per_tl = number_of_states_per_tl)
  df_FI
}
