fisher = function(df, sos = c(), w_size = 8, w_incre = 1) {
  if (length(sos) == 0) {
    sos = sost(df)
  }
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
              if (Bin[j,i] != "I" && Bin[j,i] >= tl1 && !(i %in% Bin_2)) {
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
  FI_final = cbind(FI_final, FI_means, time_windows)
  df_FI = as.data.frame(FI_final)
  write.table(df_FI, "FI.csv", sep=",", col.names = FALSE, row.names = FALSE)
  df_FI

}

df_FI = fisher(df, sos = sos)
