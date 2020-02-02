df = read.csv("sample_data.csv", header=F)
df[is.na(df)] = 0

sost = function(df, w_size = 8) {
  sos = c()
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

sos = sost(df)
