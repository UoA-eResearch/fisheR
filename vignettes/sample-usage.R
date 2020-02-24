## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- fig.show='hold', fig.cap = "Fisher's Information calculated from the Sample Data"----
library(fisheR)
df = read.csv("../sample_data.csv", header=F)
sos = sost(df)
df_FI = fisher(df, sos = sos, display_plot = TRUE)

