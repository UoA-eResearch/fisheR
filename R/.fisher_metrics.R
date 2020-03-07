#' Metric functions
#'
#' These metric functions, currently in development, are made to collapse the Fisher's Information output into single number metrics
#'
#' @param fisher_df A fisher output containing FI_means and time_windows
#' @param winds Window size
#' @param std Number of standard deviations to look beyond...
#'
#' @return a single row tibble with 72 metrics


library(moments)
library(psych)
library(tidyverse)
library(changepoint)
library(changepoint.np)
library(matrixStats)

# sample_data <- read.csv("C:/Users/qase352/Dropbox/QuinnAsenaPhD/R/fishers_information/fi_ahmad_2016/FI_Scripts_v2.00/sample_data_FI.csv")
# head(sample_data)
# sample_data <- sample_data %>%
#   select(c(-X))
# sample_fish <- fisher(sample_data, display_plot = T)
#
# qdata <- readRDS("C:/Users/qase352/Desktop/1_complete_output.RData")
# sub_qdata <- as.data.frame(qdata[[1]][[4]][1:500, 1:50])
# time <- 1:nrow(sub_qdata)
# sub_qdata <- cbind(time,sub_qdata)
# head(sub_qdata)
# qfish <- fisher(sub_qdata, display_plot = T)
# plot(qfish$time_windows, qfish$FI_means, type="l", col="blue", xlab = "Time Step", ylab = "Fisher Information")
# lines(qfish$time_windows, qfish$FI_smth, type="l", col="red")


##### basic summary statistics
summary_stats_func <- function(fisher_df) {
    tib <- tibble(
    max_fish = max(fisher_df$FI_means),
    min_fish = min(fisher_df$FI_means),
    range_fish = max_fish - min_fish,
    mid_range = (max_fish + min_fish) / 2,
    mu = mean(fisher_df$FI_means),
    std = sd(fisher_df$FI_means),
    std_mean_ratio <- std/mu,
    variance = var(fisher_df$FI_means),
    skew = skewness(fisher_df$FI_means),
    kurt = kurtosis(fisher_df$FI_means),
    mad = median(abs(fisher_df$FI_means - median(fisher_df$FI_means))), # Median Absolute Deviation

    max_percent_diff_from_median = max(fisher_df$FI_means - median(fisher_df$FI_means)) / median(fisher_df$FI_means),
    min_percent_diff_from_median = min(fisher_df$FI_means - median(fisher_df$FI_means)) / median(fisher_df$FI_means),
    abs_percent_diff_from_median = max(abs(fisher_df$FI_means - median(fisher_df$FI_means))) / median(fisher_df$FI_means),

    max_percent_diff_from_mean = max(fisher_df$FI_means - mean(fisher_df$FI_means)) / mean(fisher_df$FI_means),
    min_percent_diff_from_mean = min(fisher_df$FI_means - mean(fisher_df$FI_means)) / mean(fisher_df$FI_means),
    abs_percent_diff_from_mean = max(abs(fisher_df$FI_means - mean(fisher_df$FI_means))) / mean(fisher_df$FI_means),
    )
    tib
}
#x <- summary_stats_func(fisher_df)


##### percent rate of change and slope
slope_rates_func <- function(fisher_df) {
  percent_roc = 100 * diff(fisher_df$FI_means) / fisher_df[-nrow(fisher_df), ]$FI_means
  slope = diff(fisher_df$FI_means) / diff(fisher_df$time_windows)

  tib <- tibble(
    # percent rate of change
    #percent_roc = 100 * diff(fisher_df$FI_means) / fisher_df[-nrow(fisher_df), ]$FI_means,
    abs_roc = abs(max(percent_roc)),
    max_pos_roc = max(percent_roc),
    min_neg_roc = min(percent_roc),
    timeof_max_pos_roc = fisher_df$time_windows[which.max(percent_roc)],
    timeof_min_neg_roc = fisher_df$time_windows[which.min(percent_roc)],
    timeof_abs_pos_roc = fisher_df$time_windows[which.max(abs(percent_roc))],
    max_run_of_roc_increase = max(rle(sign(percent_roc))[[1]][rle(sign(percent_roc))[[2]] == 1]),
    max_run_of_roc_decrease = max(rle(sign(percent_roc))[[1]][rle(sign(percent_roc))[[2]] == -1]),

    #slope:
    #slope = diff(fisher_df$FI_means) / diff(fisher_df$time_windows),
    max_slope = max(slope),
    min_slope =  min(slope),
    max_abs_slope =  max(abs(slope)),
    timeof_max_slope = fisher_df$time_windows[which.max(slope)],
    timeof_min_slope = fisher_df$time_windows[which.min(slope)],
    timeof_abs_slope = fisher_df$time_windows[which.max(abs(slope))],
    max_run_of_slope_increase = max(rle(sign(slope))[[1]][rle(sign(slope))[[2]] == 1]),
    max_run_of_slope_decrease = max(rle(sign(slope))[[1]][rle(sign(slope))[[2]] == -1]),

    run_of_increase = max(rle(sign(diff(fisher_df$FI_means)))[[1]][rle(sign(diff(fisher_df$FI_means)))[[2]] == 1]),
    run_of_decrease = max(rle(sign(diff(fisher_df$FI_means)))[[1]][rle(sign(diff(fisher_df$FI_means)))[[2]] == -1])
  )
  tib
}
#x1 <- slope_rates_func(fisher_df)


##### Changepoint parameters for normal distribution and empirical distribution
# Changepoints: http://members.cbio.mines-paristech.fr/~thocking/change-tutorial/RK-CptWorkshop.html
changepoint_func <- function(fisher_df){
    # multiple changepoints, penalty default MBIC. Changes in mean and variance
    cpt_pelt <- cpt.meanvar(qfish$FI_means, method = 'PELT')
    #cpts(cpt_pelt)
    # Mean parameters
    cpt_pelt_param <- param.est(cpt_pelt)
    no_cpt_pelt <- length(cpts(cpt_pelt))
    max_mean_cpt <- max(cpt_pelt_param[[1]])
    min_mean_cpt <- min(cpt_pelt_param[[1]])
    range_mean_cpt <- max_mean_cpt - min_mean_cpt
    max_pos_mean_cpt_in_seq <- max(diff(cpt_pelt_param[[1]]))
    max_neg_mean_cpt_in_seq <- min(diff(cpt_pelt_param[[1]]))
    # Variance parameters
    max_var_cpt <- max(cpt_pelt_param[[2]])
    min_var_cpt <- min(cpt_pelt_param[[2]])
    range_var_cpt <- max_var_cpt - min_var_cpt
    max_pos_var_cpt_in_seq <- max(diff(cpt_pelt_param[[2]]))
    max_neg_var_cpt_in_seq <- min(diff(cpt_pelt_param[[2]]))
    #plot(cpt_pelt)

    # multiple change points empirical distribution. General change in distribution
    cpt_pelt_np <-  cpt.np(qfish$FI_means, nquantiles = 4*log(length(qfish$FI_means)))
    #cpts(cpt_pelt_np)
    no_cpt_np <- length(cpts(cpt_pelt_np))
    #plot(cpt_pelt_np)
    tib <- tibble(
      no_cpt_pelt = no_cpt_pelt,
      max_mean_cpt = max_mean_cpt,
      min_mean_cpt = min_mean_cpt,
      range_mean_cpt = range_mean_cpt,
      max_pos_mean_cpt_in_seq = max_pos_mean_cpt_in_seq,
      max_neg_mean_cpt_in_seq = max_neg_mean_cpt_in_seq,
      max_var_cpt = max_var_cpt,
      min_var_cpt = min_var_cpt,
      range_var_cpt = range_var_cpt,
      max_pos_var_cpt_in_seq = max_pos_var_cpt_in_seq,
      max_neg_var_cpt_in_seq = max_neg_var_cpt_in_seq,
      no_cpt_np = no_cpt_np)
}
#x2 <- changepoint_func(fisher_df)


#####

magnitude_change_func <- function(fisher_df, wind = 10) {
  largest_change_lag_multi <- mapply(diff, lag = 1:(wind-1), MoreArgs = list(fisher_df$FI_means), SIMPLIFY = F)
  max_inc <- max(unlist(lapply(largest_change_lag_multi, max)))
  max_dec <- min(unlist(lapply(largest_change_lag_multi, min)))
  tibble(max_inc = max_inc,
         max_dec = max_dec)
}
# x3 <- magnitude_change_func(fisher_df)
# x3

##### State information from fishers, should mean states be within 1 sd of mean?
fisher_metrics <- function(fisher_df) {
  tib <- tibble(
    max_time_windows_in_median_states = max(rle(fisher_df$median_no_states)[[1]]),
    median_states_at_max_time_windows = rle(fisher_df$median_no_states)[[2]][max(rle(fisher_df$median_no_states)[[1]])],

    max_time_windows_in_mean_states = max(rle(fisher_df$mean_no_states)[[1]]),
    mean_states_at_max_time_windows = rle(fisher_df$mean_no_states)[[2]][max(rle(fisher_df$mean_no_states)[[1]])],
  )
  tib
}
#x4 <-fisher_metrics(fisher_df)
    # rate at longest run of increasing numbers
    # rate at longest run of decreasing numbers
    # slope at longest run of increasing numbers
    # slope at longest run of decreasing numbers




########## From astrophysics:
# https://isadoranun.github.io/tsfeat/FeaturesDocumentation.html#Range-of-a-cumulative-sum-$R_{cs}$
### https://arxiv.org/pdf/1506.00010.pdf        same as above^
# https://iopscience.iop.org/article/10.1088/0004-637X/735/2/68/meta#apj391424fd11
# https://academic.oup.com/mnras/article/400/4/1897/1079176
# https://iopscience.iop.org/article/10.1088/0004-637X/733/1/10/meta#apj387614s2
# https://www.aanda.org/articles/aa/pdf/2016/03/aa27188-15.pdf

# Following metrics are NOT included (mostly applicable to periodic features / waves):
### Stetson K, Stetson AC, Stetson J, Stetson L (https://academic.oup.com/mnras/article/400/4/1897/1079176, does provide single band stetson indices),
### slotted autocorrelation, Amplitude, Con (consecutive points), Anderson-Darling (normality), Linear trend, CAR (continuous time auto regressive model),
### color, ΨCS (cumulative sum), period fit, Lomb periodic features, Ψη (psi von neumann),

##### Range of cumulative sum / N * sd
rcs_func <- function(fisher_df) {
    s <- cumsum(fisher_df$FI_means - mean(fisher_df$FI_means)) / (length(fisher_df$FI_means) * sd(fisher_df$FI_means)) # Range of cumulative sum (over N*sd apparently...)
    tib <- tibble(rcs_max = max(s),
                  rcs_min = min(s),
                  rcs = rcs_max - rcs_min)
    tib
}
#x5 <- rcs_func(fisher_df)



##### Von_Neumann variance index η
von_neumann_func <- function(fisher_df) {
  tibble(
    # von newmann variance index
    # vi_von_neumann <- sum(diff(x)^2)/((length(x)-1)*var(x)) ## Eta. Removed, it is not invariant to time sampling: https://iopscience.iop.org/article/10.1088/0004-637X/735/2/68/meta#apj391424f10
    vi_von_neumann = mssd(fisher_df$FI_means)/var(fisher_df$FI_means)
    # vi_von_neumann <- 1/((length(x)-1)*var(x)) * sum(diff(x)^2) # https://github.com/isadoranun/FATS/blob/master/FATS/FeatureFunctionLib.py
)
}
#x6 <- von_neumann_func(fisher_df)
# Although η is a powerful index for quantifying variability characteristics of a time series, it does not take
# into account unequal sampling.
#Thus we use η^e, which is defined as:
##### Von Neumann variance index time invariant
von_neumann_time_invariant_func <- function(fisher_df) {
  wi <- 1 / (diff(fisher_df$time_windows) ^ 2)
  w_mean <- mean(wi)
  N = length(fisher_df$time_windows)
  sig2 <- var(fisher_df$FI_means)
  s1 <-  sum(wi * (diff(fisher_df$FI_means) ^ 2))
  s2 <- sum(wi)

  ne <- w_mean * (((length(fisher_df$time_windows)-1) - fisher_df$time_windows[1]) ^ 2)  * s1 / (sig2 * s2 * N ^ 2)
  #w_mean *  ((length(fisher_df$time_windows)-1 - fisher_df$time_windows[1])^2)  * s1 / (sig2 *s2)
  tibble(von_neumann_time_invariant = ne)
}
#x7 <- von_neumann_time_invariant_func(fisher_df)


##### Autocorrelation length
acf_length_func <- function(fisher_df) {
  lag = 0
  #print(lag)
  k = NA
  #print(k)
  # ac = acf(fisher_df$FI_means, lag.max = lag)
  # ac_lags <- ac[[1]]
  # k <- which(ac_lags < exp(-1))[1]-1 #acf function indexes from 0 so don't forget the -1 on the index!
  #k <- ac[[1]][2:length(ac[[1]])][which(ac[[1]] > exp(-1))[1]]
  #print(k)
  while (is.na(k)) {
    lag = lag + 100
    # print(lag)
    ac = acf(fisher_df$FI_means, lag.max = lag)
    #print(ac[[1]])
    ac_lags <- ac[[1]]
    k <-  which(ac_lags < exp(-1))[1] - 1 #acf function indexes from 0 so don't forget the -1 on the index!
    #k <- ac[[1]][2:length(ac[[1]])][which(ac[[1]] > exp(-1))[1]]
    k
    #print(k)
  }
  tibble(lag_k = k)
}
#x8 <- acf_length_func(fisher_df)


##### Median buffer range percentage
# ORIGINAL
# MedianBRP Median buffer range percentage
# def fit(self, data):
#   magnitude = data[0]
# median = np.median(magnitude)
# amplitude = (np.max(magnitude) - np.min(magnitude)) / 10
# n = len(magnitude)
#
# count = np.sum(np.logical_and(magnitude < median + amplitude,
#                               magnitude > median - amplitude))
#
# return float(count) / n
# R version:
median_BRP_func <- function(fisher_df) {
  med <- median(fisher_df$FI_means)
  amp <- (max(fisher_df$FI_means) - min(fisher_df$FI_means)) / 10
  #count_above <- fisher_df$FI_means > med + amp
  count_below <- fisher_df$FI_means < med + amp
  count <- table(count_below)[["TRUE"]] / length(fisher_df$FI_means)
  tibble(median_brp = count)
}
#x9 <- median_BRP_func(fisher_df)


##### Beyond 1 SD
# Percentage of points beyond one df from mean. Original uses weighted mean by photometric error
beyond_std_func <- function(fisher_df, sds = 1) {
  mu <- mean(fisher_df$FI_means)
  std <- sd(fisher_df$FI_means) * sds
  count_above <- fisher_df$FI_means > mu + std
  #count_below <- fisher_df$FI_means < med + amp
  percent_greater_than_sd <- table(count_above)[["TRUE"]] / length(fisher_df$FI_means)
  tibble(beyond_sd = percent_greater_than_sd)

}
#x10 <- beyond_std_func(fisher_df)



##### Fraction of increasing and decreasing first differences, called Pair slope trend, calculated for last 30 measurements in original
# ORIGINAL
# Considering the last 30 (time-sorted) measurements of source magnitude,
# the fraction of increasing first differences minus the fraction of
# decreasing first differences.
# """
#
# def __init__(self):
#     self.Data = ['magnitude']
#
# def fit(self, data):
#     magnitude = data[0]
#     data_last = magnitude[-30:]
#
#     return (float(len(np.where(np.diff(data_last) > 0)[0]) -
#             len(np.where(np.diff(data_last) <= 0)[0])) / 30)

#Fraction of increasing and decreasing first differences
# altered to total fraction rather than lase 30 points

frac_func <- function(fisher_df){
  frac <- table(diff(fisher_df$FI_means) > 0)
  frac_inc <- frac[["TRUE"]] / length(fisher_df$FI_means)
  frac_dec <- frac[["FALSE"]] / length(fisher_df$FI_means)
  (frac[["TRUE"]] - frac[["FALSE"]]) / length(fisher_df$FI_means)
  frac_ratio <- frac_dec / frac_inc
  tibble(frac_inc = frac_inc,
         frac_dec = frac_dec,
         frac_ratio = frac_ratio)
}
#x11 <- frac_func(fisher_df)



##### Small kurtosis. Kurtosis calculation for small samples. Different from the kurtosis() function from moments package used in summary statistics func
small_kurtosis_func <- function(fisher_df) {
  n <- length(fisher_df$FI_means)
  mean <- mean(fisher_df$FI_means)
  std <- sd(fisher_df$FI_means)
  S = sum(((fisher_df$FI_means - mean) / std) ** 4)
  c1 = (n * (n + 1)) / ((n - 1) * (n - 2) * (n - 3))
  c2 = (3 * (n - 1) ** 2) / ((n - 2) * (n - 3))
  result <- c1 * S - c2
  tibble(small_kurt = result)
}
#x12 <- small_kurtosis_func(fisher_df)



##### Flux percentile ratios:
### mid 20
flux_percentile_rario_mid_20 <- function(fisher_df) {
  f_60 <- quantile(fisher_df$FI_means, .60)[["60%"]]
  f_40 <- quantile(fisher_df$FI_means, .40)[["40%"]]
  f_05 <- quantile(fisher_df$FI_means, .05)[["5%"]]
  f_95 <- quantile(fisher_df$FI_means, .95)[["95%"]]

  f_40_60 <- f_60 - f_40
  f_5_95 <- f_95 - f_05
  f_mid_20 <- f_40_60 / f_5_95
  tibble(f_mid_20 = f_mid_20)
}
#f_mid_20 <- flux_percentile_rario_mid_20(fisher_df)


### mid 35:
flux_percentile_rario_mid_35 <- function(fisher_df) {
  f_325 <- quantile(fisher_df$FI_means, .325)[["32.5%"]]
  f_675 <- quantile(fisher_df$FI_means, .675)[["67.5%"]]
  f_05 <- quantile(fisher_df$FI_means, .05)[["5%"]]
  f_95 <- quantile(fisher_df$FI_means, .95)[["95%"]]

  f_325_675 <- f_675 - f_325
  f_5_95 <- f_95 - f_05
  f_mid_35 <- f_325_675 / f_5_95
  tibble(f_mid_35 = f_mid_35)
}
#f_mid_35 <- flux_percentile_rario_mid_35(fisher_df)


### mid 50:
flux_percentile_rario_mid_50 <- function(fisher_df) {
  f_25 <- quantile(fisher_df$FI_means, .25)[["25%"]]
  f_75 <- quantile(fisher_df$FI_means, .75)[["75%"]]
  f_05 <- quantile(fisher_df$FI_means, .05)[["5%"]]
  f_95 <- quantile(fisher_df$FI_means, .95)[["95%"]]

  f_25_75 <- f_75 - f_25
  f_5_95 <- f_95 - f_05
  f_mid_50 <- f_25_75 / f_5_95
  tibble(f_mid_50 = f_mid_50)
}
#f_mid_50 <- flux_percentile_rario_mid_50(fisher_df)


### mid 65:
flux_percentile_rario_mid_65 <- function(fisher_df) {
  f_175 <- quantile(fisher_df$FI_means, .175)[["17.5%"]]
  f_825 <- quantile(fisher_df$FI_means, .825)[["82.5%"]]
  f_05 <- quantile(fisher_df$FI_means, .05)[["5%"]]
  f_95 <- quantile(fisher_df$FI_means, .95)[["95%"]]

  f_175_825 <- f_825 - f_175
  f_5_95 <- f_95 - f_05
  f_mid_65 <- f_175_825 / f_5_95
  tibble(f_mid_65 = f_mid_65)
}
#f_mid_65 <- flux_percentile_rario_mid_65(fisher_df)


### mid 80:
flux_percentile_rario_mid_80 <- function(fisher_df) {
  f_10 <- quantile(fisher_df$FI_means, .10)[["10%"]]
  f_90 <- quantile(fisher_df$FI_means, .90)[["90%"]]
  f_05 <- quantile(fisher_df$FI_means, .05)[["5%"]]
  f_95 <- quantile(fisher_df$FI_means, .95)[["95%"]]

  f_10_90 <- f_90 - f_10
  f_5_95 <- f_95 - f_05
  f_mid_80 <- f_10_90 / f_5_95
  tibble(f_mid_80 = f_mid_80)
}
#f_mid_80 <- flux_percentile_rario_mid_80(fisher_df)


### percent difference flux percentile
percent_diff_flux_percentile <- function(fisher_df) {
  f_05 <- quantile(fisher_df$FI_means, .05)[["5%"]]
  f_95 <- quantile(fisher_df$FI_means, .95)[["95%"]]
  f_5_95 <- f_95 - f_05

  percent_diff <- f_5_95 / median(fisher_df$FI_means)
  tibble(percent_diff = percent_diff)
}
#percent_diff <- percent_diff_flux_percentile(fisher_df)


##### Q3−1
# Q3−1 is the difference between the third quartile, Q3, and the first quartile, Q1, of a raw light curve.
# Q1 is a split between the lowest 25% and the highest 75% of data.
# Q3 is a split between the lowest 75% and the highest 25% of data.
q31_func <- function(fisher_df) {
  q31 <- quantile(fisher_df$FI_means, .75)[["75%"]] - quantile(fisher_df$FI_means, .25)[["25%"]]
  tibble(q31 = q31)
}
#x <- q31_func(fisher_df)



### Needs correcting
# con_func <- function(fisher_df, wind) {
#   N = length(fisher_df$FI_means)
#   if (N < wind) {
#     return(0)
#   }
#
#   sigma = sd(fisher_df$FI_means)
#   m = mean(fisher_df$FI_means)
#   #count = 0
#
#   for (i in (1:(N - wind))) {
#     flag = c()
#     for (j in (1:wind)) {
#       # print(fisher_df$FI_means[i])
#       # print(fisher_df$FI_means[j])
#       if (fisher_df$FI_means[i] + fisher_df$FI_means[j] > m + 2 * sigma |
#           fisher_df$FI_means[i] + fisher_df$FI_means[j] < m - 2 * sigma) {
#         flag[i] = 1
#         #flag <- c(flag, flag[i])
#         #print(i)
#         #print(flag)
#       } else {
#         flag[i] = 0
#         #print(flag)
#       }
#       flag <- c(flag[i])
#     }
#     print(flag)
#   }
#   # break
#   # if flag:count = count + 1
#   # return count * 1.0 / (N - self.consecutiveStar + 1)
#   flag
# }


metrics_function <- function(fisher_df) {
  tib <- bind_cols(
    summary_stats <- summary_stats_func(fisher_df),
    slope_rate <- slope_rates_func(fisher_df),
    changepoints <- changepoint_func(fisher_df),
    mag_change <- magnitude_change_func(fisher_df),
    state_metrics <- fisher_metrics(fisher_df),
    rcs <- rcs_func(fisher_df),
    vi_vn <- von_neumann_func(fisher_df),
    vi_time_inv <- von_neumann_time_invariant_func(fisher_df),
    acf_length <- acf_length_func(fisher_df),
    med_brp <- median_BRP_func(fisher_df),
    beyond_std <- beyond_std_func(fisher_df),
    fracs <- frac_func(fisher_df),
    small_kurt <- small_kurtosis_func(fisher_df),
    percentile_ratio_mid_20 <- flux_percentile_rario_mid_20(fisher_df),
    percentile_ratio_mid_35 <- flux_percentile_rario_mid_35(fisher_df),
    percentile_ratio_mid_50 <- flux_percentile_rario_mid_50(fisher_df),
    percentile_ratio_mid_65 <- flux_percentile_rario_mid_65(fisher_df),
    percentile_ratio_mid_80 <- flux_percentile_rario_mid_80(fisher_df),
    percent_diff_percentile <- percent_diff_flux_percentile(fisher_df),
    q31 <- q31_func(fisher_df)
  )
  tib
}
# metrics_df <- metrics_function(qfish)
#
#
#
# qgam <- gam(qfish$FI_means ~ s(qfish$time_windows), method = "REML")
# qgam
# summary(qgam)
#
#
# ggplot(qfish, aes(x = time_windows, y = FI_means)) +
#   geom_line(colour = "blue") +
#   geom_line(aes(x = time_windows, y = FI_smth), colour = "red") +
#   geom_smooth(method = "gam", formula = y ~ s(x))
#
# ggplot(qfish, aes(x = time_windows, y = FI_means)) + geom_point() + geom_smooth(method = "loess")



