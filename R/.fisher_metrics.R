fisher_metrics function(fisher_df) {

    max_fish = max(fisher_df$FI_means)
    min_fish = min(fisher_df$FI_means)
    range_fish = max_fish - min_fish
    mid_point = max_fish / 2
    mu = mean(fisher_df$FI_means)
    std = sd(fisher_df$FI_means)
    variance = var(fisher_df$FI_means)
    percent_roc = 100 * diff(fisher_df$FI_means) / fisher_df[-nrow(fisher_df), ]$FI_means
    max_pos_roc = max(percent_roc)
    min_pos_roc = min(percent_roc)
    timeof_max_pos_roc = fisher_df$time_windows[which.max(percent_roc)]
    timeofmin_pos_roc = fisher_df$time_windows[which.min(percent_roc)]



  longest run of increasing numbers
  longest run of decreasing numbers
  rate at longest run of increasing numbers
  rate at longest run of decreasing numbers

  longest run of increasing roc
  longest run of decreasing roc
  number of trend reversals

  max diff from ave
  greatest change on y axis... (which greatest ROC at longest run pos/neg?)
)
  # In Fisher
  number of states
  max/min probability (at time)(most ordered/disordered point or window)
  probability over time?
  longest period of single state
  shortest period of single state
}


sample_data <- read.csv("C:/Users/qase352/Dropbox/QuinnAsenaPhD/R/fishers_information/fi_ahmad_2016/FI_Scripts_v2.00/sample_data_FI.csv")
head(sample_data)

qdata <- readRDS("C:/Users/qase352/Desktop/1_complete_output.RData")
sub_qdata <- as.data.frame(qdata[[1]][[4]][1:500, 1:50])
time <- 1:nrow(sub_qdata)
sub_qdata <- cbind(time,sub_qdata)
head(sub_qdata)
qfish <- fisher(sub_qdata, display_plot = T)
plot(qfish$time_windows, qfish$FI_means, type="l", col="blue", xlab = "Time Step", ylab = "Fisher Information")
lines(qfish$time_windows, qfish$FI_smth, type="l", col="red")

qgam <- gam(qfish$FI_means ~ s(qfish$time_windows), method = "REML")
qgam

ggplot(qfish, aes(x = time_windows, y = FI_means)) +
  geom_line(colour = "blue") +
  geom_line(aes(x = time_windows, y = FI_smth), colour = "red") +
  geom_smooth(method = "gam", formula = y ~ s(x))

ggplot(qfish, aes(x = time_windows, y = FI_means)) + geom_point() + geom_smooth(method = "loess")


Rolling window what is the max/min/slope/trend/breakpoints




## Nile data with one breakpoint: the annual flows drop in 1898
## because the first Ashwan dam was built
data("Nile")
plot(Nile)

fs.nile <- Fstats(Nile ~ 1)
plot(fs.nile)
breakpoints(fs.nile)
lines(breakpoints(fs.nile))

## or
bp.nile <- breakpoints(Nile ~ 1)
summary(bp.nile)

## the BIC also chooses one breakpoint
plot(bp.nile)
breakpoints(bp.nile)

## fit null hypothesis model and model with 1 breakpoint
fm0 <- lm(Nile ~ 1)
fm1 <- lm(Nile ~ breakfactor(bp.nile, breaks = 1))
plot(Nile)
lines(ts(fitted(fm0), start = 1871), col = 3)
lines(ts(fitted(fm1), start = 1871), col = 4)
lines(bp.nile)

## confidence interval
ci.nile <- confint(bp.nile)
ci.nile
lines(ci.nile)


#####
bp.qfish <- breakpoints(qfish$FI_means ~ 1, h = 8)
summary(bp.qfish)
plot(bp.qfish)
breakpoints(bp.qfish)

fm0 <- lm(qfish$FI_means ~ 1)
fm1 <- lm(qfish$FI_means ~ breakfactor(bp.qfish))
plot(qfish$FI_means, type = 'l')
lines(ts(fitted(fm0), start = 0), col = 3)
lines(ts(fitted(fm1), start = 0), col = 4)
lines(bp.qfish)






