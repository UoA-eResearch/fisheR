library(fisheR)
library(profvis)

small_test_df <- data.frame(time = 1:100,
                            v1 = rnorm(100),
                            v2 = rnorm(100, 1),
                            v3 = rnorm(100,2,2),
                            v4 = rnorm(100),
                            v5 = rnorm(100,1,2))




medium_test_df <- data.frame(time = 1:1000,
                            v1 = rnorm(1000),
                            v2 = rnorm(1000, 1),
                            v3 = rnorm(1000,2,2),
                            v4 = rnorm(1000),
                            v5 = rnorm(1000),
                            v6 = rnorm(1000, 1),
                            v7 = rnorm(1000,2,2),
                            v8 = rnorm(1000),
                            v9 = rnorm(1000,1,2))


#####################################################################################

#####################################################################################
###################################   plotting    ###################################
#####################################################################################

#### To plot the time DECREASE using larger window increments
timing_function2 <- mapply(
  fisher,
  w_incre = c(1:5),
  MoreArgs = list(small_test_df),
  SIMPLIFY = F
)

# taking timing manually from printed output
time_df <- data.frame(time = c(14.82, 7.60, 4.14, 2.77, 2.13),
                      run = 1:5)
plot(x = time_df$time, y = time_df$run)

###################################
#### To plot the time INCREASE using larger window sizes
timing_function2 <- mapply(
  fisher,
  w_size = c(8:12),
  MoreArgs = list(small_test_df),
  SIMPLIFY = F
)

# taking timing manually from printed output
time_df2 <- data.frame(time = c(13.49, 13.38, 14.08, 15.74, 16.96),
                      run = 1:5)
plot(x = time_df2$time, y = time_df2$run)


#####################################################################################
###################################   profiling   ###################################
#####################################################################################

small_prof <- profvis::profvis({
  small_test_fi <- fisheR::fisher(small_test_df, w_size = 8, w_incre = 1)
})


medium_prof <- profvis::profvis({
  medium_test_fi <- fisheR::fisher(medium_test_df, w_size = 8, w_incre = 1)
})




