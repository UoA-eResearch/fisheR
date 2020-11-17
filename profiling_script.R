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
plot(x = time_df$run, y = time_df$time)

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
plot(x = time_df2$run, y = time_df2$time)


#####################################################################################
###################################   profiling   ###################################
#####################################################################################

small_prof <- profvis::profvis({
  small_test_fi <- fisheR::fisher(small_test_df, w_size = 8, w_incre = 1)
})


medium_prof <- profvis::profvis({
  medium_test_fi <- fisheR::fisher(medium_test_df, w_size = 8, w_incre = 1)
})


#####################################################################################
################################## version testing ##################################
#####################################################################################
# THIS SCRIPT IS AN ATTEMPT TO TEST TWO VERSIONS OF fisheR
# THE SCRIPT IS NOT WORKING YET, DESPITE SPECIFYING THE COMMIT CODE, PACKAGE VERSION IS NOT CHANGING
# library(testthat)
# 
# set.seed(1984)
# x <- data.frame(time = 1:1000,
#                 v1 = c(rnorm(500, 1, 0.1), rnorm(200, 10, 1),  rnorm(100, 100, 50), rnorm(200, 1, 0.1)),
#                 v2 = c(rnorm(500, 1, 0.1), rnorm(200, 10, 1),  rnorm(100, 100, 50), rnorm(200, 1, 0.1)),
#                 v3 = c(rnorm(500, 1, 0.1), rnorm(200, 10, 1),  rnorm(100, 100, 50), rnorm(200, 1, 0.1)),
#                 v4 = c(rnorm(500, 1, 0.1), rnorm(200, 10, 1),  rnorm(100, 100, 50), rnorm(200, 1, 0.1)),
#                 v5 = c(rnorm(500, 1, 0.1), rnorm(200, 10, 1),  rnorm(100, 100, 50), rnorm(200, 1, 0.1)),
#                 v6 = c(rnorm(500, 1, 0.1), rnorm(200, 10, 1),  rnorm(100, 100, 50), rnorm(200, 1, 0.1))
# )
# 
# test_that("string required?", {
#   
#   devtools::install_github("UoA-eResearch/fisheR@09f9324460f9797114a6a8e37b1248a90af82277", force = T)
#   fish_x_old <- fisheR::fisher(x, w_size = 8)
#   
#   devtools::install_github("UoA-eResearch/fisheR@d7ca554576f6947c3c819744d6f0082430e7c088", force = T)
#   fish_x_new <- fisheR::fisher(x, w_size = 8)
#   
#   expect_identical(fish_x_old, fish_x_new)
#   expect_equal(fish_x_old, fish_x_new)
# })





