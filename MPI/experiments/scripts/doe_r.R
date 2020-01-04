#!/usr/bin/env Rscript

pow2 <- function(n) {
    2 ^ n
}
number_of_runs <- 200
image_size <- pow2(13)
max_image_size <- 15
max_block_size <- 10
max_processors <- 16
# 
# data_frame <- lhs.design(
#                          type = "maximin",
#                          nruns = number_of_runs,
#                          nfactors = 3,
#                          digits = 0,
#                          seed = 17440,
#                          factor.names = list(image_size = c(min_image_size, max_image_size),
#                                              col_block_size = c(1, max_block_size),
#                                              processors = c(1, max_processors)))
# 
# data_frame$image_size <- image_size
# data_frame$row_block_size <- 1
# data_frame$processors <- pow2(data_frame$processors)
# # data_frame$col_block_size <- data_frame$image_size / (2 * data_frame$processors)
# 
# data_frame$col_block_size <- seq(2, image_size) 

df <- data.frame("row_block_size" = seq(2, image_size))

df <- data.frame("row_block_size" = df[image_size %% df$row_block_size == 0,])
df$image_size <- image_size
df$col_block_size <- image_size / (2 * max_processors)
df$processors <- max_processors


write.csv(file = "../data/design_experiment_r.csv", df)


