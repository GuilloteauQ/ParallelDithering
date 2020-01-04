#! /usr/bin/env Rscript

if (!require("DoE.wrapper")) install.packages("DoE.wrapper");

library("DoE.wrapper")

number_of_runs <- 200
min_image_size <- 12
max_image_size <- 20
max_processors <- 6

data_frame <- lhs.design(
                         type = "maximin",
                         nruns = number_of_runs,
                         nfactors = 2,
                         digits = 0,
                         seed = 17440,
                         factor.names = list(image_size = c(min_image_size, max_image_size),
                                             processors = c(1, max_processors)))

pow2 <- function(n) {
    2 ^ n
}

data_frame$image_size <- pow2(data_frame$image_size)
data_frame$processors <- pow2(data_frame$processors)
data_frame$col_block_size <- data_frame$image_size / (2 * data_frame$processors)
data_frame$row_block_size <- pow2(as.integer(log2(((data_frame$image_size / data_frame$col_block_size) - (data_frame$processors + 0.5)) / (data_frame$processors - 0.5))))

data_frame <- data_frame[(data_frame$image_size / data_frame$processors) %% data_frame$row_block_size == 0,]

head(data_frame)

write.csv(file = "../data/design_experiment.csv", data_frame)

# number_of_runs <- 200
# min_image_size <- 9
# max_image_size <- 13
# max_block_size <- 8
# max_processors <- 3
# 
# data_frame_rpi <- lhs.design(
#                          type = "maximin",
#                          nruns = number_of_runs,
#                          nfactors = 4,
#                          digits = 0,
#                          seed = 17440,
#                          factor.names = list(image_size = c(min_image_size, max_image_size),
#                                              row_block_size = c(1, max_block_size),
#                                              col_block_size = c(1, max_block_size),
#                                              processors = c(1, max_processors)))
# 
# data_frame_rpi$image_size <- pow2(data_frame_rpi$image_size)
# data_frame_rpi$col_block_size <- pow2(data_frame_rpi$col_block_size)
# data_frame_rpi$row_block_size <- pow2(data_frame_rpi$row_block_size)
# data_frame_rpi$processors <- pow2(data_frame_rpi$processors)
# 
# head(data_frame_rpi)
# 
# write.csv(file = "../data/design_experiment_rpi.csv", data_frame_rpi)
# 
