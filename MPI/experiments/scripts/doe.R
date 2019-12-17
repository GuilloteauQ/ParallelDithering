#! /usr/bin/env Rscript

if (!require("DoE.wrapper")) install.packages("DoE.wrapper");

library("DoE.wrapper")

number_of_runs <- 2000
min_image_size <- 15
max_image_size <- 20
max_block_size <- 10
max_processors <- 6

data_frame <- lhs.design(
                         type = "maximin",
                         nruns = number_of_runs,
                         nfactors = 4,
                         digits = 0,
                         seed = 17440,
                         factor.names = list(image_size = c(min_image_size, max_image_size),
                                             row_block_size = c(1, max_block_size),
                                             col_block_size = c(1, max_block_size),
                                             processors = c(1, max_processors)))

pow2 <- function(n) {
    2 ^ n
}

data_frame$image_size <- pow2(data_frame$image_size)
data_frame$col_block_size <- pow2(data_frame$col_block_size)
data_frame$row_block_size <- pow2(data_frame$row_block_size)
data_frame$processors <- pow2(data_frame$processors)

head(data_frame)

write.csv(file = "../data/design_experiment.csv", data_frame)
