#! /usr/bin/env Rscript


setwd("../..") # Root of the project
experiments <- read.csv(file = "experiments/data/design_experiment.csv")

valid_experiments <- experiments[(experiments$image_size / experiments$row_block_size) %% experiments$processors == 0,]


# Build the project
system("make > /dev/null")

mandrill_size <- 512

commands <- paste("mpirun -np", valid_experiments$processors, "./dithering_bw", "Image/mandrill.pgm",  valid_experiments$col_block_size, valid_experiments$row_block_size,(valid_experiments$image_size / mandrill_size) , sep = " ")

begin <- as.numeric(commandArgs(TRUE)[1])

end <- as.numeric(commandArgs(TRUE)[2])


for (i in begin:end){
    system(commands[i])
}
