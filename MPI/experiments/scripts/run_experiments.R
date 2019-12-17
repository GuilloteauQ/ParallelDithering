#! /usr/bin/env Rscript


setwd("../..") # Root of the project
experiments <- read.csv(file = "experiments/data/design_experiment.csv")

valid_experiments <- experiments[(experiments$image_size / experiments$row_block_size) %% experiments$processors == 0,]

head(valid_experiments)

# Build the project
system("make")

mandrill_size <- 512
tmp_image <- "Images/tmp.pgm"

duplicates_commands <- paste("./duplicate_image", "Images/mandrill.pgm", tmp_image, (valid_experiments$image_size / mandrill_size), sep = " ")

commands <- paste(duplicates_commands, ";", "mpirun -np", valid_experiments$processors, "./dithering_bw", tmp_image, valid_experiments$col_block_size, valid_experiments$row_block_size, sep = " ")


head(commands)
# system(commands[1])
