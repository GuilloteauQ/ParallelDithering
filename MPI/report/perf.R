#!/usr/bin/env Rscript

library(ggplot2)
library(tidyr)
library(dplyr)

data_file <- "../experiments/data/output_rpi.csv"

df <- read.csv(data_file, header = TRUE, sep = " ")

df$speedup <- df$seq_time / df$par_time

head(df)

# We need to compute the 95% interval for each quadruple (n, p, k, r)

df %>% group_by(c(w, p, k, r)) %>% summarise(mean = mean(speedup), sd sd(speedup), min = min(speedup), max = max(speedup))

df$conf_inter_inf <- df$mean - 2 * df$sd
df$conf_inter_sup <- df$mean + 2 * df$sd

head(df)
