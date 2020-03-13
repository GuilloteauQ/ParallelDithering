#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)

data_file <- "../experiments/data/output_rpi.csv"

df <- read.csv(data_file, header = TRUE, sep = " ")

df$speedup <- df$seq_time / df$par_time

# We need to compute the 95% interval for each quadruple (n, p, k, r)

df <- df %>% group_by(w, p, k, r) %>% summarise(mean = mean(speedup), sd = sd(speedup), min = min(speedup), max = max(speedup))

df$conf_inter_inf <- df$mean - 2 * df$sd
df$conf_inter_sup <- df$mean + 2 * df$sd

head(df)

# 
# df_p <- df %>% group_by(w, p) %>% summarise(mean = mean(speedup), sd = sd(speedup), min = min(speedup), max = max(speedup))
# 
# df_p$conf_inter_inf <- df_p$mean - 2 * df_p$sd
# df_p$conf_inter_sup <- df_p$mean + 2 * df_p$sd
# 
# df_p
# 
# ggplot(data = df_p, aes(x = w, y = mean)) +
#     theme_bw() +
#     geom_errorbar(aes(ymin = conf_inter_inf, ymax = conf_inter_sup), color = "black") +
#     geom_point()
