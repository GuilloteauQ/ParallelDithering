#!/usr/bin/env Rscript

opti_r <- function(h, w, p) {
    roots <- polyroot(z = c(h*w, 0, -h/2, 1))
    r <- roots[1]
    h / (as.integer(h/(r*p) + 1) * p)
}

opti_k <- function(w, p, r) {
    k <- w / (r * (p-0.5) + (p+0.5))
    # only ok because image is of size 2^n
    2^(as.integer(log2(k) + 0.5))
    # inc <- 1
    # sign <- 1
    # while ((w / k)%%1 != 0) {
    #     if (sign == 1) {
    #         k <- k + inc
    #     } else {
    #         k <- k - inc
    #     }
    #     sign <- 1 - sign
    #     inc <- inc + 1
    # }
    # k
}

number_of_runs <- 200
h <- 2^(13)
w <- 2^(13)
max_block_size <- 10
max_processors <- 16

df <- data.frame("row_block_size" = seq(2, h))

df <- data.frame("row_block_size" = df[h %% (df$row_block_size * max_processors) == 0,])
df$h <- h
df$w <- w
df$processors <- max_processors
df$col_block_size <- opti_k(df$w, df$p, df$row_block_size)

opti <- data.frame(h, w, max_processors)
names(opti) <- c("h", "w", "processors")
opti$row_block_size <- opti_r(h, w, max_processors)
opti$col_block_size <- opti_k(w, max_processors, opti$row_block_size)


df <- rbind(df, opti)

write.csv(file = "../data/design_experiment_opti.csv", df)


