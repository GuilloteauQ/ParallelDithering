#!/usr/bin/env Rscript

opti_r <- function(h, w, p) {
    roots <- polyroot(z = c(h*w, 0, -h/2, 1))
    r <- roots[1]
    h / (p * 2^(as.integer(log2(h / (p * r)) + 0.5)))
}

opti_k <- function(w, p, r) {
    k <- w / (r * (p-0.5) + (p+0.5))
    2^(as.integer(log2(k) + 0.5))
}

pow2 <- function(n) {
    2^n
}

h <- 2^(27)
w <- 2^(27)
max_processors <- 128

df <- data.frame("processors" = seq(1, log2(max_processors)))
df$processors <- pow2(df$processors)
df

# df <- data.frame("row_block_size" = df[h %% (df$row_block_size * max_processors) == 0,])
df$h <- h
df$w <- w

for (i in 1:log2(max_processors)) {
    df$row_block_size[i] <- opti_r(df$h[i], df$w[i], df$processors[i])
    df$col_block_size[i] <- opti_k(df$w[i], df$processors[i], df$row_block_size[i])
}
# df$row_block_size <- opti_r(df$h, df$w, df$processors)
# df$col_block_size <- opti_k(df$w, df$processors, df$row_block_size)

write.csv(file = "../data/design_experiment_final.csv", df)


