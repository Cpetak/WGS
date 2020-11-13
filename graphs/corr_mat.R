library(metaSEM)

my.full <- readFullMat("matrix_tail.out") 

corr <- cov2cor(my.full[[1]])

library('plot.matrix')
plot(corr)

par(mar=c(4.1, 4.1, 4.1, 4))
plot(corr, main="Correlation matrix output of Bayenv")
