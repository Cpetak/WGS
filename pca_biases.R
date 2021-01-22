#checking for biases and unwanted patterns

#By location
C <- as.matrix(read.table("pcangsd_covmatrix_noout.cov"))
ids <- read.table("pca_pops_noout.txt")
e <- eigen(C)

#ggplot
library(ggplot2)
library(tidyverse)
df <- data.frame(pop = ids$V1, PC1 = e$vectors[,3], PC2 = e$vectors[,4])

df= rownames_to_column(df)
ggplot(df, aes(x=PC1, y=PC2, fill=pop)) +
  geom_point(size=3, shape=21) +
  theme_bw()

#By pH variability (L=low variation, H=high variation)
ids <- read.table("ph_pca_noout.txt")
e <- eigen(C)

#ggplot
df <- data.frame(pop = ids$V1, PC1 = e$vectors[,1], PC2 = e$vectors[,2])

df= rownames_to_column(df)
ggplot(df, aes(x=PC1, y=PC2, fill=pop)) +
  geom_point(size=3, shape=21) +
  theme_bw()

#By north vs south vs middle
#Same as by east vs west vs middle as more north more west, more south more east
ids <- read.table("northsouthmiddle_pca_noout.txt")
e <- eigen(C)

#ggplot
df <- data.frame(pop = ids$V1, PC1 = e$vectors[,1], PC2 = e$vectors[,2])

df= rownames_to_column(df)
ggplot(df, aes(x=PC1, y=PC2, fill=pop)) +
  geom_point(size=3, shape=21) +
  theme_bw()

#By coverage
covdata <- read.table("covs.txt")
covdata$V2 <- ifelse(covdata$V1<6.5, "little", "lot")
hist(covdata$V1)

ids <-covdata
e <- eigen(C)

#ggplot
df <- data.frame(pop = ids$V2, PC1 = e$vectors[,1], PC2 = e$vectors[,2])

df= rownames_to_column(df)
ggplot(df, aes(x=PC1, y=PC2, fill=pop)) +
  geom_point(size=3, shape=21) +
  theme_bw()

 