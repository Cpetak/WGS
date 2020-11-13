C <- as.matrix(read.table("pcangsd_covmatrix_noout.cov"))
ids <- read.table("~/Downloads/pca_pops_noout.txt")
e <- eigen(C)
# base R
plot(e$vectors[,1:2],xlab="PC1",ylab="PC2", bg=ids$V1, pch=21)
#ggplot
library(ggplot2)
df <- data.frame(pop = ids$V1, PC1 = e$vectors[,1], PC2 = e$vectors[,2])

#df= rownames_to_column(df)
ggplot(df, aes(x=PC1, y=PC2, fill=pop)) +
  geom_point(size=3, shape=21) +
  theme_bw()

for (s in e$values) {
  print(s / sum(e$values)*100)
}


caps = subset(df, pop == "CAP") 
caps = rownames_to_column(caps)

ggplot(fogs, aes(x=PC1, y=PC2, fill=rowname)) +
  geom_point(size=3, shape=21) +
  theme_bw()


