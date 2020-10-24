library(ggplot2)

data <- read.csv("all_mapping_stats.csv", header = TRUE)
data$ID <- row.names(data)

plot(data$Paired.Num.reads)

ggplot(data, aes(x=ID, y=Paired.Num.reads)) + 
  geom_point() +
  ylim(0, 1) +
  xlab("Samples") +
  ylab("Prop. correctly mapped (Paired on same chr)") +
  theme(axis.text.x = element_text(size = 5, angle = 90))

ggplot(data, aes(x=ID, y=Coverage)) + 
  geom_point() +
  ylim(0, 10) +
  xlab("Samples") +
  ylab("Average coverage") +
  theme(axis.text.x = element_text(size = 5, angle = 90))

ggplot(data, aes(x=ID, y=Mapped.Tot)) + 
  geom_point() +
  ylim(0, 1) +
  xlab("Samples") +
  ylab("Mapped reads/Total reads") +
  theme(axis.text.x = element_text(size = 5, angle = 90))
