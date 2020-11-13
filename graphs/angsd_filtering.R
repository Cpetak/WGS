df <- read.csv("angsd_filtering.csv", header=TRUE)

library(ggplot2)
library(ggpubr)
theme_set(theme_pubr())

ggplot(df, aes(x = Pop, y = Freq))+
  geom_bar(
    aes(fill = Filtering), stat = "identity", color = "white",
    position = position_dodge(0.9)
  )+
  facet_wrap(~Site) + 
  fill_palette("jco") +
  xlab("Population") +
  ylab("Number of sites") +
  scale_y_continuous(limits = c(0,5e+08), expand = c(0, 0))
