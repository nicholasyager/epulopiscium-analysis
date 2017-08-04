library(plyr)
library(dplyr)
library(reshape2)

trials <- read.csv("densities.csv")
colnames(trials) <- c("analysis", 'cell', 'chromosomes', 'volume', 'density')
trials$daughter = sapply(strsplit(as.character(trials$analysis),"offspring"), "[[", 2)
trials$analysis = sapply(strsplit(as.character(trials$analysis),"_offspring"), "[[", 1)

trials.counts <- dcast(trials, analysis ~ daughter, value.var = 'chromosomes')
trials.counts$diff <- trials.counts$`1` - trials.counts$`2`

hist(trials.counts$diff)

trials.counts <- dcast(trials, analysis ~ daughter, value.var = 'density')
trials.counts$diff <- trials.counts$`1` - trials.counts$`2`

hist(trials.counts$diff)
plot(density(trials.counts$diff))

qqnorm(trials.counts$diff)
qqline(trials.counts$diff)

mean(trials.counts$diff)
wilcox.test(trials.counts$diff, alternative = "two.sided", mu = 0)
