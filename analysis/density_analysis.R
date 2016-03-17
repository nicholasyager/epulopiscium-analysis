# Cleanup old variables and load backages if necessary
rm(list=ls())
par(mfrow=c(1,1))
packages <- c('vioplot')
for (index in 1:length(packages)){
  package = packages[index]
  exists = require(package,character.only = T)
  if (exists == FALSE) {
    install.packages(package)
    require(package,character.only = T)
  }
}

#data <- read.csv("densities.csv")

data1 <- read.csv("A_Early stages of daughter cell formation//densities_raw.csv")
data2 <- read.csv("B_Fully engulfed, small daughter cells//densities_raw.csv")
data3 <- read.csv('C_Fully engulfed, full length daughter cells//densities_raw.csv')

data <- rbind(data1,data2)
data <- rbind(data,data3)

stageStrings = c("A","B","C")

data <- read.csv("aggregate_data.csv")

# Separate data by life cycle stage
stages = c()
for (row in 1:nrow(data)) {
  stageString = strsplit(as.character(data$analysis[row]),"_")[[1]][1]
  stage = which(stageString == stageStrings)
  stages = append(stages, stage)
}
data <- cbind(data, stage = stages)
stages = split(data, data$stage)


for (index in 1:3) {
print(paste(index,"-- mean:", 1/mean(stages[[index]]$density),
            "sd:",1/sd(stages[[index]]$density)))
}

vioplot(subset(data, stage=="1")$density,
        subset(data, stage=="2")$density,
        subset(data, stage=="3")$density,
        col="lightblue",
        names = stageStrings)
title( ylab = "Chromosomes per Cubic Micrometer",
       xlab = "Life Cycle Stage",
       main = "Distribution of Chromosome Density of Euplopiscium by Stage")

kruskal.test(density~stage,data=data)

library("PMCMR")
posthoc.kruskal.nemenyi.test(density~stage,data=data, dist="Tukey")

# Boxplot of the data by stage

boxplot(stages[[1]]$density,
        stages[[2]]$density,
        stages[[3]]$density,
        ylab = "Chromosome per Cubic Micrometers",
        xlab = "Life Cycle Stage",
        main = "Distribution of Chromosome Density of Euplopiscium by Stage",
        names = stageStrings)
abline(h=1.9,col="blue")

vioplot(stages[[1]]$density,
        stages[[2]]$density,
        stages[[3]]$density,
        col="lightblue",
        names = stageStrings)
title( ylab = "Chromosomes per Cubic Micrometer",
       xlab = "Life Cycle Stage",
       main = "Distribution of Chromosome Density of Euplopiscium by Stage")

# Number of chromosomes vs. volume colored by stage
plot(data$volume, data$chromosomes,col=rainbow(max(data$stage))[data$stage],pch=19,
     ylab= "Number of Chromosomes",
     xlab= "Cell Volume (um^3)",
     main= "Number of Chromosomes vs. Volume by Life Cycle Stage")
legend("topleft",legend=stageStrings,fill=rainbow(max(data$stage))[as.numeric(names(stages))])
model <- lm(data$chromosomes~data$volume-1)
summary(model)
abline(model,col="red")
text(120000,1000,paste("Density:",round(model$coefficients[1],4),"genomes per um^3."))
cat(paste("Density:",model$coefficients[1],"genomes per um^3."))

write.csv(data,"analysis.csv",row.names=F)
