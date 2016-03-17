library("ggplot2")

stage_directories <- dir(pattern = "A|B|C")

for (stageDir in stage_directories) {
    stage = substr(stageDir, 1,1)

    analysisDirectories <- dir(path = stageDir, pattern = "_analysis")
    for (analysisDir in analysisDirectories) {
        analysisParts <- strsplit(analysisDir, "_")
        analysis <- paste(analysisParts[[1]][1:3],collapse="_")
        points <- read.csv(paste(stageDir,"/",analysisDir,"/points.csv",sep=""))
        p <- ggplot(subset(subset(points, dist > 0),dist < 6), aes(x,y, colour=dist)) +
            geom_point(size=rel(2.5), alpha=0.75) +
            xlim(0,108) +
            ylim(0,108) +
            scale_colour_gradientn(colours=rainbow(10), name="Avg.\nDistance (μm)") +
            xlab("x (μm)")+
            ylab("y (μm)") +
            theme( panel.background = element_rect(fill = 'white', colour = 'white'))
        ggsave(filename = paste("plots/",analysis,".png",sep=""), plot=p)
    }
}
