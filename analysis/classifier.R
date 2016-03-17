# Cleanup old variables and load backages if necessary
rm(list=ls())
par(mfrow=c(1,1))
packages <- c('rgl','cluster','plotrix','geometry','mclust','MASS')
for (index in 1:length(packages)){
    package = packages[index]
    exists = require(package,character.only = T)
    if (exists == FALSE) {
        install.packages(package)
        require(package,character.only = T)
    }
}

dot <- function(a,b) {
    product <- 0
    for (index in 1:length(a)) {
        product = product + (a[index] * b[index])
    }
    return(product)
}

vectorLength <- function(v) {
    sumOfSquares <- 0
    for (index in 1:length(v)) {
        sumOfSquares = sumOfSquares + v[index]^2
    }
    return(sqrt(sumOfSquares))
}

projection <- function(u,v) {
    # Compute the projection of u on v.
    return((dot(u,v)/vectorLength(u)^2)*u)
}

# Function declarations
distanceFromLine <- function(a,b,p) {
    # Calculate the distance from a point to a line
    normal = c(-a[2],a[1])
    delta = b - p
    displacementVector = projection(normal, delta)
    distance = vectorLength(displacementVector)
    return(distance)
}

normalize <- function(data) {
    return ((data - min(data))/(max(data) - min(data)))
}

distanceFromLine(c(1,1), c(0,0), c(1,0))

directories <- dir(pattern="_analysis")

for (dirIndex in 1:length(directories)) {

    directory = directories[dirIndex]

    data <- read.csv(paste(directory,"/points.csv",sep=""))
    data <- subset(data, dist > 0)
    # Look for duplicates and sort based on brightness.
    if (file.exists(paste(directory,"/classification",sep=""))){
        classification <- read.csv(paste(directory,"/classification",sep=""))
    } else {
        cat("Classification missing! Manual input required.\n")
        plot(data[,1:3],col=rainbow(100)[normalize(data$dist)*100])
        cat("Classification:")
        class = as.numeric(readline())
        cat("Number of clusters:")
        nClusters = as.numeric(readline())
        cat("Allignment necessary (0/1):")
        align = as.numeric(readline())
        classification <- data.frame(nClusters=nClusters,class=class,align=align)
        cat("Saving classification.")
        write.csv(classification,paste(directory,"/classification",sep=""),row.names=F)
    }

    par(mfrow=c(1,2))
    ref = hist(data$dist, breaks=500)
    maxDist <- 6
    abline(v=median(data$dist),col="red")
    abline(v=maxDist,col="blue")
    data <- data[which(data$dist < maxDist),]
    plot(data[,1:2],col=rainbow(100)[normalize(data$dist)*100])
par(mfrow=c(1,1))

    # Choose a clustering technique based on classification and number of clusters
    # Class 1: Two linear daughter cells.
    # Class 2: Two parallel daughter cells in the XY-plane.
    # Class 3: Two parallel daughter cells parallel to the X-axis.
    # Class 4: Two parallel daughter cells perpendicular to the XY-plane.
    # Class 5: One mother cell.

    dim = list(1:2,2:3,c(1,3),1:2,1:3)
    nClust <- classification$nClusters
    coordinates <- data[,dim[[classification$class]]]

    #nClust =

    data = data[data$dist <= maxDist,]
    coordinates <- data[,dim[[classification$class]]]

    if (classification$class == 1) {

        # Remove mother cell chromosomes

        plot(coordinates)

        fit <- Mclust(coordinates, G = nClust)
        data <- cbind(data,cluster = fit$classification)

    } else if (classification$class == 2) {
        # If Class 2, examine rotation about the Z axis. Correct if necessary
        if (classification$align == 1){
            par(mfrow=c(2,2))

            plot(data[,1:2])

            model <- lm(data$y~data$x)
            abline(model,col="red")

            a = c(1,model$coefficients[2])
            b = c(median(data$x),median(data$y))

            bestA <- a
            bestB <- b
            bestHits <- Inf
            t = seq(-100,100,1)
            sampleSet <- sample(1:nrow(data),1000)
            governor = 1
            for (step in 1:20) {

                a <- bestA
                b <- bestB

                # Modify the line.
                if (runif(1) <=0.5) {
                    a[1] <- a[1] + runif(1,-0.1,0.1) * governor
                } else {
                    a[2] <- a[2] + runif(1,-0.1,0.1) * governor
                }

                # Calculate distance from line
                distances <- c()
                for (index in sampleSet) {
                    distances <- append(distances, distanceFromLine(a, b, data[index,1:2]))
                }
                hits = which(distances < 0.75)
                lines(a[1]*t+b[1], a[2]*t+b[2], col="red")
                points(b[1],b[2],col="red",pch=19)
                print(paste(step,length(hits)))

                if (length(hits) < bestHits) {
                    bestA <- a
                    governor <- (100-step)/100
                    bestHits <- length(hits)
                }
            }

            a <- bestA


            angle = atan(a[2]/a[1])
            intercept = b[2] + a[2]*(-b[1]/a[1])
            rotData = matrix(0,nrow=nrow(data),ncol=3)
            abline(b=a[2]/a[1],a=intercept,col="red")
            # Attempt to rotate the points about an axis.
            for (row in 1:nrow(data)) {

                rotData[row,] <- c(data[row,1] * cos(-angle) - data[row,2] * sin(-angle) ,
                                   data[row,1] * sin(-angle) + data[row,2] * cos(-angle),data[row,3])
            }
            plot(rotData[,1:2])
            plot(coordinates)
            coordinates <- rotData[,dim[[classification$class]]]
            plot(coordinates)
            par(mfrow=c(1,1))
        }

        fit <- kmeans(coordinates, centers=nClust)
        data <- cbind(data,cluster = fit$cluster)

    } else if (classification$class == 3) {
        # If Class 2, examine rotation about the Z axis. Correct if necessary
        if (classification$align == 1){
            par(mfrow=c(2,2))

            plot(data[,dim[[classification$class]]])

            cols = dim[[classification$class]]

            model <- lm(data$y~data$x)
            abline(model,col="red")

            bestDistance = 0
            a = c(1,model$coefficients[2])
            b = c(median(data$x),median(data$y))
            t = seq(-100,100,1)

            bestA <- a
            bestB <- b
            bestHits <- Inf
            t = seq(-100,100,1)
            sampleSet <- sample(1:nrow(data),1000)
            governor = 1
            for (step in 1:20) {

                a <- bestA
                b <- bestB

                # Modify the line.
                if (runif(1) <=0.5) {
                    a[1] <- a[1] + runif(1,-0.1,0.1) * governor
                } else {
                    a[2] <- a[2] + runif(1,-0.1,0.1) * governor
                }

                # Calculate distance from line
                distances <- c()
                for (index in sampleSet) {
                    distances <- append(distances, distanceFromLine(a, b, data[index,1:2]))
                }
                hits = which(distances < 0.75)
                lines(a[1]*t+b[1], a[2]*t+b[2], col="red")
                points(b[1],b[2],col="red",pch=19)
                print(paste(step,length(hits)))

                if (length(hits) < bestHits) {
                    bestA <- a
                    governor <- (100-step)/100
                    bestHits <- length(hits)
                }
            }

            a <- bestA

            angle = atan(a[2]/a[1])
            intercept = b[2] + a[2]*(-b[1]/a[1])
            rotData = matrix(0,nrow=nrow(data),ncol=3)
            abline(b=a[2]/a[1],a=intercept,col="red")
            # Attempt to rotate the points about an axis.
            for (row in 1:nrow(data)) {

                rotData[row,] <- c(data[row,1] * cos(-angle) - data[row,2] * sin(-angle) ,
                                   data[row,1] * sin(-angle) + data[row,2] * cos(-angle),data[row,3])
            }
            plot(rotData[,1:2])
            plot(coordinates)
            coordinates <- rotData[,dim[[classification$class]]]
            plot(coordinates)
            par(mfrow=c(1,1))
        }

        fit <- kmeans(coordinates, centers=nClust)
        data <- cbind(data,cluster = fit$cluster)

    }
    else if (classification$class == 5) {
        # If Class 2, examine rotation about the Z axis. Correct if necessary

        fit <- kmeans(coordinates, centers=nClust)
        data <- cbind(data,cluster = fit$cluster)
    }

    #hcl <- hclust(dist(coordinates))
    #clusters = cutree(hcl,nClust)
    #fit <- clara(coordinates,nClust,metric="manhattan")

    #fit.DR <- MclustDR(fit)


    #data <- cbind(data,cluster = clusters)

    plot(data[,1:3],col=data$cluster)

    # Find the centers of each cluster
    centers = matrix(0,nrow = nClust, ncol=3)
    for (cluster in 1:nClust) {
        clustData = data[data$cluster == cluster,]
        centers[cluster,] <-kmeans(clustData[,1:3],1)$centers
    }

    data["clusterDist"] <- NA

    for (index in 1:length(data$cluster)) {
        cluster = data$cluster[index]

        centerX = centers[cluster,1]
        centerY = centers[cluster,2]
        centerZ = centers[cluster,3]

        #print(paste(centerX, centerY, centerZ))

        dist = sqrt((data$x[index] - centers[cluster,1])^2 / sd(data$x)^2 +
                        (data$y[index] - centers[cluster,2])^2 / sd(data$y)^2+
                        (data$z[index] - centers[cluster,3])^2/ sd(data$z)^2)

        data$clusterDist[index] = dist
    }
    data <- data[order(data$clusterDist),]

    # Trim chromosomes that are too far away
    #data <- data[data$clusterDist < median(data$clusterDist) + sd(data$clusterDist)*2.5,]

    colFunc <- colorRampPalette(c("blue","orange","red"))
    colors <- rainbow(max(data$clusterDist)*10)[data$clusterDist*10]
    #col=rainbow(100)[colors],pch=19,
    plot(data[,1:2],cex=0.5,col=colors,pch=19)


    # # Plot the cross sections of the model
    # nPlots = length(seq(1,max(data$z),0.5))
    # mfrow=c(ceiling(sqrt(nPlots)),ceiling(sqrt(nPlots)))
    # par(mfrow=mfrow)
    # for (zindex in seq(1,max(data$z),0.5)) {
    #   plot(data[data$z == zindex,1:2],cex=0.5,col=colors[data$z == zindex],pch=19)
    # }
    # par(mfrow=c(1,1))

    # Load previous densities if possible

    if (file.exists("densities.csv")){
        densities <- read.csv("densities.csv")
    } else {
        densities <- data.frame()
    }

    # Calculate the density for the chromosomes

    for (cluster in 1:max(data$cluster)) {
        hull = convhulln(data[data$cluster==cluster,1:3], options = "FA")
        volume = hull$vol

        #Tm <- delaunayn(data[data$cluster==cluster,1:3])
        #tris <- surf.tri(data[data$cluster==cluster,1:3],Tm)
        #tetramesh(Tm, data[data$cluster==cluster,1:3], alpha=0.5)

        nChrome = length(which(data$cluster == cluster))
        density = nChrome / volume
        cat(paste("Cluster",cluster,":\n"))
        cat(paste("    Chromosomes:",nChrome,"\n"))
        cat(paste("    Volumes:",volume,"um^3\n"))
        cat(paste("    Density:",density,"chromosome/um^3\n"))
        cat(paste("           :",1/density,"um^3 per genome\n"))
        densities <- rbind(densities,data.frame(analysis=directory,
                                                chromosomes=nChrome,
                                                volume=volume,
                                                density=density))
    }
    write.csv(densities,"densities.csv",row.names=F)
}
