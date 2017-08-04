library(ggplot2)
library(plot3D)
library(stats)

coordinates <- read.csv("data/points/a/A_1.31_5.11.csv")
colnames(coordinates) <- c("x", "y", "z")

distance <- function(row1, row2) {
    actual_distance = 0
    for (index in 1:3) {
        actual_distance = actual_distance + (row1[index] - row2[index])^2
    }
    return(sqrt(actual_distance))
}

mean_distances <- c()
for (index1 in 999:nrow(coordinates)) {
    distances = c()
    for (index2 in 1:nrow(coordinates)) {
        if (index1 != index2) {
            distances <- append(distances, distance(coordinates[index1,], coordinates[index2,]))
        }
    }
    hist(sort(distances)[1:10])
    break
}

distances <- unlist(lapply(seq_len(nrow(coordinates)), function(i) median(sqrt(colSums((coordinates[i, ] - t(coordinates))^2)))))

coordinates$avg_distance <- distances

ggplot(coordinates, aes(x, y, col=avg_distance)) +
    geom_point()+
    scale_colour_gradientn(colours=rainbow(4))

ggplot(coordinates, aes(avg_distance)) +
    geom_density()

#scatter3D(coordinates$x, coordinates$y, coordinates$z, colvar = coordinates$avg_distance)

#dist(coordinates)
