setwd("./")

library(fpc)
library(cluster)
 
# CV <- read.csv("Windows_concatenate_PR.dat")
# CV <- read.csv("prova2.dat")
 CV <- read.csv("CVS_clustering_for_kmeans.dat")
 
  k=CENTERS #number of centroids

start <-matrix(c(